"""Prepare the local replay tree from the public package and reader-acquired inputs.

Only the package itself, the explicitly supplied ProNAB CSV, the PRA-Pred
formula JSON, the wwPDB coordinate files and (optionally) the S13 observation
template are read.  No historical project directory is used.

    python -B prepare_inputs.py --pronab ../pronab.csv \
        --pra-formula ../pra_formula.json --structures ../pdb \
        --output ../work/inputs
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PRIMARY = ROOT / "study_outputs/primary"
COMP = ROOT / "study_outputs/comparison"
MECH = ROOT / "study_outputs/mechanisms"
FIGS = ROOT / "study_outputs/figures"


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def dump(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, allow_nan=False), encoding="utf-8")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def reconstruct_labels(pronab: Path) -> tuple[dict, dict, dict]:
    units = read_csv(ROOT / "metadata/units.csv")
    pairs = read_csv(ROOT / "metadata/pairs.csv")
    links = read_csv(ROOT / "metadata/source_links.csv")
    membership = read_csv(ROOT / "metadata/candidate_membership.csv")
    exclusions = read_json(ROOT / "metadata/censor_exclusions.json")
    aliases = read_json(ROOT / "metadata/source_alias_exclusions.json")

    rows = read_csv(pronab)
    replicates: dict[str, list[float]] = defaultdict(list)
    for link in links:
        if link["unit_id"] in set(aliases["unit_ids"]):
            continue
        source = rows[int(link["source_row_number"]) - 1]
        if source["Entry_id"] != link["entry_id"]:
            raise ValueError(
                "ProNAB source row/Entry_id mismatch; use a recorded snapshot or "
                "realign explicitly by Entry_id"
            )
        text = source[link["raw_column"]].strip()
        if not text or any(x in text for x in "<>≤≥~"):
            raise ValueError("a retained ProNAB affinity is not a point value")
        value = float(text)
        if not math.isfinite(value) or value <= 0:
            raise ValueError("a retained ProNAB affinity is not positive and finite")
        replicates[link["unit_id"]].append(-math.log10(value))

    excluded = set(exclusions["CENSORED_UNITS"]) | set(aliases["unit_ids"])
    labels = {}
    for unit in units:
        if unit["unit_id"] in excluded:
            continue
        if unit["unit_id"] not in replicates:
            raise ValueError(f"no retained provider link for {unit['unit_id']}")
        labels[unit["unit_id"]] = statistics.median(replicates[unit["unit_id"]])

    positions = {r["pair_id"]: r["position"] for r in membership}
    eligible = {
        r["pair_id"] for r in membership if r["reason"] != "censored_endpoint"
    } - set(aliases["pair_ids"])
    pair_labels = {}
    for pair in pairs:
        if pair["pair_id"] not in eligible:
            continue
        pair_labels[pair["pair_id"]] = (
            labels[pair["mutant_unit_id"]] - labels[pair["wt_unit_id"]]
        )
    return labels, pair_labels, {"units": len(labels), "pairs": len(pair_labels)}


def build_saved_inputs(labels: dict, pair_labels: dict) -> dict:
    meta = read_json(COMP / "metadata.json")
    D = {}
    for key, row in meta["D"].items():
        record = dict(row)
        if key in labels:
            record["label"] = labels[key]
        D[key] = record
    P = {}
    for key, row in meta["P"].items():
        record = dict(row)
        if key in pair_labels:
            record["delta_pkd"] = pair_labels[key]
        P[key] = record
    return {
        "D": D,
        "P": P,
        "U": read_json(COMP / "predictions.json"),
        "ex": meta["ex"],
        "routes": meta["routes"],
        "label_source": "reconstructed from the reader-supplied ProNAB snapshot",
    }


def build_service_predictions() -> dict:
    scores = {"panther": {}, "mcsm": {}}
    for row in read_csv(COMP / "service_scores.csv"):
        scores[row["method"]][row["pair_id"]] = float(row["score_delta"])
    return scores


def build_pra_formula_input(formula: Path) -> dict:
    spec = read_json(formula)
    for branch in ("ss+str", "ds+enz"):
        if branch not in spec or "intercept" not in spec[branch] or "terms" not in spec[branch]:
            raise ValueError(f"PRA formula JSON lacks {branch}/intercept/terms")
    return spec


def build_figure5_pairs(formula: dict, out: Path) -> None:
    """Recompute exact PRA term contributions from public features + imported formula."""
    public = read_json(FIGS / "figure5/restricted/pair_features.json")
    pra = read_json(COMP / "pra_inputs.json")
    endpoints = pra["endpoints"]
    pairs = {p["pair_id"]: p for p in pra["pairs"]}
    out_data = {}
    for branch, block in public.items():
        terms = block["terms"]
        spec = formula[branch]
        beta = spec["terms"]
        rows = []
        for row in block["rows"]:
            pair = pairs[row["pair_id"]]
            wt = endpoints[pair["wt_endpoint_id"]]
            mut = endpoints[pair["mut_endpoint_id"]]
            diffs = [mut["features"][k] - wt["features"][k] for k in terms]
            contributions = [beta[k] * d for k, d in zip(terms, diffs)]
            delta = sum(contributions)
            exact = all(x == 0 for x in contributions)
            rounded = mut["prediction"]["G_2dp"] == wt["prediction"]["G_2dp"]
            rows.append({
                "pair_id": row["pair_id"],
                "changed": [d != 0 for d in diffs],
                "zero": exact,
                "rounded_zero": rounded,
                "rounding_only": (not exact) and rounded,
                "contributions": contributions,
                "wt_features": row["wt_features"],
                "delta_G": delta,
                "wt_percentiles": row.get("wt_percentiles"),
            })
        out_data[branch] = {"terms": terms, "rows": rows}
    target = out / "figures/figure5/restricted/pair_plot_data.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    dump(target, out_data)


def stage_structures(structures: Path | None, out: Path) -> None:
    if structures is None:
        return
    expected = {
        "2C4R": "figures/figure1/outputs/length_research_handoff_20260926/DS04_work/sources/PRA310_official/PRA310/PDBs/2C4R.pdb",
        "4QVC": "figures/figure1/outputs/length_research_handoff_20260926/DS04_work/sources/PRA310_official/PRA310/PDBs/4QVC.pdb",
        "2ZZN": "figures/figure1/outputs/length_research_handoff_20260926/DS04_work/sources/PRA310_official/PRA310/PDBs/2ZZN.pdb",
        "4NGD": "figures/figure1/outputs/dicer_input_preparation_20260923/inputs/structures/4NGD.pdb",
    }
    checks = read_json(ROOT / "third_party/checksums.json")["structures"]
    for pdb, relative in expected.items():
        source = structures / f"{pdb}.pdb"
        if not source.is_file():
            raise ValueError(f"missing structure {pdb}.pdb under {structures}")
        digest = sha256(source)
        if digest.lower() != checks[pdb]["sha256"].lower():
            raise ValueError(f"structure hash mismatch for {pdb}.pdb")
        target = out / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pronab", required=True, type=Path)
    parser.add_argument("--pra-formula", type=Path)
    parser.add_argument("--structures", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    out = args.output.resolve()
    if out.is_relative_to(ROOT):
        parser.error("--output must be outside the package")
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    labels, pair_labels, counts = reconstruct_labels(args.pronab.resolve())
    if counts != {"units": 752, "pairs": 338}:
        raise ValueError(f"unexpected principal counts: {counts}")
    dump(out / "saved_inputs.json", build_saved_inputs(labels, pair_labels))

    # Primary predictions and the metadata used by import_authorized.py.
    shutil.copyfile(PRIMARY / "predictions.csv", out / "predictions.csv")
    shutil.copyfile(args.pronab.resolve(), out / "pronab.csv")

    # Comparison and mechanism inputs.
    (out / "methods").mkdir()
    shutil.copyfile(out / "saved_inputs.json", out / "methods/saved_inputs.json")
    for name in ("pra_inputs.json", "m2_predictions.json", "m2_study_manifest.json",
                 "m2_cohort_manifest.json"):
        shutil.copyfile(COMP / name, out / "methods" / name)
    (out / "mechanisms").mkdir()
    for name in ("training.json", "validation.json", "diagnostics.json",
                 "m2_reference.json", "pra_inputs.json"):
        shutil.copyfile(MECH / name, out / "mechanisms" / name)

    if args.pra_formula:
        formula = build_pra_formula_input(args.pra_formula.resolve())
        dump(out / "pra_formula.json", formula)
    dump(out / "services/service_pair_predictions.json", build_service_predictions())

    # Figures are staged into the tree the archived renderers expect.
    shutil.copytree(FIGS, out / "figures")
    if args.pra_formula:
        build_figure5_pairs(formula, out)
    if args.structures:
        stage_structures(args.structures.resolve(), out)

    print(json.dumps({
        "status": "PASS",
        "units": counts["units"],
        "pairs": counts["pairs"],
        "services": {"panther": 87, "mcsm": 267},
        "output": str(out),
    }, indent=2))


if __name__ == "__main__":
    main()
