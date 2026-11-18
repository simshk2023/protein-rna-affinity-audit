"""Replay comparisons and mechanism diagnostics from study outputs and a local provider export."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
import tempfile
from collections import defaultdict

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "study_outputs/comparison"


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def frozen_input(path, identity):
    if not path.is_file():
        raise ValueError(f"Missing {identity['path']}; see data_access.md")
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    if path.stat().st_size != identity["bytes"] or digest != identity["sha256"]:
        raise ValueError(f"Wrong frozen input for {identity['path']}; see data_access.md")


def materialize(provider, directory):
    metadata = read(DATA / "membership.json")
    with provider.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    labels = defaultdict(list)
    for record in metadata["source_links"]:
        source = rows[int(record["source_row_number"]) - 1]
        if source["Entry_id"] != record["entry_id"]:
            raise ValueError("ProNAB snapshot differs from the study's source-row/Entry-ID mapping")
        value = float(source[record["raw_column"]].strip())
        if not math.isfinite(value) or value <= 0:
            raise ValueError("Invalid retained provider affinity")
        labels[record["unit_id"]].append(-math.log10(value))
    units = {
        key: dict(family_id=row["family_id"], mutation_tokens=row["mutation_tokens"],
                  PDB="1JBR" if row["position_control"] else "", label=statistics.median(labels[key]))
        for key, row in metadata["units"].items()
    }
    pairs = {
        key: dict(row, delta_pkd=units[row["mutant_unit_id"]]["label"] - units[row["wt_unit_id"]]["label"])
        for key, row in metadata["pairs"].items()
    }
    record = dict(D=units, P=pairs, U=read(DATA / "endpoint_predictions.json"),
                  ex=metadata["censor_exclusions"], alias_exclusions=metadata["alias_exclusions"])
    (directory / "saved_inputs.json").write_text(json.dumps(record, allow_nan=False), encoding="utf-8")
    for name in ("pra_inputs.json", "m2_predictions.json", "m2_cohort_manifest.json"):
        shutil.copyfile(DATA / name, directory / name)


def run(script, *args):
    subprocess.run([sys.executable, "-B", str(ROOT / script), *map(str, args)], check=True)


def three_arm(output):
    ref = read(DATA / "three_arm_figure7_reference.json")
    cells = [cell for scope in ref["matrix"].values() for metric in scope.values() for cell in metric.values()]
    for cell in cells:
        values = cell["seeds"]
        if (abs(statistics.fmean(values) - cell["mean"]) > 1e-10
                or abs(statistics.stdev(values) - cell["sd"]) > 1e-10):
            raise ValueError("Archived three-arm seed summary differs from Figure 7 reference")
    (output / "three_arm_seed_aggregate_check.json").write_text(
        json.dumps(dict(status="PASS", cells=len(cells), scope="published seed summaries")), encoding="utf-8"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pronab", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--services", type=Path, help="Matching authorized historical response directory")
    parser.add_argument("--mechanisms", action="store_true")
    parser.add_argument("--pra-formula", type=Path, help="Matching independently acquired PRA-Pred formula JSON")
    args = parser.parse_args()
    if args.pra_formula and not args.pra_formula.is_file():
        parser.error("Acquire the PRA-Pred coefficients and format them as described in data_access.md.")
    out = args.output.resolve()
    if out.is_relative_to(ROOT):
        parser.error("--output must be outside the package")
    frozen = read(ROOT / "input_manifest.json")
    try:
        frozen_input(args.pronab, next(r for r in frozen if r["path"] == "pronab.csv"))
        if args.services:
            for identity in frozen:
                if identity["path"].startswith("services/"):
                    frozen_input(args.services / identity["path"].removeprefix("services/"), identity)
    except ValueError as error:
        parser.error(str(error))
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="copra_compare_") as tmp:
        work = Path(tmp)
        materialize(args.pronab, work)
        service_args = []
        if args.services:
            run("replay_services.py", "--inputs", args.services.resolve(), "--output", work / "services")
            service_args = ["--service-predictions", work / "services/service_pair_predictions.json"]
        run("replay_methods.py", "--inputs", work, "--output", out, "--subset", *service_args)
        if args.pra_formula:
            run("replay_pra_mechanism.py", "--inputs", DATA, "--output", out / "pra",
                "--subset", "--formula", args.pra_formula.resolve())
        if args.mechanisms:
            run("replay_mechanisms.py", "--inputs", DATA, "--methods", work, "--output", out / "m2", "--m2-subset")
            three_arm(out)
    methods = read(out / "method_results.json")
    print(json.dumps(dict(status="PASS", scope="release prediction subset", compared_scopes=list(methods),
                          prorb_replay=False, historical_services=args.services is not None,
                          mechanisms=args.mechanisms, pra_formula_check=args.pra_formula is not None, model_execution_count=0)))


if __name__ == "__main__":
    main()
