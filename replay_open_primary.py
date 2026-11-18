"""Score the released primary subset using a reader-obtained ProNAB snapshot."""

import argparse
import csv
import hashlib
import json
import math
import shutil
import statistics
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STUDY = ROOT / "study_outputs/primary"
SEEDS = ["20260916", "20260917", "20260918"]
sys.path.insert(0, str(ROOT / "toolkit"))
import evaluate


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def replay(pronab, output):
    identity = json.loads((STUDY / "study_identity.json").read_text(encoding="utf-8"))
    provider = identity["pronab"]
    if pronab.stat().st_size != provider["bytes"] or sha256(pronab) != provider["sha256"]:
        raise ValueError("ProNAB input differs from the study snapshot. See study_outputs/primary/ACCESS.md.")
    for file, digest in identity["files"].items():
        if sha256(STUDY / file) != digest:
            raise ValueError(f"Primary input SHA-256 mismatch: {file}")

    source_rows = read_csv(pronab)
    links = read_csv(STUDY / "source_links.csv")
    replicates = defaultdict(list)
    for link in links:
        source = source_rows[int(link["source_row_number"]) - 1]
        if source["Entry_id"] != link["entry_id"]:
            raise ValueError("ProNAB Entry ID does not match its source row")
        try:
            value = float(source[link["raw_column"]])
        except ValueError:
            raise ValueError("A retained ProNAB Kd is not a point value") from None
        if not math.isfinite(value) or value <= 0:
            raise ValueError("A retained ProNAB Kd is not positive and finite")
        replicates[link["unit_id"]].append(-math.log10(value))
    units = [
        dict(unit_id=r["unit_id"], family=r["family_id"], observed=statistics.median(replicates[r["unit_id"]]))
        for r in read_csv(STUDY / "units.csv")
    ]
    observed = {r["unit_id"]: r["observed"] for r in units}
    pairs = [
        dict(pair_id=r["pair_id"], family=r["family_id"], landscape=r["landscape_id"],
             position=r["position"], wt_unit_id=r["wt_unit_id"], mutant_unit_id=r["mutant_unit_id"],
             observed_delta=observed[r["mutant_unit_id"]] - observed[r["wt_unit_id"]])
        for r in read_csv(STUDY / "pairs.csv")
    ]
    methods = {
        name: dict(type="endpoint", seeds=SEEDS, identity="frozen controlled family OOF")
        for name in ("baseline", "copra")
    }
    methods["zero"] = dict(type="zero", seeds=["fixed"], identity="zero-change reference")
    config = dict(units="units.csv", pairs="pairs.csv", predictions="predictions.csv", unit="pKd",
                  direction="mutant minus WT; higher is stronger", reference="baseline", methods=methods)
    with tempfile.TemporaryDirectory(prefix="copra_primary_") as tmp:
        work = Path(tmp)
        write_csv(work / "units.csv", units)
        write_csv(work / "pairs.csv", pairs)
        shutil.copyfile(STUDY / "predictions.csv", work / "predictions.csv")
        (work / "config.json").write_text(json.dumps(config), encoding="utf-8")
        results = evaluate.run(work / "config.json")
    counts = dict(
        units=results["absolute"]["common_units"],
        pairs=results["coverage"]["common_pairs"],
        source_links=len(links),
    )
    if any(identity[key] != value for key, value in counts.items()):
        raise ValueError("Primary replay coverage differs from the released subset")
    summary = dict(status="PASS", scope=identity["scope"], **counts,
                   local_pairs=results["coverage"]["local_pairs"], model_execution_count=0)
    output.mkdir(parents=True, exist_ok=True)
    results["replay"] = summary
    (output / "primary_subset_replay.json").write_text(
        json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary))


def main():
    parser = argparse.ArgumentParser(description=__doc__, epilog='Obtain the frozen ProNAB export using study_outputs/primary/ACCESS.md and pass --pronab PATH.')
    parser.add_argument("--pronab", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if not args.pronab.is_file():
        parser.error('Obtain the frozen ProNAB export using study_outputs/primary/ACCESS.md and pass --pronab PATH.')
    output = args.output.resolve()
    if output.is_relative_to(ROOT):
        parser.error("--output must be outside the release source directory")
    replay(args.pronab.resolve(), output)


if __name__ == "__main__":
    main()
