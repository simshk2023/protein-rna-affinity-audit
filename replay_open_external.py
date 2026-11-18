"""Replay S13 from the released author predictions and reader-obtained observations."""

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

from replay_external import NAMES, SEEDS, read_csv, replay_family_inputs

ROOT = Path(__file__).resolve().parent
STUDY = ROOT / "study_outputs/external"
NM_PER_UNIT = {"M": 1e9, "mM": 1e6, "uM": 1e3, "µM": 1e3, "μM": 1e3, "nM": 1.0, "pM": 1e-3}


def read_observations(path, expected_ids):
    columns = ("unit_id", "status", "kd_value", "kd_unit")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != columns:
            raise ValueError("Observation columns must match observations_template.csv")
        rows = list(reader)
    if len(rows) != len(expected_ids) or {r["unit_id"] for r in rows} != expected_ids:
        raise ValueError("Observations must contain exactly one row for each source member")

    pkd = {}
    for row in rows:
        member = row["unit_id"]
        if None in row or any(value is None for value in row.values()):
            raise ValueError(f"Invalid CSV row for {member}")
        status, value, unit = (row[key].strip() for key in columns[1:])
        if status == "censored_or_nonbinding" and not value and not unit:
            pkd[member] = None
        elif status == "point" and unit in NM_PER_UNIT:
            try:
                kd_nm = float(value) * NM_PER_UNIT[unit]
            except ValueError:
                raise ValueError(f"Invalid point Kd for {member}") from None
            if not math.isfinite(kd_nm) or kd_nm <= 0:
                raise ValueError(f"Point Kd must be positive and finite for {member}")
            pkd[member] = 9 - math.log10(kd_nm)
        else:
            raise ValueError(f"Complete the observation for {member} using the instructions in ACCESS.md")
    return pkd


def replay(observations, output):
    identity = json.loads((STUDY / "study_identity.json").read_text(encoding="utf-8"))
    for file, digest in identity["files"].items():
        if hashlib.sha256((STUDY / file).read_bytes()).hexdigest() != digest:
            raise ValueError(f"External input SHA-256 mismatch: {file}")

    members = read_csv(STUDY / "source_members.csv")
    observations_by_unit = read_observations(observations, {r["unit_id"] for r in members})
    pairs = read_csv(STUDY / "pairs.csv")
    predictions = defaultdict(lambda: defaultdict(dict))
    for row in read_csv(STUDY / "endpoint_predictions.csv"):
        family = row["family"]
        task = (row["method"], int(row["seed"]), int(row["fold"]))
        member = row["unit_id"]
        prediction = float(row["prediction"])
        task_predictions = predictions[family][task]
        if member in task_predictions or not math.isfinite(prediction):
            raise ValueError("Duplicate or nonfinite external endpoint prediction")
        task_predictions[member] = prediction

    if set(predictions) != set(NAMES):
        raise ValueError("External prediction families differ from the S13 cohort")
    expected_tasks = {
        (method, seed, fold)
        for method in ("copra", "baseline") for seed in SEEDS for fold in range(5)
    }
    data = {}
    for name in NAMES:
        labels = {
            r["unit_id"]: dict(token=r["source_token"], pkd=observations_by_unit[r["unit_id"]])
            for r in members if r["family"] == name
        }
        family_pairs = [r for r in pairs if r["family"] == name]
        tasks = predictions[name]
        if set(tasks) != expected_tasks or any(set(rows) != set(labels) for rows in tasks.values()):
            raise ValueError(f"Incomplete external task/member coverage for {name}")
        data[name] = labels, family_pairs, tasks

    replay_family_inputs(data, output, verified_files=len(identity["files"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__, epilog='Copy study_outputs/external/observations_template.csv, obtain observations using ACCESS.md, and pass --observations PATH.')
    parser.add_argument("--observations", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if not args.observations.is_file():
        parser.error('Copy study_outputs/external/observations_template.csv, obtain observations using ACCESS.md, and pass --observations PATH.')
    output = args.output.resolve()
    if output.is_relative_to(ROOT):
        parser.error("--output must be outside the release source directory")
    replay(args.observations.resolve(), output)


if __name__ == "__main__":
    main()
