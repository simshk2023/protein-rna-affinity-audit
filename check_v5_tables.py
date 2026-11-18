"""Check the v5 display tables from the archived-result numerical replay."""
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parent
REFERENCES = ROOT / "manuscript_v5"
NUMERIC = re.compile(r"[-+]?\d+(?:\.\d+)?")


def numbers(text):
    return [float(value) for value in NUMERIC.findall(text.replace("$-$", "-"))]


def rows(name):
    for line in (REFERENCES / f"{name}.tex").read_text(encoding="utf-8").splitlines():
        if "&" in line and line.endswith(r"\\") and not line.startswith((r"\textbf{", "Measure &")):
            yield [cell.strip() for cell in line.removesuffix(r"\\").split("&")]


def close(displayed, actual):
    if len(displayed) != len(actual) or any(
        abs(expected - value) > 1e-10 for expected, value in zip(displayed, actual)
    ):
        raise ValueError(f"Manuscript mismatch: displayed {displayed}, reproduced {actual}")


def check_primary(replay):
    legacy = json.loads((replay / "checks/table2_keyed_check.json").read_text())
    data = {row["key"].split("/copra")[0]: row["expected"] for row in legacy}
    for name, keys in (
        ("table2", ["absolute/mae", "ordinary/ordinary_mae", "local/ordinary_mae",
                    "local/cmae", "local/position_mean_spearman", "local/sd_ratio",
                    "local/amplitude_distance_from_one"]),
        ("table_s36", ["ordinary/linked_absolute_mae", "local/linked_absolute_mae"]),
    ):
        displayed = [numbers(" ".join(cells[1:])) for cells in rows(name)]
        if len(displayed) != len(keys):
            raise ValueError(f"{name} has {len(displayed)} rows, expected {len(keys)}")
        for values, key in zip(displayed, keys):
            close(values, data[key])

    summary = json.loads((replay / "primary/replay_summary.json").read_text())
    local = summary["scopes"]["local"]
    boot = local["bootstrap"]
    families = boot["family_order"]
    draws = np.random.default_rng(boot["seed"]).integers(
        0, len(families), size=(boot["draws"], len(families))
    )
    differences = []
    for metric in ("ordinary_mae", "cmae"):
        delta = np.array([
            local["methods"]["zero"]["family_values"][family][metric]
            - local["methods"]["copra"]["family_values"][family][metric]
            for family in families
        ])
        differences.extend([delta.mean(), *np.quantile(delta[draws].mean(axis=1), [.025, .975])])
    reference = (REFERENCES / "table_s36.tex").read_text(encoding="utf-8").split("contrasts were ")[1]
    close(numbers(reference.split("See Tables")[0]), [round(x, 4) for x in differences])


def check_comparisons(replay, output):
    methods_path = replay / "methods/method_results.json"
    raw = json.loads(methods_path.read_text())
    names = ["deeppnap", "omni_official", "omni_controlled", "prdelta", "mcsm",
             "gbnl", "svsbi", "prorb", "pra", "m2_aligned"]
    display = list(rows("table3"))
    if len(display) != len(names):
        raise ValueError("Unexpected main Table 3 configuration count")
    for cells, method in zip(display, names):
        _, scope, comparator, *reported = cells
        block = raw[scope]
        reference_key = {"Sequence": "seq_baseline", "CoPRA": "baseline"}[comparator]
        reference, candidate = [block["methods"][key] for key in (reference_key, method)]
        coverage = [block.get(k) for k in ("units", "pairs", "bshape_pairs")]
        numerics = [
            None if key not in candidate else round(reference[key] - candidate[key], 4)
            for key in ("absolute", "ordinary", "cmae")
        ]
        for printed, actual in zip(reported, coverage + numerics):
            if actual is None:
                if printed != "---":
                    raise ValueError(f"Unexpected reported absolute metric for {method}")
            else:
                close(numbers(printed), [actual])

    subprocess.run([
        sys.executable, "-B", str(ROOT / "check_table3.py"), str(methods_path),
        "--manuscript", str(REFERENCES / "table_s37.tex"),
        "--output", str(output / "table_s37_check.json"),
    ], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("replay", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    check_primary(args.replay)
    check_comparisons(args.replay, args.output)
    result = {"status": "PASS", "table_2_rows": 7, "table_3_configurations": 10,
              "table_s36_linked_rows": 2, "table_s36_zero_contrasts": 2,
              "table_s37_rows": 37}
    (args.output / "manuscript_tables.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result))


if __name__ == "__main__":
    main()

