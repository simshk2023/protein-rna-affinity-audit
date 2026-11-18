"""Recalculate selected manuscript results from bundled aggregate statistics."""
import argparse
import csv
import json
from pathlib import Path
from statistics import fmean

import numpy as np

from check_v5_tables import REFERENCES, check_comparisons, close
from plot_external import NAMES, make_figure

METHODS = ("copra", "baseline")
ERRORS = ("absolute_MAE", "ordinary_MAE", "cMAE")


def external_summary(reference, seed_rows):
    groups = {(family, method): {} for family in NAMES for method in METHODS}
    for row in seed_rows:
        block = groups[row["family"], row["method"]]
        if row["seed"] in block:
            raise ValueError("Duplicate external family, method and seed")
        block[row["seed"]] = row
    seed_ids = {20260916, 20260917, 20260918}
    if any(set(block) != seed_ids for block in groups.values()):
        raise ValueError("Missing or unexpected external seeds")

    family_metrics = reference["family_metrics"]
    for (family, method), block in groups.items():
        means = {metric: fmean(row[metric] for row in block.values()) for metric in ERRORS}
        close(list(family_metrics[family][method][m] for m in ERRORS), list(means.values()))
        family_metrics[family][method].update(means)

    macro = {
        metric: {m: fmean(family_metrics[f][m][metric] for f in NAMES) for m in METHODS}
        for metric in reference["macro"]
    }
    for metric, methods in macro.items():
        close([reference["macro"][metric][m] for m in METHODS], list(methods.values()))

    draws = np.random.default_rng(20261109).integers(0, len(NAMES), size=(10000, len(NAMES)))
    contrasts = {}
    for metric, expected in reference["contrasts"].items():
        values = np.array([
            family_metrics[f]["baseline"][metric] - family_metrics[f]["copra"][metric]
            for f in NAMES
        ])
        mean = float(values.mean())
        bounds = np.quantile(values[draws].mean(axis=1), [.025, .975]).tolist()
        close([expected["mean"], *expected["bootstrap95"]], [mean, *bounds])
        contrasts[metric] = dict(
            mean=mean, bootstrap95=bounds, family_values=dict(zip(NAMES, values.tolist()))
        )
    return dict(coverage=reference["coverage"], family_metrics=family_metrics,
                macro=macro, contrasts=contrasts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    output = parser.parse_args().output
    output.mkdir(parents=True, exist_ok=True)

    table = check_comparisons(REFERENCES / "matched_comparison_aggregates.json", output)
    with (output / "table3.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)

    external = json.loads((REFERENCES / "external_aggregates.json").read_text(encoding="utf-8"))
    seeds = json.loads((REFERENCES / "external_seed_error_aggregates.json").read_text(encoding="utf-8"))
    results = external_summary(external, seeds)
    (output / "external_S13.json").write_text(json.dumps(results, indent=2, allow_nan=False), encoding="utf-8")
    make_figure(results, seeds, output / "figureE1.pdf")
    report = dict(status="PASS", input_grain="cohort and external family/seed aggregates",
                  table3_configurations=len(table), table_s37_rows=37,
                  external_seed_error_means_checked=len(NAMES) * len(METHODS) * len(ERRORS),
                  external_macro_cells_checked=sum(len(x) for x in results["macro"].values()),
                  external_bootstrap_contrasts=len(results["contrasts"]),
                  external_quantitative_pairs=sum(x["point_pairs"] for x in results["coverage"].values()),
                  figure_e1_rendered=True)
    (output / "aggregate_verification.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
