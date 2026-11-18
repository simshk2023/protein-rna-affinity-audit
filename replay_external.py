"""Independent numerical review of the three new source-defined evaluations.

Reads archived endpoint predictions; never executes a model or changes source data.
Uses standard-library sums, centered residuals and tied ranks.
"""
import argparse
import csv
import hashlib
import json
import math
import re
import statistics
from pathlib import Path

import numpy as np

from plot_external import make_figure

REFERENCE = Path(__file__).resolve().parent / 'manuscript_v5'
SEEDS = (20260916, 20260917, 20260918)
NAMES = ('sf3a1', 'sepsecs', 'npl3')


def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))


def avg(values):
    return statistics.fmean(values)


def ranks(values):
    ordered = sorted(range(len(values)), key=values.__getitem__)
    out = [0.0] * len(values)
    i = 0
    while i < len(ordered):
        j = i + 1
        while j < len(ordered) and values[ordered[j]] == values[ordered[i]]:
            j += 1
        for index in ordered[i:j]:
            out[index] = (i + j - 1) / 2
        i = j
    return out


def spearman(a, b):
    a, b = ranks(a), ranks(b)
    ac = [x - avg(a) for x in a]
    bc = [x - avg(b) for x in b]
    denominator = math.sqrt(math.fsum(x*x for x in ac) * math.fsum(x*x for x in bc))
    return math.fsum(x*y for x, y in zip(ac, bc)) / denominator if denominator else None


def load_archived_family(inputs, name):
    directory = inputs / 'restricted' / (name + '_tranche')
    label_rows = json.loads((directory / 'labels.json').read_text())
    labels = {r['unit_id']: r for r in label_rows}
    assert len(labels) == len(label_rows), (name, 'duplicate units')
    all_pairs = read_csv(directory / 'external_pairs.csv')
    files = list((directory / 'predictions').glob('*/predictions_units.csv'))
    assert len(files) == 30, (name, len(files))
    endpoint_predictions = {}
    for path in files:
        method, fold, seed = re.fullmatch(r'(copra|baseline)_final_f([0-4])_s(\d+)', path.parent.name).groups()
        rows = read_csv(path)
        assert len(rows) == len(labels) and {r['unit_id'] for r in rows} == set(labels)
        assert all(r['task_id'] == path.parent.name and int(r['seed']) == int(seed) for r in rows)
        values = {r['unit_id']: float(r['prediction']) for r in rows}
        assert all(math.isfinite(v) for v in values.values())
        key = (method, int(seed), int(fold))
        assert key not in endpoint_predictions
        endpoint_predictions[key] = values
    return labels, all_pairs, endpoint_predictions


def replay_family_inputs(data, output, verified_files):
    expected = json.loads((REFERENCE / 'external_aggregates.json').read_text(encoding='utf-8'))
    families, coverage, seed_rows = {}, {}, []
    max_difference = 0.0
    cells = 0
    for name in NAMES:
        labels, all_pairs, endpoint_predictions = data[name]
        files_count = len(endpoint_predictions)
        points = {u: r['pkd'] for u, r in labels.items() if r['pkd'] is not None}
        pairs = [r for r in all_pairs
                 if r['mutant_unit_id'] in points and r['wt_unit_id'] in points]
        positions = [int(re.search(r'\d+', labels[r['mutant_unit_id']]['token']).group()) for r in pairs]
        assert len(set(positions)) == len(positions), (name, 'position weighting needs grouping')
        families[name] = {}
        for method in ('copra', 'baseline'):
            metrics = []
            for seed in SEEDS:
                predictions = {u: avg(endpoint_predictions[method, seed, f][u] for f in range(5)) for u in labels}
                truth = [points[r['mutant_unit_id']] - points[r['wt_unit_id']] for r in pairs]
                response = [predictions[r['mutant_unit_id']] - predictions[r['wt_unit_id']] for r in pairs]
                residual = [p-y for p, y in zip(response, truth)]
                row = dict(absolute_MAE=avg(abs(predictions[u]-points[u]) for u in points),
                           ordinary_MAE=avg(abs(x) for x in residual),
                           cMAE=avg(abs(x-avg(residual)) for x in residual),
                           zero_ordinary_MAE=avg(abs(y) for y in truth),
                           zero_cMAE=avg(abs(y-avg(truth)) for y in truth),
                           spearman=spearman(response, truth),
                           SD_ratio=statistics.pstdev(response)/statistics.pstdev(truth),
                           direction=avg(((p > 0)-(p < 0)) == ((y > 0)-(y < 0)) for p, y in zip(response, truth)))
                metrics.append(row)
                seed_rows.append(dict(family=name, method=method, seed=seed, **row))
            combined = {k: avg(row[k] for row in metrics) for k in metrics[0]}
            families[name][method] = combined
            for key, value in combined.items():
                diff = abs(value - expected['family_metrics'][name][method][key])
                assert diff < 1e-10, (name, method, key, diff)
                max_difference = max(max_difference, diff)
                cells += 1
        coverage[name] = dict(
            checkpoints=files_count, input_endpoints=len(labels), point_endpoints=len(points),
            attempted_pairs=len(all_pairs), point_pairs=len(pairs),
            censored_or_nonbinding=len(all_pairs) - len(pairs), local_positions=len(positions)
        )
        if coverage[name] != expected['coverage'][name]:
            raise ValueError(f'External coverage mismatch for {name}')
    macro = {k: {m: avg(families[n][m][k] for n in NAMES) for m in ('copra', 'baseline')}
             for k in families[NAMES[0]]['copra']}
    for key, methods in macro.items():
        for method, value in methods.items():
            diff = abs(value - expected['macro'][key][method])
            assert diff < 1e-10, (key, method, diff)
            max_difference = max(max_difference, diff)
            cells += 1
    draws = np.random.default_rng(20261109).integers(0, len(NAMES), size=(10000, len(NAMES)))
    contrasts = {}
    for key, archived in expected['contrasts'].items():
        values = np.array([families[name]['baseline'][key] - families[name]['copra'][key] for name in NAMES])
        mean = float(values.mean())
        interval = np.quantile(values[draws].mean(axis=1), [.025, .975]).tolist()
        if not np.allclose([mean, *interval], [archived['mean'], *archived['bootstrap95']], rtol=0, atol=1e-10):
            raise ValueError(f'External contrast mismatch for {key}')
        contrasts[key] = dict(mean=mean, bootstrap95=interval, family_values=dict(zip(NAMES, values.tolist())))
    report = dict(status='PASS', families=3, coverage=coverage, cells_checked=cells,
                  maximum_absolute_difference=max_difference,
                  contrasts_checked=len(contrasts), frozen_files_verified=verified_files,
                  family_metrics=families, macro=macro, contrasts=contrasts, model_execution_count=0)
    output.mkdir(parents=True, exist_ok=True)
    (output / 'external_replay.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    make_figure(report, seed_rows, output / 'figureE1.pdf')
    print(json.dumps({k: v for k, v in report.items() if k not in ('family_metrics', 'macro', 'contrasts')}, indent=2))



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    frozen = json.loads((REFERENCE / 'external_input_manifest.json').read_text())
    for relative, digest in frozen.items():
        actual = hashlib.sha256((args.inputs / relative).read_bytes()).hexdigest()
        if actual != digest:
            raise ValueError(f'External input SHA-256 mismatch: {relative}')
    data = {name: load_archived_family(args.inputs, name) for name in NAMES}
    replay_family_inputs(data, args.output, len(frozen))


if __name__ == '__main__':
    main()
