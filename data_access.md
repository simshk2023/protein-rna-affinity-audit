# Reader acquisition and replay

The package runs from author-generated predictions and aggregate inputs. Four
things are obtained separately and passed explicitly; nothing is read from a
historical project directory.

| Input | Obtain | Pass to | Needed for |
|---|---|---|---|
| ProNAB measurement CSV | <https://web.iitm.ac.in/bioinfo2/pronab/>; match by `Entry_id` + `raw_column`; see [ACQUISITION.md](ACQUISITION.md) | `--pronab` | Table 2, S36, all per-record comparisons |
| PRA-Pred fitted equations | <https://github.com/Harinikannan98/PRA-Pred>; JSON schema in `third_party/pra_formula_schema.json` | `--pra-formula` | Figure 5, PRA term path |
| wwPDB coordinates | `python -B third_party/fetch_structures.py --output ../imports/pdb` | `--structures` | Figure 1 |
| S13 observations | Fill `study_outputs/external/observations_template.csv` from the three articles; see `study_outputs/external/ACCESS.md` | `--observations` | Figure E1, S13 |

## Start here

```sh
python -B reproduce.py prepare --pronab ../imports/pronab.csv \
    --pra-formula ../imports/pra_formula.json --structures ../imports/pdb \
    --output ../work
```

`prepare` verifies the package manifest, builds the local label table, stages
the figures, and writes `work/saved_inputs.json`, `work/pra_formula.json`,
`work/services/service_pair_predictions.json` and the `work/figures/` tree.
The ProNAB file is copied into `work/pronab.csv`; outputs stay outside the
package.

## The four routes

```sh
python -B reproduce.py primary --pronab ../imports/pronab.csv --output ../results/primary
python -B reproduce.py full --pronab ../imports/pronab.csv \
    --pra-formula ../imports/pra_formula.json --structures ../imports/pdb \
    --output ../results/full
python -B reproduce.py external --observations ../imports/observations.csv \
    --output ../results/external
python -B reproduce.py all --pronab ../imports/pronab.csv \
    --pra-formula ../imports/pra_formula.json --structures ../imports/pdb \
    --observations ../imports/observations.csv --output ../results/all
```

- **primary**: reconstructs the 752 labels through `metadata/units.csv`,
  `metadata/pairs.csv` and the 827 retained `source_links`; asserts
  752/338/220; recomputes Table 2 (9 rows, 44 cells) and S36.
- **full**: uses the bundled method prediction blocks and service scores to
  recompute Table 3 (10 configurations), S37 (37 rows), PRA (128 endpoints,
  106 main pairs, 42 exact zeros, 13 rounding ties), M2 (36 condition rows, 72
  metrics against the archived reference) and the three-arm summaries (30
  cells, 90 seed values); renders Figures 1–7 and S14.
- **external**: averages the five prepared folds within each family, method and
  seed, scores the three seeds, checks S13's 64 cells and five paired
  contrasts, and renders Figure E1.
- **all**: full + external + `frozen/reference.json` comparison.

## What each route reads

| Route | Package data | Imported |
|---|---|---|
| primary | `study_outputs/primary/predictions.csv`, `metadata/*` | ProNAB |
| full | `study_outputs/comparison/predictions.json`, `metadata.json`, `service_scores.csv`, `study_outputs/mechanisms/*`, `study_outputs/figures/*` | ProNAB, PRA formula, PDB |
| external | `study_outputs/external/endpoint_predictions.csv`, `pairs.csv`, `source_members.csv` | observations |

No route opens `outputs/reproducibility_delivery_20261109/`,
`outputs/code_data_release_20261115/` or any other author working directory.
The clean-room check in `VERIFICATION_MATRIX.md` enforces this with a Python
audit hook and a network deny.

## Output discipline

All outputs are written below the `--output` directory, outside the package.
Keep acquired observations and any outputs derived from them under their source
terms. The author outputs bundled here are covered by `LICENSE_SCOPE.md`.
