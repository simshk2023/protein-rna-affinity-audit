# Input contract and provenance

The package is self-contained for every author-generated value. Four external
artifacts are imported explicitly; see `ACQUISITION.md`. No historical project
directory is needed.

## Primary reconstruction

`metadata/` is the as-fitted and selected source-identity layer:

| File | Rows / records | Fields | Observed values? |
|---|---:|---|---|
| `units.csv` | 771 as-fitted units | `unit_id`, `family_id`, `fold`, `side`, `landscape_id`, `mutation_tokens`, `mutation_kind`, `structure_key`, `input_key` | **No** |
| `pairs.csv` | 354 candidate pairs | `pair_id`, `landscape_id`, `family_id`, `fold`, `wt_unit_id`, `mutant_unit_id`, `structure_key` | **No** |
| `source_links.csv` | 846 selected links | `unit_id`, `source_record_id`, `source_row_number`, `entry_id`, `raw_column` | **No** |
| `candidate_membership.csv` | 354 rows | `pair_id`, `landscape_id`, `family_id`, `position`, `old_shape`, `new_shape`, `reason`, `new_position_weight` | **No** |
| `censor_exclusions.json` | 7 units / 6 pairs | frozen bound exclusions | **No** |
| `source_alias_exclusions.json` | 12 units / 10 pairs | Hfq duplicate-alias rule | **No** |

The importer reads each retained link's ProNAB `Kd_wild (M)` or `Kd_mutant (M)`
column, converts to pKd = −log10(Kd/M), and takes the median within a unit. The
current evaluation is 771 − 7 censored − 12 aliases = **752 units**;
354 − 6 censored − 10 alias pairs = **338 pairs**; the B-shape subset is
**220 pairs** (21 landscapes, 20 families). Training manifests remain the
as-fitted 771-unit layer; the corrections change evaluation only.

`study_outputs/primary/` mirrors the current evaluation:

| File | Rows | Notes |
|---|---:|---|
| `predictions.csv` | 4,512 | `method, seed, member_id, status, value, checkpoint`; two methods × three seeds × 752 units |
| `units.csv` | 752 | evaluation unit identities |
| `pairs.csv` | 338 | evaluation pairs, no observed delta |
| `source_links.csv` | 827 | retained links after the censor/alias rules |
| `study_identity.json` | — | counts, ProNAB bytes/SHA-256, file hashes |

## Comparison outputs

`study_outputs/comparison/`:

- `predictions.json` — 17 method blocks (`copra`, `baseline`, `svsbi`, `gbnl`,
  `prorb`, `prdelta`, `pra`, `m2`, `m2_ck0..4`, `deeppnap`,
  `seq_baseline`, `omni_official`, `omni_controlled`), each with its seed keys
  and unit coverage.
- `metadata.json` — `D` (771 unit records) and `P` (354 pair records) with the
  observed `label` / `delta_pkd` fields removed; `ex` and `routes` retained.
- `pra_inputs.json` — 128 endpoints with features and author predictions; the
  fitted `formula` strings and absolute source paths are removed.
- `m2_predictions.json`, `m2_study_manifest.json`,
  `m2_cohort_manifest.json` — saved M2 endpoint outputs, conditions,
  checkpoint identities, pair/endpoint manifests.
- `three_arm_predictions.json`, `three_arm_figure7_reference.json` — three-arm
  prediction and seed-summary layers.
- `service_scores.csv` — 87 PANTHER and 267 mCSM pair scores with request IDs,
  input/response SHA-256 and cache flags; no response bodies.
- `service_manifest.json` — counts, request IDs, provider links and the
  historical-version rule.

`prepare_inputs.py` merges `predictions.json` and `metadata.json`, adds the
reconstructed labels/deltas, and writes the `saved_inputs.json` shape used by
`replay_methods.py`. It never invents a label; a missing or mismatched link
raises.

## Mechanism inputs

`study_outputs/mechanisms/`:

- `training.json`, `validation.json` — 451/200/124 and 164/70/54 scopes,
  arm/seed metrics and family arrays; checkpoint paths are reduced to logical
  basenames.
- `diagnostics.json` — rank/amplitude cross-checks and the 200/124-pair
  decomposition.
- `m2_reference.json` — 211/131-pair scopes, five checkpoints, condition
  metrics, intervention effects.
- `pra_inputs.json` — same formula-free PRA layer as the comparison directory.

## External S13

`study_outputs/external/`:

| File | Content |
|---|---|
| `endpoint_predictions.csv` | 930 rows: family, method, seed, fold, unit_id, prediction |
| `pairs.csv` | 28 attempted pairs across SF3A1, SEPSECS and Npl3 |
| `source_members.csv` | 31 members with PMID and exact source table/figure |
| `observations_template.csv` | blank reader input: unit_id, status, kd_value, kd_unit |
| `study_identity.json` | family/member/pair counts and file hashes |

## Figures

`study_outputs/figures/` contains aggregate plot inputs only. No file contains a
per-pair observed value or an invertible row-level residual. Figure 5's
per-pair *features* are bundled; its *contributions* are recomputed at prepare
time from the imported PRA equation.

## Prediction schema

`toolkit/SCHEMA.md` defines the endpoint prediction table. A valid row has a
method, a seed/fixed label, a unit `member_id`, `status=success`, a finite pKd
value and a checkpoint identity. The package's `predictions.csv` and
`predictions.json` follow that schema; `replay_methods.py` intersects methods by
their available members and never imputes a missing prediction.
