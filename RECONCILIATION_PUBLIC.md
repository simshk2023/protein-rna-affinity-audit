# Reconciliation: authorized archive → fixed public package

The fixed package `v2026.11.17` restores the full paper scope. The prior
`v2026.11.16.1` candidate released 672/297/199 because an undocumented v1
identifier screen removed 80 units; it was not a provider term or a field
restriction. The full machine-readable tables are `reconciliation.csv`,
`reconciliation_field_level.csv` and, for the local internal copy,
`restricted/reconciliation_omitted_{units,pairs}.csv`.

## 1. Population changes

| Population | Archive | Fixed public | Prior open | Difference | Explanation |
|---|---:|---:|---:|---:|---|
| Absolute units | 771 as fitted | **752** | 672 | +80 | 7 censored + 12 Hfq aliases = 752; the prior 80-unit screen is removed |
| Candidate pairs | 354 as fitted | **338** | 297 | +41 | 6 censored + 10 Hfq aliases = 338; 41 pairs restored |
| Local pairs | 230 candidate B-shape | **220** | 199 | +21 | 10 alias pairs leave 220; 21 restored |
| Prediction rows | 4,512 | **4,512** | 4,032 | +480 | 80 units × 2 methods × 3 seeds |

The 80-unit difference is exactly the set listed by the prior v1 screen:
**80 units, 26 landscapes, 15 families**, 54 mutant and 26 WT records, all on the protein side. All 21 omitted local pairs are a subset of the 41 omitted ordinary pairs. 40 of the 41 pairs have both endpoints in the 80; 1 has one endpoint in the 80. No pair was omitted by a changed B-shape threshold and no unit was omitted for a second reason.

The screen file is `open_release_screen_internal.json` (8,880 bytes, SHA-256 `2b988d3c6a39eff9e11499f1685d968b43bf24b3ef9efdb6a19e118eccb2de81`).
It has exactly two keys, `unit_ids` and `pair_ids`; it has no provider
term, no author instruction, no reason field, no field
marking a record as controlled, and no source-level access condition. The
records are ordinary public ProNAB entries from 15 studies; the predictions
are author-generated. Treating that list as a restriction would expand the
withheld scope without evidence. Genuine restrictions are handled separately
and are listed in section 3.

### Restored categories

| Protein | Source DOI(s) | Landscapes | Units restored |
|---|---|---:|---:|
| tRNA  | http://dx.doi.org/10.1038/nsmb.3282 | 1 | 22 |
| Nucleoprotein  | http://dx.doi.org/10.1074/jbc.M114.604678 | 2 | 14 |
| Genome polyprotein  | http://dx.doi.org/10.1128/JVI.03283-14 | 2 | 12 |
| Gag polyprotein  | http://dx.doi.org/10.1021/bi047639q; http://dx.doi.org/10.1038/nature02944 | 7 | 7 |
| Non-structural protein 1  | http://dx.doi.org/10.1038/cr.2008.288 | 1 | 5 |
| Gag-Pol polyprotein  | http://dx.doi.org/10.1006/jmbi.2000.3979; http://dx.doi.org/10.1038/nature13709 | 4 | 4 |
| Genome polyprotein Cleaved into: P3; Protein 3AB; P2; P1; Capsid protein VP0  | http://dx.doi.org/10.1093/nar/gkz1170 | 1 | 4 |
| Polymerase acidic protein | http://dx.doi.org/10.1038/nature20594 | 1 | 3 |
| Suppressor of silencing 2b  | http://dx.doi.org/10.1038/embor.2008.118 | 1 | 3 |
| 30S ribosomal protein S8 | http://dx.doi.org/10.1093/nar/gku743 | 1 | 1 |
| Capsid protein p27, alternate cleaved 1  | http://dx.doi.org/10.1016/j.jmb.2005.04.046 | 1 | 1 |
| Non-structural protein 3  | http://dx.doi.org/10.1016/s0092-8674(01)00632-8 | 1 | 1 |
| Polymerase cofactor VP35 | http://dx.doi.org/10.1128/JVI.01452-13 | 1 | 1 |
| Polymerase cofactor VP35  | http://dx.doi.org/10.1038/nsmb.1765 | 1 | 1 |
| RNA silencing suppressor p19  | http://dx.doi.org/10.1016/s0092-8674(03)00984-x | 1 | 1 |

## 2. File and field status

The archive-to-public mapping has 233 rows. Summary of the A-tree status:

| Status | Files | Meaning |
|---|---:|---|
| included | 43 | bundled unchanged or with path-only sanitization |
| split | 9 | author fields bundled; third-party labels, coefficients or response bodies removed |
| external | 169 | not redistributed; provider route recorded |
| not_applicable | 0 | not part of the public route |

Key splits:

- `methods/saved_inputs.json`: all 17 method prediction blocks are bundled;
  `D[*].label` and `P[*].delta_pkd` are removed and reconstructed from the
  reader's ProNAB import.
- `methods/pra_inputs.json`: features and author predictions are bundled;
  the fitted `formula` strings are removed and imported separately.
- `figures/figure5/restricted/pair_plot_data.json`: per-pair features are
  bundled; contributions are recomputed at prepare time from the imported
  equation (42 exact zeros, 13 rounding ties verified).
- `pronab.csv`, mCSM/PANTHER response bodies and the four PDB files are
  external; the package carries the mapping, numeric scores and hashes.

## 3. Genuine restrictions retained

| Material | Why not bundled | Legal reader route | Effect on results |
|---|---|---|---|
| ProNAB rows | third-party observations; snapshot redistribution terms not established | provider download + 827-link mapping | none; labels reconstruct exactly |
| mCSM-NA/PANTHER bodies | individual historical service responses; redistribution not established | extracted scores bundled; archive route for holders | none; Table 3/S37/S11 recompute |
| ProRB checkpoint | checkpoint license not verified | upstream repo/HF; predictions bundled | none; ProRB predictions are author outputs |
| M2 weights | third-party weights | upstream repo/HF; outputs bundled | none; M2 outputs are author outputs |
| PRA-Pred coefficients | third-party fitted parameters | public repo + import schema | Figure 5/PRA term path only |
| S13 values | source-article terms | fill the blank template | none; S13/E1 recompute |
| wwPDB files | structural inputs; CC0 but not needed in the package | fetch helper with SHA-256 | Figure 1 rendering |

## 4. Prior open package changes

| Prior item | Fixed item | Change |
|---|---|---|
| 672 units / 297 pairs / 199 local | 752 / 338 / 220 | restored the unsupported screen |
| Table 2 route needed the original controlled archive | `reproduce.py primary` | per-record predictions + ProNAB import; 9 rows/44 cells |
| Table 3 route covered S3/S4/S6/S9/S10 | `reproduce.py full` | S3/S4/S6/S7/S8/S9/S10/S11, ProRB + service scores |
| Figures 2/5/7 needed an authorized figure tree | bundled aggregates + formula import | no author-only cache |
| `endpoint_predictions.json` 10 methods | `predictions.json` 17 blocks | full author prediction inventory |
| `membership.json` subset | `metadata.json` 771/354 + runtime labels | full provenance |

The figure route now uses the v5 raw renderer pages. Figure 3 is the full
A–D page (504.6 × 564.1 pt), not the older two-panel composition. Figures
1–7 match the manuscript page size and text token multiset; Figure 4's PDF
byte stream differs without a layout change.
