# Figure reproduction

The manuscript uses the archived renderers' own page output. `collect_figures.py`
copies those pages as `figure1.pdf` … `figure7.pdf`; `compose_v5.py` is retained
only for the supplementary S14 layout. This is the v5 layout: the raw Figure 3
page is 504.6 × 564.1 pt with panels A–D, not the older two-panel composition.

```sh
python -B reproduce.py full --pronab ../imports/pronab.csv \
    --pra-formula ../imports/pra_formula.json --structures ../imports/pdb \
    --output ../results/full
```

| Figure | Renderer | Inputs | Reader import |
|---|---|---|---|
| 1 | `figure_code/figure1/make_framing_grounded_20261107.py` | archived ribbon/contact PNGs and the 1A9N complex render | four wwPDB files for the structure parser |
| 2 | `figure_code/figure2/render_figure1_wide.py` + `render_response_wide.py` | `family_aggregates.json` (full 113-family table), `bootstrap_regions.json`, `plot_coordinates.json` | ProNAB-derived labels not needed at render time; the frozen aggregate is bundled |
| 3 | `figure_code/figure3/render_figure2_wide.py` | full `plot_data.json`: binned masses, 52-family arrays, residual KDE grid, 2,000 bootstrap draws | ProNAB labels are not embedded as rows; the aggregate layer is bundled |
| 4 | `figure_code/figure4/render_figure3.py` | `plot_data.json`: 220-pair family arrays, supervision components, joint bootstrap | none |
| 5 | `figure_code/figure5/render_figure4.py` | `plot_data.json`, `covariance_data.json`, `ellipse_data.json`, and `restricted/pair_plot_data.json` built by `prepare_inputs.py` from the bundled features and the imported PRA equation | PRA-Pred fitted equation |
| 6 | `figure_code/figure6/render_figure5_short.py` | `plot_data.json`, `plot_data_D.json` | none |
| 7 | `figure_code/figure7/render_figure6.py` | `plot_data.json`, `family_seed_metrics.csv`, `landscape_seed_metrics.csv` | none |
| S14 | `compose_v5.py` from the Figure 6 renderer | same Figure 6 inputs | none |
| E1 | `plot_external.py` via `reproduce.py external` | external endpoint predictions + observation template | S13 source values |

## What changed from the previous candidate

- Figure 2 no longer requires an "authorized figures" directory: the full
  `family_aggregates.json` (697,106 bytes, SHA-256 `6dfeefd6…`) is bundled.
- Figure 5 still separates the third-party equation. The package bundles the
  per-pair feature vectors and the branch statistics; `prepare_inputs.py`
  recomputes the exact term contributions and the A/B tracks from the imported
  equation. This reproduces 42 exact zeros and 13 rounding-only ties.
- Figure 7 bundles the family/landscape seed CSVs and the full `plot_data.json`;
  no archived cache is requested from the reader.
- Figure 3 uses the full A–D page. The previous two-panel composition is not
  the v5 manuscript figure.
- Figure 1 needs the four PDB files only for the renderer's structure parser;
  the archived PNG layers are bundled.

## Verification notes

`verify_frozen.py` checks the page size and the full extracted text-token
multiset against the manuscript PDF for every main figure. The reference run
matched Figures 1, 2, 3, 5, 6 and 7 byte-for-byte and Figure 4 in page size and
text (the PDF byte stream differs because of font/compression metadata without a
visual change). Platform or font differences can change raster pixels; the
numeric values and layouts are the acceptance target.
