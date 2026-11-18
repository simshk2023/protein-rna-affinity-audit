# Verification matrix

Target: fixed package `v2026.11.17`, paper v5 (15 November 2026). The reference
run is CPU-only; it executes no training, inference or service request.

Environment: 64-bit CPython 3.11.5, Windows; NumPy 1.26.4, SciPy 1.16.3,
Matplotlib 3.10.8, pandas 2.3.3, PyMuPDF 1.26.7, Pillow 12.1.0; Arial.

| # | Check | Command / evidence | Result |
|---:|---|---|---|
| 1 | Package bytes match `manifest.json` | `reproduce.py` before any route | PASS, every file |
| 2 | Synthetic tests | `reproduce.py smoke` | PASS, 9 tests |
| 3 | ProNAB reconstruction | `reproduce.py primary` | PASS, **752 units / 338 pairs / 827 links**; max label difference 0.0 vs frozen archive |
| 4 | Main Table 2 | `check_table2.py` | PASS, **9 rows / 44 numeric cells**, 4-decimal manuscript agreement |
| 5 | Supplementary Table S36 | `check_v5_tables.check_primary` | PASS, 2 linked rows + 2 zero contrasts |
| 6 | Main Table 3 | `check_table3.py` | PASS, **10 configurations**, 30 displayed cells |
| 7 | Supplementary Table S37 | `check_table3.py` | PASS, **37 rows** |
| 8 | PRA term path | `replay_pra_mechanism.py --formula` | PASS, **128 endpoints**; 106 main pairs / 121 units / 15 landscapes / 14 families; 42 exact zeros, 55 rounded zeros, 13 rounding-only; max formula error 0.0; max contribution-sum error 2.0e-15 |
| 9 | M2 condition statistics | `replay_mechanisms.py` | PASS, **36 condition rows / 72 metrics**; max difference 3.33e-16 |
| 10 | Three-arm summaries | `replay_mechanisms.py` | PASS, **30 cells / 90 seed values**; max difference 0.0 |
| 11 | Figures 1–7 | `build_figures.py` + `collect_figures.py` | PASS, page size and full text-token multiset match the manuscript; Figures 1, 2, 3, 5, 6, 7 byte-identical, Figure 4 same size/text with a different PDF byte stream |
| 12 | Supplementary S14 | `compose_v5.py --figures S14` | PASS, rendered |
| 13 | S13 / Figure E1 | `reproduce.py external` | PASS, **3 families, 31 members, 28 attempted pairs, 24 quantitative pairs, 64 cells, 5 paired contrasts**; max difference 1.11e-16 |
| 14 | Frozen-value comparison | `verify_frozen.py` | PASS, **160 checks, 0 failed** |
| 15 | Clean-room unpack and replay | `clean_room.py` on the fixed ZIP | PASS, fresh temp extraction; audit hook active in every Python process; zero reads under the original project tree; network denied |
| 16 | Output discipline | clean-room audit | all outputs under the selected `--output`; no package write |

## Figure notes

- The v5 Figure 3 is the raw renderer page (504.6 × 564.1 pt, panels A–D),
  not the older two-panel composition. The package's `collect_figures.py`
  copies the renderer page directly.
- Figure 4 has the same page size, dimensions and extracted text as the
  manuscript; its PDF byte stream differs because of font/compression metadata.
  Raster pixels can differ with fonts/platform; numeric and layout checks pass.
- Figure 1's renderer parses the four wwPDB files; the package verifies their
  bytes/SHA-256 through `third_party/fetch_structures.py`.

## What is not claimed

- No new model training or inference is performed.
- No new mCSM-NA/PANTHER request is made; the historical service version cannot
  be re-queried exactly.
- A fresh ProNAB export with different content is matched by `Entry_id` and the
  recorded column; the importer reports a mismatch rather than silently
  accepting it.
- The exact historical mCSM-NA/PANTHER response bodies are not redistributed;
  their extracted numeric scores are sufficient for the paper comparisons.
