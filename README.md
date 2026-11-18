# Published study archive

Versioned release: https://github.com/simshk2023/protein-rna-affinity-audit/releases/tag/v2026.11.17

The release attachment is the unchanged approved archive. See [MANUSCRIPT_NUMBERING.md](MANUSCRIPT_NUMBERING.md) for correspondence with the full and compact manuscripts. The bundled manifest describes the frozen archive; publication documentation added here is outside that manifest.

# CoPRA paper replay · v2026.11.17

Fixed public replay for the v5 manuscript (15 November 2026). This version
supersedes the screened `v2026.11.16.1` candidate. It restores the **full**
author-generated prediction and metadata scope used by the paper:

| Population | Units / pairs |
|---|---:|
| Absolute affinity | **752 units** |
| Ordinary mutation effects | **338 pairs** |
| Local mutation effects | **220 pairs** (21 landscapes, 20 families) |
| Prediction rows (two methods × three seeds) | **4,512** |

The package contains author-generated endpoint predictions, unit/pair/source
metadata, comparison and intervention outputs, PRA/M2/three-arm intermediates,
and the aggregate plot inputs needed for Figures 1–7. It does **not** contain
ProNAB measurement rows, individual mCSM-NA/PANTHER response bodies,
third-party model weights, PRA-Pred fitted coefficients, or the S13 observation
values. Readers obtain those from the recorded provider routes and import them;
the supplied code then recomputes the paper results from per-record predictions
rather than from a pre-computed summary.

## Install and verify

Use 64-bit CPython 3.11.5 and the pinned CPU dependencies. The reference replay
used Windows, Arial, and the versions in `requirements.lock.txt`.

```sh
python -m venv .venv
# Windows: .venv\Scripts\activate    POSIX: source .venv/bin/activate
python -m pip install -r requirements.lock.txt
python -B reproduce.py smoke --output ../results/smoke
```

Verify the package bytes before running anything:

```sh
python -B reproduce.py prepare --pronab /path/to/pronab.csv \
    --pra-formula /path/to/pra_formula.json \
    --structures /path/to/pdb --output ../work
```

`prepare` verifies every file in `manifest.json`, reconstructs the 752 labels
from the imported ProNAB snapshot through `metadata/source_links.csv`, stages
the figures, and writes a local `work/` tree. The ProNAB file is not bundled.

## Reproduce the paper

```sh
python -B reproduce.py primary --pronab ../imports/pronab.csv --output ../results/primary
python -B reproduce.py full --pronab ../imports/pronab.csv \
    --pra-formula ../imports/pra_formula.json --structures ../imports/pdb \
    --output ../results/full
python -B reproduce.py external --observations ../imports/observations.csv \
    --output ../results/external
```

`all` runs the full route, the independent S13 route and the frozen-value
comparison in one call:

```sh
python -B reproduce.py all --pronab ../imports/pronab.csv \
    --pra-formula ../imports/pra_formula.json --structures ../imports/pdb \
    --observations ../imports/observations.csv --output ../results/all
```

`primary` recomputes Table 2 and Supplementary Table S36. `full` recomputes
Table 3 and S37, the PRA feature path, the M2 condition statistics, the
three-arm seed summaries, and Figures 1–7/S14. `external` recomputes S13 and
Figure E1 from the three source-matched panels. `verify_frozen.py` compares every
numeric result and every figure page against `frozen/reference.json`.

## What is and is not in this package

Bundled, author-generated:

- 4,512 endpoint predictions for the full 752-unit scope;
- 771 as-fitted unit identities, 354 candidate pairs, 846 source links (827
  retained), the censor/alias ledgers and the current membership rules;
- all method prediction blocks used by Table 3/S37, including ProRB, PRA-Pred,
  M2-PRNet and its five checkpoints;
- PRA features/predictions, M2 condition and intervention outputs, three-arm
  predictions and seed summaries, plot aggregates and renderer code;
- extracted per-pair numeric mCSM-NA/PANTHER scores with request IDs and
  SHA-256 identities, but no response bodies.

Not bundled, with an acquisition route in `ACQUISITION.md`:

- the ProNAB measurement CSV (third-party observations);
- mCSM-NA and PANTHER response bodies (third-party service outputs);
- ProRB, M2-PRNet, ESM2 and RiNALMo weights/checkpoints;
- PRA-Pred fitted equations;
- the S13 observation values (source articles);
- four wwPDB coordinate files used by Figure 1.

Those are genuine source restrictions or third-party terms. They are not a
restriction on the authors' own predictions. In particular, the previous
`screened` 672/297/199 scope was not supported by a provider term or a source
field; `RECONCILIATION_PUBLIC.md` documents the 80-unit/41-pair/21-local-pair
restoration and the evidence.

## Layout

```
metadata/                 as-fitted units, pairs, source links, exclusions
study_outputs/primary/    full prediction table + 752/338/220 tables
study_outputs/comparison/ all method prediction blocks + service scores
study_outputs/mechanisms/ PRA/M2/three-arm archived intermediates
study_outputs/figures/    aggregate plot inputs for Figures 1-7
study_outputs/external/   S13 predictions, pairs, members, blank template
manuscript_v5/            table snapshots and alignment drawing
frozen/                   frozen reference values and figure page identities
third_party/              acquisition identities, PDB fetch helper, formula schema
figure_code/              archived renderers
toolkit/                  prediction/observation schema and scoring core
```

## License and terms

The authors license their original code under MIT (`LICENSE`). The bundled
study outputs are released for reuse and redistribution under the scope in
`LICENSE_SCOPE.md`; cite the paper and the version. Third-party materials
remain with their providers (`THIRD_PARTY_ACCESS.md`, `NOTICE.md`).
