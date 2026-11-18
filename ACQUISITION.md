# Acquisition and import contract

This file records both routes for every non-author artifact:

1. **Local verification route** — the read-only copy used to validate this
   release. It is outside the public ZIP and is not redistributed.
2. **Reader route** — the official public or application route a reader can use,
   with the version/commit, identity fields and import command.

A provider landing page is not treated as proof of redistribution permission.
A current service page is not treated as the historical response version.

## 1. ProNAB measurements

| | |
|---|---|
| Local verification copy | `authorized_inputs/pronab.csv`, read-only; 17,608,829 bytes; SHA-256 `132cae5018808e73e4585324e7013e30f292fdceb08386e5a20b0e93582d81a3`; 21,206 rows |
| Reader route | Search/export or complete-snapshot download: <https://web.iitm.ac.in/bioinfo2/pronab/>. Cite Harini et al., DOI [10.1093/nar/gkab848](https://doi.org/10.1093/nar/gkab848). |
| Fields used | `Entry_id`; `Kd_wild (M)`; `Kd_mutant (M)`; the column named in `metadata/source_links.csv:raw_column` |
| Stable match | `Entry_id` **and** `source_row_number` **and** `raw_column`. `source_row_number` is a convenience for the frozen snapshot; if the provider export changes order, realign explicitly by `Entry_id` and verify the row/column value, never by row number alone. |
| Snapshot hash | The provider does not publish this study's snapshot SHA-256. The hash above identifies our frozen export; a reader with a different export can still match by `Entry_id` and the recorded column, then report any differing rows. |
| Import | `reproduce.py prepare --pronab <file>` or `reproduce.py primary --pronab <file>` |
| Rule | No ProNAB measurement row is bundled. The package ships only the selected link table (unit → source record/row/column), with no values. |
| Coverage | 846 selected links, 827 retained by the current membership rules; median −log10(Kd/M) per unit reproduces the 752 labels to 0.0 difference against the frozen archive. |

If a reader cannot obtain the exact historical snapshot, the reproducible path
is: download the current export, select rows by `Entry_id`, confirm
`metadata/source_links.csv:raw_column` exists and is a positive point value,
and run the importer. The importer rejects row/Entry-ID mismatches and bound
symbols rather than silently accepting a different record.

## 2. PRA-Pred fitted equations

| | |
|---|---|
| Local verification copy | extracted from the executed application (`outputs/pra_feature_mechanism_20261020/restricted/source/pdpredict_command_rev.py`); not bundled |
| Reader route | <https://github.com/Harinikannan98/PRA-Pred> — obtain the fitted equations from the source used by the paper. Executed source SHA-256: `667157288de06e6ef57c98f65c4c0ea2ec8be95a17d1c2cdb566046e79ebd88f`. |
| Why not bundled | The coefficients are third-party fitted parameters. The package publishes the author-extracted features and author-computed predictions; it does not redistribute the equation. |
| Import schema | `third_party/pra_formula_schema.json`: two branches `ss+str`, `ds+enz`; numeric `intercept` and `terms` map. The printed source term `charged` is the feature column named `charged_float` in the feature export. |
| Import | `--pra-formula <json>` on `prepare`, `full`, `all` and `reproduce.py full` |
| Verification | With the imported equation, all 128 endpoint values reproduce to `<1e-10` and the contribution sums to `2.0e-15`. Without it, PRA term contributions and Figure 5 cannot be rebuilt; the author predictions and the 42 exact-zero/13 rounding-tie counts are still present. |

## 3. mCSM-NA and PANTHER historical scores

| | mCSM-NA | PANTHER |
|---|---|---|
| Local verification | `authorized_inputs/services/mcsm/` (54 request TSVs, request manifest, scope table); read-only | `authorized_inputs/services/panther/` (110 JSON responses, ledger, cohort pairs); read-only |
| Provider | <https://biosig.lab.uq.edu.au/mcsm_na/>; paper DOI [10.1093/nar/gkx236](https://doi.org/10.1093/nar/gkx236), CC BY 4.0 article | <https://nova.disfarm.unimi.it/panther/>; author article <https://pmc.ncbi.nlm.nih.gov/articles/PMC12810179/> |
| Historical version | The recorded 2026 service responses. The current server code/version is not re-queried; a fresh request may differ and is **not** an exact historical route. | Same. The original site's current availability and response terms could not be verified; the author article documents the method, not a byte-identical historical endpoint. |
| Bundled | Per-pair numeric scores in `study_outputs/comparison/service_scores.csv`, request IDs, input SHA-256, response SHA-256, 14 historical-cache flags | Same table, 110 endpoint responses → 87 pairs |
| Not bundled | Response bodies, HTTP headers, cookies, credentials | Same |
| Rebuild if you have the archive | `replay_services.py --inputs <services> --output <dir>` recreates the pair scores and asserts the request/response identities | Same |
| Rebuild from the public package | `prepare_inputs.py` reads `service_scores.csv` and writes the same `service_pair_predictions.json` used by the comparison route | Same |
| Rule | The numeric scores are published as author-extracted analysis inputs; the response bodies are not redistributed because their terms were not established. No scientific POST is made by this package. |

The 267 mCSM pair scores and 87 PANTHER pair scores are sufficient for Table 3,
S37, S11 and S1/S3/S8/S15. They are not a substitute for the raw historical
response files if a reader wants to re-parse the provider's own text.

## 4. ProRB checkpoint

| | |
|---|---|
| Repository | <https://github.com/ymaa19/ProRB/tree/6b18f2d144251fc5b8a5ac465528afe07c25da08> |
| Model host | <https://huggingface.co/ymaaaa/ProRB/tree/4b0091890c761fcb610a4bbf941e88307a3a0c36> |
| Checkpoint | `author_bundle/model.pt`; 1,733,918,683 bytes; SHA-256 `8037368b9eba020021a23db8ca396da4fd693a8b95d9b65285721a326da5ba36` |
| License | The checkpoint/model-card license could not be verified from the checked sources. The checkpoint is **not** bundled. |
| What is bundled | The authors' own ProRB predictions (`predictions.json:prorb`, 406 units) and the comparison aggregates computed from them. These are author-generated outputs, judged separately from the checkpoint's redistribution terms. |
| Reader route | Obtain the checkpoint under the upstream terms if new inference is required. The paper replay does not load it; it consumes the saved predictions. |
| Identity | `model_identities/prorb.json` records the repo/HF revisions, checkpoint hash and RNA-FM adaptation hashes. |

## 5. M2-PRNet checkpoints and inputs

| | |
|---|---|
| Repository | <https://github.com/CSUBioGroup/M2-PRNet/tree/8353164d2d4f9244c6fce6ea6bab4e35710c08cd> |
| Model host | <https://huggingface.co/NikoWz/M2-PRNet/tree/1b895b2330ed056bc973d55ad64a3a01af418cf0> |
| Checkpoints | Five fixed folds; SHA-256 recorded in `study_outputs/comparison/m2_cohort_manifest.json:checkpoints` |
| License | No weight redistribution is asserted. Weights are not bundled. |
| Bundled | Author-generated M2 endpoint predictions for the original, mask-aligned and all-indices-aligned conditions; checkpoint deltas, fusion interventions, integrated-gradient readout, M2 reference metrics and the 211/131-pair scopes |
| Reader route | Obtain the five checkpoints from the hosts above for new inference. The released replay uses the saved outputs and performs no model execution. |

## 6. Other model assets

| Method | Recorded source | Redistributed? |
|---|---|---|
| CoPRA upstream | <https://github.com/hanrthu/CoPRA/tree/948debf7914db48faf16f7bcdce280601dd25415> | No upstream weights; the study's controlled fits are author outputs. Predictions bundled. |
| ESM2 | <https://huggingface.co/facebook/esm2_t33_650M_UR50D> | No weights |
| RiNALMo | <https://zenodo.org/records/15043668> | No weights |
| DeePNAP | <https://github.com/IntegrativeStructuralBiologyLab/DeePNAPWebsite/tree/911790b3cc5d43ecc9387996f5813288f0a59baa> | No weights; saved OOF predictions bundled |
| OmniBioTE | <https://github.com/nyuolab/OmniBioTE> | No weights; official/adapted predictions bundled |
| SVSBI | <https://github.com/WeilabMSU/SVS> | No weights; author-approved reconstruction predictions bundled |
| GBNL | <https://github.com/mzia-s/GBNL> | No weights; author-approved reconstruction predictions bundled |
| PRdeltaGPred | <https://doi.org/10.1002/prot.26503> | No source/parameters; the study's documented port predictions are bundled |
| NABE | <http://nabe.denglab.org/> | No data |
| UniProt | <https://www.uniprot.org/help/license> | No data |

## 7. wwPDB structures for Figure 1

| PDB | Bytes | SHA-256 | Official URL |
|---|---:|---|---|
| 2C4R | 376,488 | `f7f838f4b63cb6bb4e12800f5099df70a10f9dfba40fb63919948b30bd304125` | <https://files.rcsb.org/download/2C4R.pdb> |
| 2ZZN | 751,842 | `d39906841d66b60ada807df2768434e08342403136ece9fee40957bf9766b44a` | <https://files.rcsb.org/download/2ZZN.pdb> |
| 4QVC | 291,276 | `e3ec7d71e9b59712ba25a72319ac2146ff5417b377e1a242718fd4ba3e07166d` | <https://files.rcsb.org/download/4QVC.pdb> |
| 4NGD | 693,036 | `deb9caacc44be0155d7044a662ca8fb41005c7e4cb35bc2127311be0b22d5f6d` | <https://files.rcsb.org/download/4NGD.pdb> |

PDB archive coordinate files are CC0 1.0 under the wwPDB policy
(<https://www.rcsb.org/pages/policies>). They are not bundled because the paper
uses them only as rendering inputs; the fetch helper verifies bytes and hashes:

```sh
python -B third_party/fetch_structures.py --output ../imports/pdb
```

## 8. S13 external observations

| Family | Source | Location | Members |
|---|---|---|---:|
| SF3A1 | de Vries et al., DOI [10.1073/pnas.2114092119](https://doi.org/10.1073/pnas.2114092119), PMID 35101980 | Table 2 | 11 |
| SEPSECS | Puppala et al., DOI [10.1093/nar/gkad182](https://doi.org/10.1093/nar/gkad182), PMID 36929010 | Table 1 | 13 |
| Npl3 | Moursy et al., DOI [10.1038/s41467-023-42962-4](https://doi.org/10.1038/s41467-023-42962-4), PMID 37989742 | Supplementary Figure 7 | 7 |

Copy `study_outputs/external/observations_template.csv`, fill `kd_value` and
`kd_unit` for point rows, leave censored rows blank, and run:

```sh
python -B reproduce.py external --observations ../imports/observations.csv --output ../results/external
```

The importer accepts M, mM, uM, μM, nM and pM. The three families give 24
quantitative pairs from 28 attempted pairs and 31 members; the replay checks
64 family/macro cells, five paired contrasts and renders Figure E1.

## 9. Historical-version statement

The mCSM-NA and PANTHER responses are historical service outputs from the 2026
study run. Neither provider exposes a version-pinned replay of those exact
requests. The package therefore:

- bundles the author-extracted numeric scores and request identities;
- states that a new request to the current service is not an exact historical
  retrieval;
- does not present either provider's current home page as a historical archive;
- does not make a new scientific POST.

All other third-party artifacts (ProNAB, PDB, ProRB, M2, PRA-Pred) have an
official version or commit and a stable identity as recorded above.
