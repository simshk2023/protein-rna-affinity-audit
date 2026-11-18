# Acquisition and recorded identities

These links identify upstream resources; no weight file is bundled or
downloaded. The offline replay uses the prepared author outputs and needs none
of these runtimes. A checkpoint license is not a license for the authors'
predictions, and an author prediction is not a redistribution right for the
checkpoint.

| Configuration | Acquisition | Recorded identity / role |
|---|---|---|
| Controlled CoPRA/additive | <https://github.com/hanrthu/CoPRA/tree/948debf7914db48faf16f7bcdce280601dd25415> | Five family folds × three seeds × two methods; study-controlled fits are author outputs. Predictions bundled. |
| ESM2 / RiNALMo | <https://huggingface.co/facebook/esm2_t33_650M_UR50D> ; <https://zenodo.org/records/15043668> | Native `esm2_t33_650M_UR50D.pt` and `rinalmo_giga_pretrained.pt`; weights not bundled |
| DeePNAP / sequence comparator | <https://github.com/IntegrativeStructuralBiologyLab/DeePNAPWebsite/tree/911790b3cc5d43ecc9387996f5813288f0a59baa> | Saved OOF predictions bundled; weights not |
| OmniBioTE official/adapted | <https://github.com/nyuolab/OmniBioTE> | `omnibinder-small.pt` MD5 `728c3b5482fee12061918d01f9a75a1f`; `omnibiota-multi-E-mixed-small.pt` MD5 `9c556af61281bebeb40469eae10260bd`; weights not bundled |
| SVSBI reconstruction | <https://github.com/WeilabMSU/SVS> | Author-approved reconstruction, not substituted author weights; encoder fingerprints in `svsbi.json` |
| GBNL reconstruction | <https://github.com/mzia-s/GBNL> | Author-approved reconstruction; ESM2-t36-3B SHA and implementation hashes in `gbnl.json` |
| ProRB | <https://github.com/ymaa19/ProRB/tree/6b18f2d144251fc5b8a5ac465528afe07c25da08> ; <https://huggingface.co/ymaaaa/ProRB/tree/4b0091890c761fcb610a4bbf941e88307a3a0c36> | Checkpoint `model.pt` SHA-256 `8037368b9eba020021a23db8ca396da4fd693a8b95d9b65285721a326da5ba36`; 1,733,918,683 bytes; checkpoint license not verified; **not bundled**. The 406 author predictions are bundled. |
| PRdeltaGPred | <https://doi.org/10.1002/prot.26503> | Local port retained released parameters and recorded DSSP substitution; saved outputs and method scorer supplied; no unverified source revision assigned |
| PRA-Pred | <https://github.com/Harinikannan98/PRA-Pred> | Executed source SHA-256 `667157288de06e6ef57c98f65c4c0ea2ec8be95a17d1c2cdb566046e79ebd88f`; features/predictions bundled; fitted equations not bundled; new feature extraction requires FoldX/3DNA separately |
| M2-PRNet | <https://github.com/CSUBioGroup/M2-PRNet/tree/8353164d2d4f9244c6fce6ea6bab4e35710c08cd> ; <https://huggingface.co/NikoWz/M2-PRNet/tree/1b895b2330ed056bc973d55ad64a3a01af418cf0> | Five fixed checkpoints with SHA-256 in `m2_cohort_manifest.json`; weights not bundled; outputs and alignment/fusion/IG intermediates bundled |
| mCSM-NA / PANTHER | <https://biosig.lab.uq.edu.au/mcsm_na/> ; <https://nova.disfarm.unimi.it/panther/> | Historical 2026 request identities; numeric scores bundled; response bodies not bundled; current server versions are not exact historical routes |
| wwPDB | <https://files.rcsb.org/download/> ; <https://www.rcsb.org/pages/policies> | Four CC0 coordinate files; `third_party/fetch_structures.py` verifies bytes/SHA-256 |
| NABE / UniProt | <http://nabe.denglab.org/> ; <https://www.uniprot.org/help/license> | Provider identities only |

Upstream landing pages for OmniBioTE, SVSBI, GBNL and M2 were opened on
2026-11-10. Version identities above come from execution records, not an
assumption that today's default branch equals the executed source. Where the
historical record supplies only an implementation hash, that is retained
without inventing an upstream commit. Fresh training, preprocessing and service
availability were not tested in this delivery.

The full scope is 752 absolute units / 338 ordinary pairs / 220 local pairs;
the 771/354 as-fitted layer and the 7+12 censor/alias comparison are in
`MEMBERSHIP.md` and `RECONCILIATION_PUBLIC.md`.
