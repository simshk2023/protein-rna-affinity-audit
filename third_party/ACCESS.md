# Third-party input index

| Artifact | Official route | Recorded identity | Import |
|---|---|---|---|
| ProNAB measurements | <https://web.iitm.ac.in/bioinfo2/pronab/> | 17,608,829 bytes; SHA-256 `132cae5018808e73e4585324e7013e30f292fdceb08386e5a20b0e93582d81a3`; 21,206 rows | `--pronab` |
| PRA-Pred equations | <https://github.com/Harinikannan98/PRA-Pred> | executed source SHA-256 `667157288de06e6ef57c98f65c4c0ea2ec8be95a17d1c2cdb566046e79ebd88f` | `--pra-formula` |
| wwPDB 2C4R, 2ZZN, 4QVC, 4NGD | `third_party/fetch_structures.py` and `checksums.json` | per-file SHA-256 in `checksums.json` | `--structures` |
| S13 observations | `study_outputs/external/ACCESS.md` | three articles, 31 members, 28 attempted pairs | `--observations` |
| ProRB checkpoint | <https://github.com/ymaa19/ProRB/tree/6b18f2d144251fc5b8a5ac465528afe07c25da08>; <https://huggingface.co/ymaaaa/ProRB/tree/4b0091890c761fcb610a4bbf941e88307a3a0c36> | SHA-256 `8037368b9eba020021a23db8ca396da4fd693a8b95d9b65285721a326da5ba36` | not needed; predictions bundled |
| M2-PRNet checkpoints | <https://github.com/CSUBioGroup/M2-PRNet/tree/8353164d2d4f9244c6fce6ea6bab4e35710c08cd>; <https://huggingface.co/NikoWz/M2-PRNet/tree/1b895b2330ed056bc973d55ad64a3a01af418cf0> | five fold SHA-256 in `study_outputs/comparison/m2_cohort_manifest.json` | not needed; outputs bundled |
| mCSM-NA / PANTHER | <https://biosig.lab.uq.edu.au/mcsm_na/>; <https://nova.disfarm.unimi.it/panther/> | 54/110 historical requests; per-pair IDs and hashes in `service_scores.csv` | not needed; scores bundled |

No third-party weight or raw response body is in this directory. The full
acquisition contract is [ACQUISITION.md](../ACQUISITION.md).
