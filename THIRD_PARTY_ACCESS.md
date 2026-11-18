# Third-party terms and separation of authorship

The package separates three things that the previous candidate conflated:

1. **Author-generated outputs** — predictions, features, aggregates, figures.
   These are the authors' own results and are bundled.
2. **Third-party observations** — ProNAB measurement rows and the S13 source
   values. These are not bundled; readers import them from the providers.
3. **Third-party models and service products** — weights, fitted PRA
   coefficients, and raw mCSM-NA/PANTHER response bodies. These are not
   bundled; the authors' numeric outputs computed from them are bundled where
   the output itself is an author analysis result.

A restrictive third-party term on a checkpoint does not make the authors'
predictions from that checkpoint restricted. Conversely, an author prediction
does not grant a right to redistribute the checkpoint or the raw provider rows.

| Resource | Terms / status | Bundled | Not bundled |
|---|---|---|---|
| ProNAB | Provider describes educational/research use; the exact snapshot's redistribution terms were not established. Cite Harini et al., DOI 10.1093/nar/gkab848. | Unit/pair/source-link mapping, no values | 21,206 measurement rows |
| mCSM-NA | Paper CC BY 4.0; individual historical server-response redistribution not established. | 267 author-extracted pair scores, request IDs, hashes | 54 response TSVs, headers, cookies |
| PANTHER | Author article documents the service; current site availability and response terms could not be verified. | 87 author-extracted pair scores, request IDs, hashes, 14 cache flags | 110 response JSON bodies |
| PRA-Pred | Public code repository; fitted coefficients are third-party parameters. | Features, author predictions, branch statistics | Fitted equations; reader imports them |
| ProRB | Repo/HF revisions and checkpoint SHA recorded; checkpoint license not verified. | Author predictions for 406 units | 1.73 GB checkpoint and encoder weights |
| M2-PRNet | Repo/HF revisions and five checkpoint hashes recorded; weights not redistributed. | Predictions, condition deltas, interventions | Five checkpoints and preprocessing assets |
| CoPRA upstream | MIT upstream code; the study's controlled fits are author outputs. | Controlled predictions | Upstream/controlled weights |
| ESM2 / RiNALMo | MIT / Apache-2.0 source notices; weights not redistributed. | Predictions | Weights |
| DeePNAP / OmniBioTE / SVSBI / GBNL / PRdeltaGPred | Method-specific terms; only author-generated predictions and identities are included. | Predictions and identities | Weights, source revisions, parameters |
| wwPDB | CC0 1.0 archive policy. | Nothing | Four coordinate files (fetch helper verifies hashes) |
| S13 articles | Article terms; values transcribed by the reader. | Endpoint predictions, pairs, member map, blank template | Observation values and source PDFs |

## Service responses: what is public and why

The public `study_outputs/comparison/service_scores.csv` contains the numeric
score difference per pair, the request IDs, the provider input SHA-256 values,
and the response SHA-256 values. These are the authors' extracted analysis
inputs, not the provider's response text. The raw bodies remain in the
authors' authorized archive and are identified by request ID and hash. A reader
who has the matching archive can rerun `replay_services.py`; a reader without
it can still run the paper's comparison tables from the extracted scores.

The historical service versions cannot be re-queried exactly. A fresh request
to the current mCSM-NA or PANTHER page is not an exact historical retrieval and
is not presented as one. No new scientific POST is made by this package.

## What "not bundled" means for reproducibility

The full scientific results in Table 2, Table 3, S36, S37, Figures 1–7, S13,
E1, PRA, M2 and the three-arm study are reproducible from
(a) the bundled author predictions and aggregates and
(b) the four reader imports listed in `data_access.md`.
No result requires a private cache. The only results that require the fitted
PRA equation are the term-contribution path and Figure 5; that equation has a
public source and an import schema.

## Attribution

See `NOTICE.md`, `third_party_notices/`, and
`model_identities/{prorb,svsbi,gbnl}.json` for the recorded revisions, byte
sizes and SHA-256 values. The author code is MIT; see `LICENSE` and
`LICENSE_SCOPE.md`.
