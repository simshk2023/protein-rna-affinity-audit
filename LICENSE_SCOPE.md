# License scope

The authors license their original analysis, replay and figure-generation code
under the MIT license in [LICENSE](LICENSE). The authors make the
study-generated outputs and metadata bundled in `study_outputs/` available for
reuse and redistribution, including the 4,512 endpoint predictions, the
771/354/846 identity layer, the method comparison blocks, the PRA/M2/three-arm
intermediates and the aggregate figure inputs. Please cite the paper and the
version used.

The MIT license applies to the authors' code. It does not relicense third-party
material. In particular:

- ProNAB measurement rows are not included;
- mCSM-NA and PANTHER response bodies are not included;
- ProRB, M2-PRNet, ESM2 and RiNALMo weights are not included;
- PRA-Pred fitted coefficients are not included;
- the S13 observation values are not included;
- the four wwPDB files are not included (the fetch helper verifies them).

Numeric service scores extracted by the authors, author predictions from
third-party checkpoints, and author features are included as study outputs.
Their use is subject to the source providers' terms described in
[THIRD_PARTY_ACCESS.md](THIRD_PARTY_ACCESS.md). The inclusion of an author
prediction does not imply a redistribution right for the underlying checkpoint
or provider rows.
