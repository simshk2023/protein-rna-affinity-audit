# Acquire and import the external S13 observations

The three families are source-matched, not new experiments. Use
`source_members.csv` to map each `unit_id` to its exact source token, PMID and
table/figure.

| Family | Source | Location | Members | Point members |
|---|---|---|---:|---:|
| SF3A1 | de Vries et al., DOI [10.1073/pnas.2114092119](https://doi.org/10.1073/pnas.2114092119), PMID 35101980 | Table 2 | 11 | 8 |
| SEPSECS | Puppala et al., DOI [10.1093/nar/gkad182](https://doi.org/10.1093/nar/gkad182), PMID 36929010 | Table 1 | 13 | 12 |
| Npl3 | Moursy et al., DOI [10.1038/s41467-023-42962-4](https://doi.org/10.1038/s41467-023-42962-4), PMID 37989742 | Supplementary Figure 7 | 7 | 7 |

Copy `observations_template.csv` outside the package. Fill `kd_value` and
`kd_unit` for each `point` row; leave both empty on a
`censored_or_nonbinding` row. The importer accepts M, mM, uM, μM, nM and pM and
converts each point to pKd. Supply one row for each of the 31 source members.
Values transcribed from rounded figures can fail an exact numerical comparison;
that is a source-transcription difference, not a package failure.

```sh
python -B reproduce.py external --observations ../observations.csv --output ../results/external
```

The code averages the five prepared folds within each family, method and seed,
then scores the three seeds. The output checks S13's 24 quantitative pairs from
28 attempted pairs, 64 family/macro cells and five paired contrasts, and renders
Figure E1. No provider measurements, raw fitting inputs or weights are in this
directory.
