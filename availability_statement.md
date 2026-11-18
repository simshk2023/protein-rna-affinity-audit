# Data and code availability

## Current wording (v5 manuscript)

**Data availability.** The release package contains the study-generated
endpoint predictions for the full principal evaluation: 752 absolute-affinity
units, 338 ordinary mutation pairs and 220 local mutation pairs, together with
the unit, pair, source-link, censor and alias metadata, the method comparison
and intervention outputs, the PRA/M2/three-arm intermediates, and the aggregate
plot inputs for Figures 1–7. Primary affinity observations remain available
from ProNAB under its terms; the package ships the frozen source-link mapping
and imports the reader's ProNAB export rather than redistributing measurement
rows. The independent S13 measurements are identified by source article,
table/figure and construct; the package ships the endpoint predictions, pair
and member maps, and a blank observation template. Individual mCSM-NA and
PANTHER response bodies, third-party model weights and PRA-Pred fitted
coefficients are not redistributed. Numeric scores extracted from the
historical service responses, with request identities and hashes, are included
so the comparison tables can be recomputed; the historical service versions
themselves cannot be re-queried exactly. The replay recomputes Table 2, Table 3,
Supplementary Tables S36 and S37, Figures 1–7, Supplementary Figure E1, and the
PRA, M2 and three-arm diagnostics from the released predictions and the
declared imports. [ProNAB](https://web.iitm.ac.in/bioinfo2/pronab/)

**Code availability.** MIT-licensed analysis, replay and figure-generation code
and a fixed versioned archive are prepared for deposition. The clean-room
replay was verified on the fixed archive with the pinned CPU environment.
[Repository URL, tag and version DOI to be inserted after deposition.]

## At public deposition

Replace the bracketed sentence with the anonymous archive URL, tag, commit and
optional DOI, after the fixed ZIP is downloadable and its SHA-256 is recorded.
Do not describe the package as a "controlled subset": the released prediction
scope is the full 752/338/220 principal evaluation. The only withheld materials
are the third-party items listed above, and their acquisition routes are in
`ACQUISITION.md`.
