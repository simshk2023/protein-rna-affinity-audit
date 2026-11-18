# Acquire and import the ProNAB export

The package ships the author-curated link table, not the provider measurement
rows. Every retained row is identified by `Entry_id`, `source_row_number` and
the exact `Kd_wild (M)` / `Kd_mutant (M)` column.

1. Open <https://web.iitm.ac.in/bioinfo2/pronab/> and use **Download Data** for
   a complete export or the search/export route for individual records.
2. Keep the provider column names and `Entry_id` values.
3. Pass the CSV to the importer:

```sh
python -B reproduce.py primary --pronab /path/to/pronab.csv --output ../results/primary
```

The importer refuses:

- a row whose `Entry_id` differs from `metadata/source_links.csv`;
- a value containing `<`, `>`, `≤`, `≥` or `~` (a bound, not a point);
- a non-finite or non-positive Kd.

If the export order changed, realign by `Entry_id` and the recorded column; do
not keep the old row number. If a value differs from the frozen snapshot, report
the affected `Entry_id`/unit. The frozen export used here is 17,608,829 bytes,
SHA-256 `132cae5018808e73e4585324e7013e30f292fdceb08386e5a20b0e93582d81a3`,
21,206 rows.

The resulting 752 labels reproduce the frozen archive to 0.0 maximum absolute
difference; the 827 retained links, 7 censored units and 12 Hfq alias units are
documented in `../MEMBERSHIP.md`.
