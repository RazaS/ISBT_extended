# Python App

This is a Python recreation of the Shiny app using only Python standard library (`http.server`) plus browser-side JavaScript.

## Run

```bash
cd /Users/sraza/Library/CloudStorage/OneDrive-UniversityofToronto/CBS_work/Geno_research/Geno_database_R_Mockup
python3 python_app/app.py
```

Then open:

- http://127.0.0.1:8050

## Notes

- Reads data from:
  - `Gene_table.csv`
  - `Allele_table.csv`
  - `Variant_table.csv`
  - `Exon_table.csv`
  - `Bridge_table.csv`
- Stores feedback submissions in process memory.
- Implements:
  - A1 mode switches: `Gene Table`, `Allele Table`, `Exon Table`
  - Allele mode split layout with gene search + selectable gene list
  - Allele/Variant detail toggle and global text filtering
  - Bridge popups:
    - allele row double-click -> linked variants
    - variant row click -> linked alleles
  - Row checkboxes on all tables (main + popup), with `Download CSV` and `Clear All`
  - `Ref_allele_curated` and `Alt_allele_curated` truncation (`...`) with hover tooltip
  - Full-column display with horizontal scrolling
  - Feedback form with reactive app-state autofill
