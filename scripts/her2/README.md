# HER2 paralog synthetic-lethality analysis

Scripts behind the HER2 section of *Systematic prioritisation of context-specific
paralog pair vulnerabilities in cancer* (Daldal et al. 2026).

## Input data

The prediction matrix is not held in the repo. Download it from Figshare and place
the unzipped CSV at `data/input/her2/full_SLprediction_matrix.csv`:

<https://doi.org/10.6084/m9.figshare.31058182>

> **Use v3 or later.** v3 contains 33,419 paralog pairs; v1 and v2 contain 36,140,
> including 2,721 pairs predictable only in HeLa or PC9. The variance filter below
> keeps the top 10% of pairs *by matrix width*, so running against v1/v2 gives
> 3,612 tested pairs instead of 3,341 and a hit set differing by two pairs — the
> numbers from an earlier draft, not the published ones. `her2_ttest_top10pctvar.py`
> prints the matrix dimensions on startup; check it reports **33,419 pairs**.

Also required in `data/input/her2/`:

| File | Source |
|---|---|
| `DepmapHer2_withDai.xlsx` | Breast HER2 status (consensus call) |
| `Oesophageal_HER2_CN.xlsx` | Gastro-oesophageal HER2 status (copy-number call) |
| `kegg_erbb.txt` | KEGG ErbB pathway, hsa04012 (86 genes) |
| `reactome_erbb2.txt` | Reactome "Signaling by ERBB2", R-HSA-1227986 (50 genes) |
| `BIOGRID-GENE-108376-5.0.258.tab3.txt` | BioGRID ERBB2 interactions (839 partners) |
| `HINT-search-results-main.tsv` | HINT ERBB2 interactions (321 partners) |

`excluded_pairs.csv` in this directory lists the 2,721 pairs dropped between v2 and
v3, with their HeLa/PC9 predictions and, for the 20 that were experimentally
screened, the Parrish or CHyMErA genetic-interaction call.

## Scripts

Run `her2_ttest_top10pctvar.py` first — every other script reads its output.

```bash
python3 scripts/her2/her2_ttest_top10pctvar.py
python3 scripts/her2/erbb_pathway_enrichment.py
python3 scripts/her2/erbb2_reactome_enrichment.py
python3 scripts/her2/erbb2_interactor_gene_enrichment_efilter.py
python3 scripts/her2/oesophageal_validation.py
python3 scripts/her2/fig6_her2_boxplots.py
```

| Script | Output | Manuscript value |
|---|---|---|
| `her2_ttest_top10pctvar.py` | `her2_ttest_results_top10pctvar.csv` | **Primary hit list.** 3,341 testable pairs → **57 HER2-selective pairs** |
| `erbb_pathway_enrichment.py` | *(stdout)* | KEGG ErbB, **OR=5.3, p=0.0044** |
| `erbb2_reactome_enrichment.py` | *(stdout)* | Reactome "Signaling by ERBB2", **OR=14.1, p=2×10⁻⁵** |
| `erbb2_interactor_gene_enrichment_efilter.py` | `erbb2_interactor_gene_enrichment_efilter.csv` | **BioGRID OR=2.07, p=0.045; HINT OR=3.70, p=0.0051** |
| `oesophageal_validation.py` | `oesophageal_validation{,_summary}.csv` | 52/57 (91%) concordant, sign-test p=3×10⁻¹¹; 21 at q<0.20; 36/41 (88%) excluding amplicon pairs, p=3.9×10⁻⁷ |
| `fig6_her2_boxplots.py` | `fig6_her2_boxplots.{png,pdf}` | **Figure 6b** |

## Paths

Each script resolves the repo root from its own location, so it runs from any
working directory:

```python
BASE_DIR = os.environ.get("HER2_ROOT", <two levels above scripts/her2/>)
IN  = os.path.join(BASE_DIR, "data", "input",  "her2")
OUT = os.path.join(BASE_DIR, "data", "output", "her2")
```

Results are written to `data/output/her2/`, the figure to `figures/`. Set
`HER2_ROOT` to override. All three locations are gitignored.

## Method

1. **Input:** `full_SLprediction_matrix.csv` (1,005 cell lines × 33,419 paralog
   pairs). HER2 status from `Consensus HER2 Status`. Overlap = **42 cell lines
   (10 HER2+, 32 HER2−)**.
2. **Variance filter:** top 10% most variable pairs across the 42 lines (3,342;
   3,341 testable with ≥3 values per group).
3. **Welch's two-sided t-test** HER2+ vs HER2−, **BH-FDR** across tested pairs.
4. **Hit definition:** q < 0.20 **and** mean(HER2+) − mean(HER2−) > 0.1 → **57 pairs**.

Enrichments are gene-level Fisher's exact tests: universe = all genes in tested
pairs (**3,050** after removing ERBB2), query = genes in the 57 hit pairs (**80**
after removing ERBB2). ERBB2 is excluded from every gene set and annotation set.
The 17q12 amplicon is a known confounder, so each enrichment is reported both
including and excluding a curated amplicon gene set.

Cross-tissue validation applies the fixed 57-pair list to 42 gastro-oesophageal
lines (6 HER2+, 36 HER2−), with no variance pre-filter.

## Requirements

`pandas`, `numpy`, `scipy`, `openpyxl`; `matplotlib` and `seaborn` for the figure.
Verified against scipy 1.13.1 — `oesophageal_validation.py` supports both
`scipy.stats.binom_test` and its replacement `scipy.stats.binomtest`.
