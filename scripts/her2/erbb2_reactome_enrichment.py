"""
Pathway enrichment: are SL hit genes enriched for ERBB/ErbB pathway members
(KEGG hsa04012) relative to the tested background?

Universe   = all genes in tested pairs (top-10% var, consensus, Welch t-test).
Hit sets   = (a) q<0.20 ; (b) q<0.20 & mean_diff>0.1 (HER2+, effect-size filter).
ERBB2 excluded from gene sets (comparability with the interactor analysis).
Reported incl. and excl. amplicon genes. Fisher exact (hit vs non-hit tested).
"""
import os
import pandas as pd
from scipy import stats

# Repo root, resolved from this file's location (scripts/her2/) so the scripts run
# from any working directory; set HER2_ROOT to override.
BASE_DIR = os.environ.get("HER2_ROOT", os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)))
IN = os.path.join(BASE_DIR, "data", "input", "her2")
OUT = os.path.join(BASE_DIR, "data", "output", "her2")
os.makedirs(OUT, exist_ok=True)
AMPLICON = {"ERBB2","GRB7","STARD3","PNMT","NEUROD2","MIEN1","C17orf37","PGAP3",
            "TCAP","PPP1R1B","GSDMA","GSDMB","ORMDL3","IKZF3","ZPBP2","PSMD3",
            "MED1","CDK12","RARA","TOP2A","NR1D1","THRA","MSL1","CASC3","WIPF2","CDC6"}

pathway = {l.strip() for l in open(f"{IN}/reactome_erbb2.txt") if l.strip()} - {"ERBB2"}
print(f"Reactome Signaling by ERBB2 genes (excl ERBB2): {len(pathway)}")

sl = pd.read_csv(f"{OUT}/her2_ttest_results_top10pctvar.csv")
tdf = sl[sl["p_value"].notna()]
tested = {g for p in tdf["paralog_pair"] for g in p.split("_")} - {"ERBB2"}

hitsets = {
    "q<0.20": {g for p in tdf[tdf["q_value_BH"] < 0.20]["paralog_pair"]
               for g in p.split("_")} - {"ERBB2"},
    "q<0.20 & mean_diff>0.1": {g for p in tdf[(tdf["q_value_BH"] < 0.20) &
               (tdf["mean_diff"] > 0.1)]["paralog_pair"]
               for g in p.split("_")} - {"ERBB2"},
}

def fisher(hits, genes_sub, label):
    pw = pathway & genes_sub
    hit = hits & genes_sub
    nonhit = genes_sub - hits
    a = len(hit & pw); b = len(hit) - a
    c = len(nonhit & pw); d = len(nonhit) - c
    OR, p = stats.fisher_exact([[a, b], [c, d]])
    print(f"   {label:<16} hit {a}/{len(hit)} ({100*a/len(hit) if hit else 0:.1f}%) "
          f"vs non-hit {c}/{len(nonhit)} ({100*c/len(nonhit) if nonhit else 0:.1f}%)  "
          f"OR={OR:.2f}, p={p:.3g}")
    return hit & pw

for name, hits in hitsets.items():
    print(f"\n=== Hit set: {name}  ({len(hits)} genes) ===")
    pw_hits = fisher(hits, tested, "incl. amplicon")
    fisher(hits, tested - AMPLICON, "excl. amplicon")
    print(f"   pathway members among hits: {sorted(pw_hits)}")

# detail table for the effect-size-filtered hits
hits = hitsets["q<0.20 & mean_diff>0.1"]
det = sorted(hits & pathway)
print(f"\nERBB-pathway hit genes (effect-size-filtered, excl ERBB2): "
      f"{len(det)}\n  {det}")
