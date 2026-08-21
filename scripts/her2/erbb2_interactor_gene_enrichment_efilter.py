"""
Gene-level test: are ERBB2 interactors over-represented among SL hit genes
relative to all tested genes?

Universe   = all genes in tested pairs (top-10% var, consensus, Welch t-test).
Hit genes  = genes in significant pairs (q<0.20 AND mean prediction-score
             difference > 0.1, increased in HER2+: mean_HER2pos-mean_HER2neg>0.1).
ERBB2 is excluded from gene sets and from the interactor sets.
Interactomes: BioGRID (human physical) and HINT. Membership matched by symbol.
Enrichment = Fisher exact (hit vs non-hit tested genes x interactor vs not).
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
# ---- gene sets -------------------------------------------------------------
sl = pd.read_csv(f"{OUT}/her2_ttest_results_top10pctvar.csv")
tested_df = sl[sl["p_value"].notna()]
tested = {g for p in tested_df["paralog_pair"] for g in p.split("_")}
hit_pairs = tested_df[(tested_df["q_value_BH"] < 0.20) &
                      (tested_df["mean_diff"] > 0.1)]   # increased in HER2+
hits = {g for p in hit_pairs["paralog_pair"] for g in p.split("_")}
tested.discard("ERBB2"); hits.discard("ERBB2")
print(f"Hit pairs (q<0.20 & mean_diff>0.1, HER2+): {len(hit_pairs)}")
print(f"Tested genes (excl ERBB2): {len(tested)} | hit genes: {len(hits)}")

# ---- ERBB2 interactor sets (symbols, excl ERBB2) ---------------------------
def biogrid_partners():
    bg = pd.read_csv(f"{IN}/BIOGRID-GENE-108376-5.0.258.tab3.txt", sep="\t", low_memory=False)
    bg = bg[(bg["Organism ID Interactor A"] == 9606) &
            (bg["Organism ID Interactor B"] == 9606) &
            (bg["Experimental System Type"] == "physical")]
    s = set()
    for a, b in zip(bg["Official Symbol Interactor A"], bg["Official Symbol Interactor B"]):
        if a == "ERBB2" and b != "ERBB2": s.add(b)
        elif b == "ERBB2" and a != "ERBB2": s.add(a)
    return s

def hint_partners():
    h = pd.read_csv(f"{IN}/HINT-search-results-main.tsv", sep="\t", header=None,
                    names=["upA","upB","symA","symB","evidence","org","x","y"])
    s = set()
    for a, b in zip(h["symA"], h["symB"]):
        if a == "ERBB2" and b != "ERBB2": s.add(b)
        elif b == "ERBB2" and a != "ERBB2": s.add(a)
    return s

def run(partners, name):
    partners = partners - {"ERBB2"}
    t_in = tested & partners
    h_in = hits & partners
    nonhit = tested - hits
    a = len(h_in); b = len(hits) - a
    c = len(nonhit & partners); d = len(nonhit) - c
    OR, p = stats.fisher_exact([[a, b], [c, d]])
    print(f"\n=== {name} ({len(partners)} ERBB2 partners) ===")
    print(f"  tested genes that are ERBB2 interactors : {len(t_in)}/{len(tested)} "
          f"({100*len(t_in)/len(tested):.1f}%)")
    print(f"  hit genes that are ERBB2 interactors    : {a}/{len(hits)} "
          f"({100*a/len(hits):.1f}%)")
    print(f"  non-hit tested interactors              : {c}/{len(nonhit)} "
          f"({100*c/len(nonhit):.1f}%)")
    print(f"  Fisher OR={OR:.2f}, p={p:.3g}")
    print(f"  hit interactor genes: {sorted(h_in)}")
    return name, len(t_in), a, OR, p

bg = biogrid_partners(); hint = hint_partners()
rows = [run(bg, "BioGRID physical"), run(hint, "HINT")]
print(f"\nBioGRID∩HINT partners: {len(bg & hint)} | union: {len(bg | hint)}")
pd.DataFrame(rows, columns=["interactome", "tested_interactors",
                            "hit_interactors", "odds_ratio", "fisher_p"]
             ).to_csv(f"{OUT}/erbb2_interactor_gene_enrichment_efilter.csv", index=False)
print(f"-> {OUT}/erbb2_interactor_gene_enrichment_efilter.csv")
