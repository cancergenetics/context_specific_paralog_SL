"""
Test the 57 HER2-selective breast paralog pairs (q<0.20 & mean_diff>0.1, higher
in HER2+) in an independent OESOPHAGEAL/gastro-oesophageal cohort, stratified by
CN-based HER2 status. Same Welch t-test on the prediction scores.

Output: oesophageal_validation.csv
"""
import os
import numpy as np, pandas as pd
from scipy import stats

# scipy >= 1.12 removed stats.binom_test in favour of stats.binomtest; support both
if hasattr(stats, "binom_test"):
    def sign_test_p(k, n):
        return stats.binom_test(k, n, 0.5, alternative="greater")
else:
    def sign_test_p(k, n):
        return stats.binomtest(k, n, 0.5, alternative="greater").pvalue

# Repo root, resolved from this file's location (scripts/her2/) so the scripts run
# from any working directory; set HER2_ROOT to override.
BASE_DIR = os.environ.get("HER2_ROOT", os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)))
IN = os.path.join(BASE_DIR, "data", "input", "her2")
OUT = os.path.join(BASE_DIR, "data", "output", "her2")
os.makedirs(OUT, exist_ok=True)
AMPLICON = {"ERBB2","GRB7","STARD3","PNMT","NEUROD2","MIEN1","C17orf37","PGAP3",
 "TCAP","PPP1R1B","GSDMA","GSDMB","ORMDL3","IKZF3","ZPBP2","PSMD3","MED1",
 "CDK12","RARA","TOP2A","NR1D1","THRA","MSL1","CASC3","WIPF2","CDC6"}

# ---- breast hit pairs (HER2-selective) -------------------------------------
br = pd.read_csv(f"{OUT}/her2_ttest_results_top10pctvar.csv")
br = br[br["p_value"].notna()]
hits = br[(br["q_value_BH"] < 0.20) & (br["mean_diff"] > 0.1)].copy()
pairs = hits["paralog_pair"].tolist()
print(f"Breast HER2-selective pairs to test: {len(pairs)}")

# ---- oesophageal HER2 status -----------------------------------------------
oe = pd.read_excel(f"{IN}/Oesophageal_HER2_CN.xlsx")
status = oe.set_index("Model")["HER2 Status (CN-based)"]
status = status[status.isin(["Positive", "Negative"])]

# ---- prediction matrix -----------------------------------------------------
mat = pd.read_csv(f"{IN}/full_SLprediction_matrix.csv", index_col="DepMap_ID",
                  usecols=["DepMap_ID"] + pairs)
common = mat.index.intersection(status.index)
mat = mat.loc[common]; status = status.loc[common]
pos_id = status.index[status == "Positive"]; neg_id = status.index[status == "Negative"]
print(f"Oesophageal overlap with prediction matrix: {len(common)} cell lines "
      f"({len(pos_id)} HER2+, {len(neg_id)} HER2-)")

# ---- per-pair Welch t-test in oesophageal ----------------------------------
rows = []
for pair in pairs:
    a = mat.loc[pos_id, pair].dropna().to_numpy()
    b = mat.loc[neg_id, pair].dropna().to_numpy()
    if len(a) < 3 or len(b) < 3:
        rows.append({"paralog_pair": pair, "n_pos": len(a), "n_neg": len(b),
                     "oe_mean_diff": np.nan, "oe_p": np.nan}); continue
    t, p = stats.ttest_ind(a, b, equal_var=False)
    rows.append({"paralog_pair": pair, "n_pos": len(a), "n_neg": len(b),
                 "oe_mean_HER2pos": round(a.mean(), 4), "oe_mean_HER2neg": round(b.mean(), 4),
                 "oe_mean_diff": round(a.mean() - b.mean(), 4), "oe_p": p})
oe_df = pd.DataFrame(rows)

# BH-FDR over tested breast hits
valid = oe_df["oe_p"].notna().to_numpy()
pv = oe_df.loc[valid, "oe_p"].to_numpy(); m = len(pv)
order = np.argsort(pv); q = np.empty(m)
q[order] = np.clip(np.minimum.accumulate((pv[order]*m/np.arange(1, m+1))[::-1])[::-1], 0, 1)
oe_df.loc[valid, "oe_q_BH"] = q

# merge breast stats for comparison
out = hits[["paralog_pair", "mean_diff", "q_value_BH"]].rename(
    columns={"mean_diff": "breast_mean_diff", "q_value_BH": "breast_q"}).merge(oe_df, on="paralog_pair")
out["oe_higher_in_HER2pos"] = out["oe_mean_diff"] > 0
out["has_ERBB2"] = out["paralog_pair"].apply(lambda p: "ERBB2" in p.split("_"))
out["has_amplicon"] = out["paralog_pair"].apply(lambda p: any(g in AMPLICON for g in p.split("_")))
out = out.sort_values("oe_p").reset_index(drop=True)
out.to_csv(f"{OUT}/oesophageal_validation.csv", index=False)

# ---- amplicon-stratified concordance summary -------------------------------
def strat(df, label):
    d = df[df["oe_p"].notna()]
    n = len(d); conc = int((d["oe_mean_diff"] > 0).sum())
    return {"subset": label, "n_pairs": n, "concordant_HER2pos": conc,
            "pct_concordant": round(100*conc/n, 1),
            "sign_test_p": sign_test_p(conc, n),
            "median_oe_mean_diff": round(d["oe_mean_diff"].median(), 4),
            "wilcoxon_p_vs0": stats.wilcoxon(d["oe_mean_diff"])[1],
            "n_raw_p_lt_0.05": int((d["oe_p"] < 0.05).sum())}
summary = pd.DataFrame([
    strat(out, "all_breast_hits"),
    strat(out[~out["has_ERBB2"]], "excluding_ERBB2_pairs"),
    strat(out[~out["has_amplicon"]], "excluding_any_amplicon_pairs"),
])
summary.to_csv(f"{OUT}/oesophageal_validation_summary.csv", index=False)
print("\nConcordance summary:")
print(summary.to_string(index=False))

# ---- summary ---------------------------------------------------------------
tested = out[out["oe_p"].notna()]
n = len(tested)
conc = int((tested["oe_mean_diff"] > 0).sum())
print(f"\nTested in oesophageal: {n} / {len(pairs)} pairs")
print(f"Direction concordant (higher in oeso HER2+): {conc}/{n} ({100*conc/n:.0f}%)")
# sign test vs 50/50
sign_p = sign_test_p(conc, n)
print(f"  sign test (concordance > 50%): p={sign_p:.2g}")
print(f"raw p<0.05 in oesophageal: {int((tested['oe_p']<0.05).sum())}")
print(f"  of which higher in HER2+: {int(((tested['oe_p']<0.05)&(tested['oe_mean_diff']>0)).sum())}")
print(f"oe q<0.20: {int((tested['oe_q_BH']<0.20).sum())}")
med = tested['oe_mean_diff'].median()
print(f"median oesophageal mean_diff across breast hits: {med:+.4f} "
      f"(Wilcoxon vs 0 p={stats.wilcoxon(tested['oe_mean_diff'])[1]:.2g})")
print("\nTop 20 by oesophageal p-value:")
cols = ["paralog_pair", "breast_mean_diff", "oe_mean_diff", "oe_p", "oe_q_BH"]
with pd.option_context("display.width", 200):
    print(out[cols].head(20).to_string(index=False))
print(f"\n-> {OUT}/oesophageal_validation.csv")
