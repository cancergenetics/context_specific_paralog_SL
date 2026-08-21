"""
Welch's t-test of paralog SL predictions between HER2+ and HER2- cell lines,
restricted to the 10% MOST VARIABLE paralog pairs.

Mirrors her2_mwu_top10pctvar.py but uses Welch's two-sample t-test
(unequal variances) instead of the Mann-Whitney U test. Welch is used
because the groups are small and unequal in size (10 HER2+ vs 32 HER2-).

Variance for the top-10% selection is computed over the 42 overlapping cell
lines (HER2+ and HER2- combined), exactly as in the Mann-Whitney version, so
the same set of pairs is tested.

Output:
  - her2_ttest_results_top10pctvar.csv
"""

import os
import numpy as np
import pandas as pd
from scipy import stats

# Repo root, resolved from this file's location (scripts/her2/) so the scripts run
# from any working directory; set HER2_ROOT to override.
BASE_DIR = os.environ.get("HER2_ROOT", os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)))
IN = os.path.join(BASE_DIR, "data", "input", "her2")
OUT = os.path.join(BASE_DIR, "data", "output", "her2")
os.makedirs(OUT, exist_ok=True)
MIN_PER_GROUP = 3       # require >=3 non-NaN predictions in each group to test
TOP_VAR_FRAC = 0.10     # keep the most variable 10% of pairs

# ---- 1. HER2 status --------------------------------------------------------
her2 = pd.read_excel(f"{IN}/DepmapHer2_withDai.xlsx")
status = her2.set_index("Model")["Consensus HER2 Status"]
status = status[status.isin(["Positive", "Negative"])]
print(f"HER2 calls available: {status.value_counts().to_dict()}")

# ---- 2. Prediction matrix --------------------------------------------------
mat = pd.read_csv(f"{IN}/full_SLprediction_matrix.csv", index_col="DepMap_ID")
print(f"Prediction matrix: {mat.shape[0]} cell lines x {mat.shape[1]} pairs")

# ---- 3. Align cell lines ---------------------------------------------------
common = mat.index.intersection(status.index)
mat = mat.loc[common]
status = status.loc[common]
pos_ids = status.index[status == "Positive"]
neg_ids = status.index[status == "Negative"]
print(f"Overlap: {len(common)} cell lines "
      f"({len(pos_ids)} HER2+, {len(neg_ids)} HER2-)")

# ---- 4. Variance filter over the 42 overlapping cell lines -----------------
arr = mat.to_numpy(dtype=float)
n_obs = np.sum(~np.isnan(arr), axis=0)
with np.errstate(invalid="ignore"):
    variance = np.nanvar(arr, axis=0, ddof=1)
rank_key = np.where(n_obs >= 2, variance, -np.inf)
n_keep = int(np.ceil(TOP_VAR_FRAC * mat.shape[1]))
keep_idx = np.sort(np.argsort(rank_key)[::-1][:n_keep])
var_cut = variance[keep_idx].min()
print(f"Variance filter: keeping top {TOP_VAR_FRAC:.0%} = {n_keep} pairs "
      f"(variance >= {var_cut:.3g})")

mat = mat.iloc[:, keep_idx]
pairs = mat.columns.to_numpy()
kept_var = variance[keep_idx]
P = mat.loc[pos_ids].to_numpy(dtype=float)
N = mat.loc[neg_ids].to_numpy(dtype=float)

# ---- 5. Per-pair Welch t-test ---------------------------------------------
n_pairs = P.shape[1]
n_pos = np.zeros(n_pairs, dtype=int)
n_neg = np.zeros(n_pairs, dtype=int)
mean_pos = np.full(n_pairs, np.nan)
mean_neg = np.full(n_pairs, np.nan)
t_stat = np.full(n_pairs, np.nan)
p = np.full(n_pairs, np.nan)
cohen_d = np.full(n_pairs, np.nan)   # pooled-SD standardized mean difference

for j in range(n_pairs):
    a = P[:, j]; b = N[:, j]
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    na, nb = len(a), len(b)
    n_pos[j], n_neg[j] = na, nb
    if na:
        mean_pos[j] = a.mean()
    if nb:
        mean_neg[j] = b.mean()
    if na < MIN_PER_GROUP or nb < MIN_PER_GROUP:
        continue
    va, vb = a.var(ddof=1), b.var(ddof=1)
    if va == 0 and vb == 0:
        continue
    t, pv = stats.ttest_ind(a, b, equal_var=False)   # Welch
    t_stat[j] = t
    p[j] = pv
    # Cohen's d with pooled SD (>0 => higher in HER2+)
    sp = np.sqrt(((na - 1) * va + (nb - 1) * vb) / (na + nb - 2))
    cohen_d[j] = (a.mean() - b.mean()) / sp if sp > 0 else np.nan

out = pd.DataFrame({
    "paralog_pair": pairs,
    "variance_42cl": kept_var,
    "n_HER2pos": n_pos,
    "n_HER2neg": n_neg,
    "mean_HER2pos": mean_pos,
    "mean_HER2neg": mean_neg,
    "mean_diff": mean_pos - mean_neg,
    "t_statistic": t_stat,
    "p_value": p,
    "cohen_d": cohen_d,
})
out["direction"] = np.where(out["mean_diff"] >= 0,
                            "higher_in_HER2pos", "higher_in_HER2neg")

# ---- 6. BH-FDR over the tested (top-variance) pairs ------------------------
out["q_value_BH"] = np.nan
valid = out["p_value"].notna().to_numpy()
pv = out.loc[valid, "p_value"].to_numpy()
m = len(pv)
order = np.argsort(pv)
ranked = pv[order]
q_ranked = ranked * m / np.arange(1, m + 1)
q_ranked = np.minimum.accumulate(q_ranked[::-1])[::-1]
q = np.empty(m)
q[order] = np.clip(q_ranked, 0, 1)
out.loc[valid, "q_value_BH"] = q

out = out.sort_values("p_value", kind="mergesort").reset_index(drop=True)
out.to_csv(f"{OUT}/her2_ttest_results_top10pctvar.csv", index=False)

# ---- 7. Summary ------------------------------------------------------------
print(f"\nPairs tested (>= {MIN_PER_GROUP}/group, non-degenerate): "
      f"{m} / {len(out)}")
for thr in (0.05, 0.10, 0.25):
    print(f"  q (BH) < {thr:>4}: {int((out['q_value_BH'] < thr).sum())} pairs")
print(f"  raw p < 0.05 : {int((out['p_value'] < 0.05).sum())} pairs")
print("\nTop 20 by p-value:")
cols = ["paralog_pair", "variance_42cl", "mean_HER2pos", "mean_HER2neg",
        "mean_diff", "t_statistic", "p_value", "q_value_BH", "cohen_d",
        "direction"]
with pd.option_context("display.width", 240, "display.max_columns", None):
    print(out[cols].head(20).to_string(index=False))
print(f"\nFull results written to: {OUT}/her2_ttest_results_top10pctvar.csv")
