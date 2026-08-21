"""
Fig 6-style box plots matching fig6_breast_cancer_map.ipynb (cells 34-37):
predicted SL score in HER2+ vs HER2- breast cancer cell lines for selected
paralog pairs. Replicates the notebook's draw_cont_graph / draw_signif_line
styling (seaborn boxplots, HER2- #56B4E9 / HER2+ #BBBBBB). P-values are from
Welch's two-sided t-test, matching the Figure 6 caption and the discovery
analysis in her2_ttest_top10pctvar.py (the notebook originally used a
Mann-Whitney U test; this was changed to Welch for the revision).
Scores from full_SLprediction_matrix.csv; consensus HER2 labels.
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Repo root, resolved from this file's location (scripts/her2/) so the scripts run
# from any working directory; set HER2_ROOT to override.
BASE_DIR = os.environ.get("HER2_ROOT", os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)))
IN = os.path.join(BASE_DIR, "data", "input", "her2")
OUT = os.path.join(BASE_DIR, "data", "output", "her2")
FIGS = os.path.join(BASE_DIR, "figures")
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIGS, exist_ok=True)
PAIRS = ["EGFR_ERBB2", "AKT1_AKT2", "TRAF3_TRAF4",
         "AP1G1_AP1G2", "SAR1A_SAR1B", "ASF1A_ASF1B"]
FEATURE = "full_prediction_score"

# ---- notebook helpers (verbatim) -------------------------------------------
def draw_signif_line(ax, x0, x1, y, pval, hd=20):
    h = y / hd
    text = '$p$=%.1e' % pval if pval < 0.0001 else '$p$=%.4f' % pval
    ax.plot([x0, x0, x1, x1], [y + h, y + 2*h, y + 2*h, y + h], lw=0.8, c='#666')
    ax.text((x0 + x1) * .5, y + 2.5*h, text, ha='center', va='bottom',
            color='#666', fontsize=9)

def compute_upper_whisker(df, feature):
    Q1 = df[feature].quantile(0.25); Q3 = df[feature].quantile(0.75)
    IQR = Q3 - Q1
    return df[df[feature] <= (Q3 + 1.5*IQR)][feature].max()

def draw_cont_graph(df, feature, ax, label=None, draw_signif=True, hd=20, xticklabels=None):
    my_pal = {1: "#BBBBBB", 0: "#56B4E9"}
    mutated = (df['mut_binary'] == 1).sum(); not_mutated = (df['mut_binary'] == 0).sum()
    sns.despine(top=True, right=True, left=False, bottom=False)
    sns.boxplot(y=feature, x='mut_binary', hue='mut_binary', data=df, ax=ax,
                linewidth=0.8, saturation=0.8, showfliers=False, showmeans=False,
                palette=my_pal, order=[0, 1], medianprops={'color': 'black', 'linewidth': 0.8},
                boxprops={'edgecolor': 'black', 'alpha': 0.85})
    ax.set_ylabel(label if label else feature, fontsize=12); ax.set_xlabel('')
    ax.set_xticks([0, 1])
    if xticklabels is None:
        xticklabels = [f'WT \n(n={not_mutated})', f'M \n(n={mutated})']
    ax.set_xticklabels(xticklabels, fontsize=12)
    ax.tick_params(axis='y', rotation=0, labelsize=8)
    leg = ax.legend()
    if leg: leg.remove()
    if draw_signif:
        pval = stats.ttest_ind(df.loc[df['mut_binary'] == 1, feature],
                               df.loc[df['mut_binary'] == 0, feature],
                               equal_var=False)[1]   # Welch t-test
        uw = max(compute_upper_whisker(df.loc[df['mut_binary'] == 1], feature),
                 compute_upper_whisker(df.loc[df['mut_binary'] == 0], feature))
        draw_signif_line(ax, 0, 1, uw + 0.05, pval, hd)

# ---- data ------------------------------------------------------------------
her2 = pd.read_excel(f"{IN}/DepmapHer2_withDai.xlsx").set_index("Model")["Consensus HER2 Status"]
mut = her2.map({"Positive": 1, "Negative": 0}).dropna()
mat = pd.read_csv(f"{IN}/full_SLprediction_matrix.csv", index_col="DepMap_ID",
                  usecols=["DepMap_ID"] + PAIRS)
common = mat.index.intersection(mut.index)
mat = mat.loc[common]; mut = mut.loc[common]

# ---- figure (mirrors cell 37) ----------------------------------------------
fig, axes = plt.subplots(3, 2, figsize=(4.5, 9), sharey=True, sharex=True)
for i, pair in enumerate(PAIRS):
    ax = axes[i // 2, i % 2]
    df = pd.DataFrame({FEATURE: mat[pair].values, "mut_binary": mut.values}).dropna()
    df["mut_binary"] = df["mut_binary"].astype(int)
    n_neg = (df["mut_binary"] == 0).sum(); n_pos = (df["mut_binary"] == 1).sum()
    draw_cont_graph(df, FEATURE, ax, label='Full Prediction Score',
                    xticklabels=[f'HER2-\n(n={n_neg})', f'HER2+\n(n={n_pos})'],
                    draw_signif=True, hd=22)
    ax.set_title(pair, fontsize=10, fontweight='bold', pad=15)

plt.tight_layout(rect=[0, 0.03, 1, 0.97])
fig.savefig(f"{FIGS}/fig6_her2_boxplots.png", dpi=500, bbox_inches="tight")
fig.savefig(f"{FIGS}/fig6_her2_boxplots.pdf", bbox_inches="tight")
print("Saved fig6_her2_boxplots.png (+ .pdf)")
for pair in PAIRS:
    df = pd.DataFrame({FEATURE: mat[pair].values, "mut_binary": mut.values}).dropna()
    p = stats.ttest_ind(df.loc[df.mut_binary == 1, FEATURE],
                        df.loc[df.mut_binary == 0, FEATURE], equal_var=False)[1]
    print(f"  {pair:14s} Welch p={p:.3e}")
