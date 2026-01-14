import pandas as pd
import numpy as np
import pickle
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from tqdm import tqdm

# --------------------- CONFIGURATION ---------------------
model_path = Path("/home/administrator/SLpred/model/full_model.pickle")
input_dir  = Path("/home/administrator/SLpred/chunks_processed")
output_dir = Path("/home/administrator/SLpred/output/late")
output_dir.mkdir(parents=True, exist_ok=True)

# final matrix path
matrix_out = output_dir / "full_prediction_matrix.parquet"

pd.options.mode.copy_on_write = True

# --------------------- LOAD MODEL & FEATURES ---------------------
with open(model_path, "rb") as f:
    RF: RandomForestClassifier = pickle.load(f)
print("Model loaded.")

# use model’s training feature order if available, else your list
feature_columns = list(getattr(RF, "feature_names_in_", [
    'rMaxExp_A1A2', 'rMinExp_A1A2',
    'max_ranked_A1A2', 'min_ranked_A1A2',
    'max_cn', 'min_cn', 'Protein_Altering', 'Damaging', 
    'min_sequence_identity',
    'ranked_Essentiality_weighted_PPI', 'Expression_weighted_PPI',
    'smallest_GO_ranked_ess', 'smallest_GO_CC_ranked_ess',
    'smallest_gene_expression', 'smallest_GO_CC_expression',
    'closest', 'WGD', 'family_size',
    'cds_length_ratio', 'shared_domains', 'has_pombe_ortholog',
    'has_essential_pombe_ortholog', 'has_cerevisiae_ortholog',
    'has_essential_cerevisiae_ortholog', 'conservation_score', 'mean_age',
    'either_in_complex', 'mean_complex_essentiality', 'colocalisation',
    'interact', 'n_total_ppi', 'fet_ppi_overlap',
    'gtex_spearman_corr', 'gtex_min_mean_expr', 'gtex_max_mean_expr'
]))

# --------------------- HELPERS ---------------------
def merge_into_matrix(M: pd.DataFrame | None, P: pd.DataFrame) -> pd.DataFrame:
    """Expand M to cover P's rows/cols and overwrite only P's block."""
    if M is None:
        return P.copy()
    M = M.reindex(index=M.index.union(P.index), columns=M.columns.union(P.columns))
    # write ONLY the observed block from P (no reindex on P!)
    M.loc[P.index, P.columns] = P.to_numpy()
    return M

# --------------------- MAIN ---------------------
M = None  # DepMap_ID × genepair matrix

for chunk_file in tqdm(sorted(input_dir.glob("*.parquet")), desc="Predicting (late)"):
    needed = list(dict.fromkeys(["DepMap_ID", "genepair"] + feature_columns))
    df = pd.read_parquet(chunk_file, engine="pyarrow", columns=needed)

    if "DepMap_ID" not in df.columns or "genepair" not in df.columns:
        print(f"Skipping {chunk_file.name} — needs DepMap_ID & genepair")
        continue
    missing = [c for c in feature_columns if c not in df.columns]
    if missing:
        print(f"Skipping {chunk_file.name} — missing features: {missing}")
        continue

    X = df.loc[:, feature_columns].astype("float32", copy=False)
    ypred = RF.predict_proba(X)[:, 1].astype("float32")

    small = df.loc[:, ["DepMap_ID", "genepair"]].copy()
    small["prediction_late"] = ypred

    # tolerate any accidental duplicates
    P = small.pivot_table(
        index="DepMap_ID",
        columns="genepair",
        values="prediction_late",
        aggfunc="mean"
    ).astype("float32")

    M = merge_into_matrix(M, P)

# save once
if M is None:
    raise RuntimeError("No predictions produced. Check input chunks/features.")

M.sort_index().to_parquet(matrix_out, engine="pyarrow", index=True)
print(f"Saved matrix to {matrix_out} with shape {M.shape}")

# quick sanity report (optional)
non_null = int(M.notna().to_numpy().sum())
total = M.shape[0] * M.shape[1]
print(f"Non-null cells: {non_null:,} / {total:,} ({non_null/total:.2%})")