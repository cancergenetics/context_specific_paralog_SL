import pandas as pd
import pickle
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from tqdm import tqdm

# --------------------- CONFIGURATION ---------------------

# Path to your trained RF model
model_path = Path("/home/administrator/SLpred/model/contextualized_model.pickle")

# Directory with input files (annotated .parquet chunks)
input_dir = Path("/home/administrator/SLpred/chunks_processed")

# Where to save the output predictions
output_dir = Path("/home/administrator/SLpred/output/early")
output_dir.mkdir(parents=True, exist_ok=True)

matrix_out = output_dir / "contextualized_prediction_matrix.parquet"

# Features used by the model — must match training
feature_columns = [
    'rMaxExp_A1A2', 'rMinExp_A1A2',
    'max_ranked_A1A2', 'min_ranked_A1A2',
    'max_cn', 'min_cn', 'Protein_Altering', 'Damaging', 
    'min_sequence_identity',
    'prediction_score', 
    'weighted_PPI_essentiality', 'weighted_PPI_expression',
    'go_CC_expression', 'smallest_BP_GO_essentiality', 'smallest_BP_GO_expression',
    'smallest_CC_GO_essentiality'
]

pd.options.mode.copy_on_write = True

# --------------------- MAIN SCRIPT ---------------------

# Load model
with open(model_path, "rb") as f:
    RF = pickle.load(f)
print("Model loaded.")

# Use model's feature order if available
try:
    feature_columns = list(RF.feature_names_in_)
except AttributeError:
    raise RuntimeError("Model was not trained using a DataFrame — feature names are unavailable.")

def merge_into_matrix(M: pd.DataFrame | None, P: pd.DataFrame) -> pd.DataFrame:
    """Expand M to cover P's rows/cols and write only P's observed block."""
    if M is None:
        return P.copy()
    M = M.reindex(index=M.index.union(P.index), columns=M.columns.union(P.columns))
    # DO NOT reindex P; assign only the concrete block
    M.loc[P.index, P.columns] = P.to_numpy()
    return M

M = None  # final matrix DepMap_ID × genepair

for chunk_file in tqdm(sorted(input_dir.glob("*.parquet")), desc="Predicting"):
    # Read only what we need
    needed = list(dict.fromkeys(["DepMap_ID", "genepair"] + feature_columns))
    df = pd.read_parquet(chunk_file, engine="pyarrow", columns=needed)

    # quick checks
    if "DepMap_ID" not in df.columns or "genepair" not in df.columns:
        print(f"Skipping {chunk_file.name} — needs DepMap_ID & genepair")
        continue
    missing = [c for c in feature_columns if c not in df.columns]
    if missing:
        print(f"Skipping {chunk_file.name} — missing features: {missing}")
        continue

    # predict
    X = df.loc[:, feature_columns].astype("float32", copy=False)
    ypred = RF.predict_proba(X)[:, 1].astype("float32")

    # pivot this chunk to DepMap_ID × genepair
    small = df.loc[:, ["DepMap_ID", "genepair"]].copy()
    small["prediction"] = ypred
    P = small.pivot_table(index="DepMap_ID",
                          columns="genepair",
                          values="prediction",
                          aggfunc="mean").astype("float32")

    # merge chunk into master matrix
    M = merge_into_matrix(M, P)

# save once
if M is None:
    raise RuntimeError("No predictions produced.")
M.sort_index().to_parquet(matrix_out, engine="pyarrow", index=True)
print(f"Saved matrix to {matrix_out} with shape {M.shape}")

# Optional sanity check
non_null = int(M.notna().to_numpy().sum())
total = M.shape[0] * M.shape[1]
print(f"Non-null cells: {non_null:,} / {total:,} ({non_null/total:.2%})")