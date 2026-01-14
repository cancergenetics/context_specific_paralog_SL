import pandas as pd
from pathlib import Path
from tqdm import tqdm

# ---------- CONFIG ----------
input_dir = Path("/home/administrator/SLpred/chunks_go")
output_dir = Path("/home/administrator/SLpred/chunks_processed")
output_dir.mkdir(parents=True, exist_ok=True)

# ---------- Feature Lists ----------
drop_na_values = ['rMaxExp_A1A2', 'rMinExp_A1A2', 'max_ranked_A1A2', 'min_ranked_A1A2']

fillna_zero_cols = ['weighted_PPI_expression', 'smallest_BP_GO_expression', 'go_CC_expression']

fillna_large_cols = ['weighted_PPI_essentiality', 'smallest_BP_GO_essentiality', 'smallest_CC_GO_essentiality']

# ---------- Function ----------

pd.options.mode.copy_on_write = True

def existing_cols(df, cols):
    return [c for c in cols if c in df.columns]

def preprocess_dataset(
    df,
    dropna_cols=None,
    fillna_zero_cols=None,
    fillna_large_cols=None,
    fillna_large_value=18_000
):
    if dropna_cols:
        keep = existing_cols(df, dropna_cols)
        if keep:
            # ensure we own this df after filtering
            df = df.dropna(subset=keep, how='any', axis=0).copy()

    cols0 = existing_cols(df, fillna_zero_cols or [])
    if cols0:
        df.loc[:, cols0] = df.loc[:, cols0].fillna(0)

    colsL = existing_cols(df, fillna_large_cols or [])
    if colsL:
        df.loc[:, colsL] = df.loc[:, colsL].fillna(fillna_large_value)

    return df.reset_index(drop=True)

# ---------- Run Preprocessing ----------
chunk_files = sorted(Path(input_dir).glob("*.parquet"))

for chunk_path in tqdm(chunk_files, desc="Preprocessing", unit="file"):  # tqdm optional
    print(f"Processing {chunk_path.name}")
    df = pd.read_parquet(chunk_path, engine="pyarrow")  # faster IO

    df_cleaned = preprocess_dataset(
        df=df,
        dropna_cols=drop_na_values,
        fillna_zero_cols=fillna_zero_cols,
        fillna_large_cols=fillna_large_cols
    )

    output_path = output_dir / chunk_path.name
    df_cleaned.to_parquet(output_path, engine="pyarrow", compression="zstd", index=False)
    print(f"Saved: {output_path.name} | Rows: {len(df_cleaned)}\n")

print("Preprocessing complete for all chunks.")