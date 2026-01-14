import os
import pandas as pd
import pyarrow.parquet as pq
import numpy as np
from pathlib import Path

# Define paths using Pathlib
ranked_essentiality_file = Path("/home/administrator/SLpred/ranked/ranked_essentiality.csv")
input_chunks_dir = Path("/home/administrator/SLpred/chunks")  # Ensure this is a Path object
output_chunks_dir = Path("/home/administrator/SLpred/chunks_ranked")

output_chunks_dir.mkdir(parents=True, exist_ok=True)  # Ensure output directory exists

# Load supporting files
cell_lines = pd.read_csv("/home/administrator/SLpred/DepMap_ID.csv").rename(columns={"0": "DepMap_ID"})

############################################
def process_ranked_ess_files(ranked_essentiality_file, cell_lines):

    ranked_df = pd.read_csv(ranked_essentiality_file, index_col=0)
    filtered_ranked_df = ranked_df[ranked_df.index.isin(cell_lines['DepMap_ID'])]
    filtered_ranked_df = filtered_ranked_df.reset_index(drop=False).rename(columns={"index":"DepMap_ID"})

    # Convert wide format to long format
    melt_ranked_df = pd.melt(filtered_ranked_df, id_vars=["DepMap_ID"], var_name="entrez_id", value_name="ranked_essentiality")
    melt_ranked_df['entrez_id'] = melt_ranked_df['entrez_id'].astype(int).astype(str)

    print(f"Processed {ranked_essentiality_file}")
    return melt_ranked_df

melt_ranked_df = process_ranked_ess_files(ranked_essentiality_file, cell_lines)

############################################

def process_chunk(chunk_path, melt_ranked_df):
    """Reads a chunk, annotates features, and saves it."""
    df = pd.read_parquet(chunk_path)
    print(f"Processing {chunk_path}...")

    # Merge all PPI data at once
    df = pd.merge(df, melt_ranked_df.rename(columns={'entrez_id':'A1_entrez', 'ranked_essentiality':'A1_rank'}), on=['DepMap_ID', 'A1_entrez'], how='left')
    df_v2 = pd.merge(df, melt_ranked_df.rename(columns={'entrez_id':'A2_entrez', 'ranked_essentiality':'A2_rank'}), on=['DepMap_ID', 'A2_entrez'], how='left')

    df_v2['max_ranked_A1A2'] = df_v2[['A1_rank', 'A2_rank']].max(axis=1)
    df_v2['min_ranked_A1A2'] = df_v2[['A1_rank', 'A2_rank']].min(axis=1)

    # Save the processed chunk
    output_path = output_chunks_dir / chunk_path.name
    df_v2.to_parquet(output_path)
    print(f"Saved to {output_path}")

# Process all chunks using pathlib
chunk_files = sorted(input_chunks_dir.iterdir())

for chunk_file in chunk_files:
    if chunk_file.suffix == ".parquet":  # Only process Parquet files
        process_chunk(chunk_file, melt_ranked_df)

print("Annotation completed for all chunks.")