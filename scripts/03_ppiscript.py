import os
import pandas as pd
import numpy as np
from pathlib import Path

# === Path setup ===
ppi_folder = Path("/home/administrator/SLpred/PPI/")
input_chunks_dir = Path("/home/administrator/SLpred/chunks_ranked")
output_chunks_dir = Path("/home/administrator/SLpred/chunks_ppi")
output_chunks_dir.mkdir(parents=True, exist_ok=True)

print("Loading DepMap ID list...")
cell_lines = pd.read_csv("/home/administrator/SLpred/DepMap_ID.csv").rename(columns={"0": "DepMap_ID"})
print("Loaded", len(cell_lines), "DepMap IDs.\n")

# === Step 1: Load and melt PPI files ===
def process_ppi_files(ppi_files, cell_lines):
    results = []
    for file in ppi_files:
        print("Reading PPI file:", file.name)
        ppi_df = pd.read_parquet(file)

        # Convert index (DepMap_ID) into a column
        ppi_df["DepMap_ID"] = ppi_df.index
        ppi_df = ppi_df.reset_index(drop=True)

        # Filter to valid DepMap_IDs
        ppi_df = ppi_df[ppi_df["DepMap_ID"].isin(cell_lines["DepMap_ID"])]

        # Melt to long format
        melt_ppi = ppi_df.melt(id_vars=["DepMap_ID"], var_name="genepair", value_name="value")
        melt_ppi = melt_ppi.rename(columns={"value": file.stem})

        print("Finished processing", file.name, "with", len(melt_ppi), "rows.\n")
        results.append(melt_ppi)

    return results

# === Load PPI files ===
ppi_files = sorted(ppi_folder.glob("*.parquet"))
print("Found", len(ppi_files), "PPI files.\n")

results = process_ppi_files(ppi_files, cell_lines)

print("Merging PPI data...")
merged_ppi = pd.merge(results[0], results[1], on=["DepMap_ID", "genepair"], how="left")
print("Merged PPI shape:", merged_ppi.shape, "\n")
print('')
merged_ppi.head()
print('')

# Rename columns for clarity
rename_dict = {
    'combined_weighted_PPI_expression_new': 'Expression_weighted_PPI',
    'combined_weighted_PPI_essentiality_new': 'ranked_Essentiality_weighted_PPI',
}
merged_ppi = merged_ppi.rename(columns=rename_dict)

# === Step 2: Annotate each chunk ===
def process_chunk(chunk_path, merged_ppi):
    print("Processing chunk:", chunk_path.name)
    df = pd.read_parquet(chunk_path)
    before = len(df)

    df = pd.merge(df, merged_ppi, on=['DepMap_ID', 'genepair'], how='left')

    # Fill missing values by DepMap_ID group
    ppi_columns = ['ranked_Essentiality_weighted_PPI', 'Expression_weighted_PPI']

    after = len(df)
    print("Annotated chunk with", before, "rows. Final row count after fill:", after)

    output_path = output_chunks_dir / chunk_path.name
    df.to_parquet(output_path)
    print("Saved annotated chunk to:", output_path, "\n")

# === Process all chunk files ===
chunk_files = sorted(input_chunks_dir.glob("*.parquet"))
print("Found", len(chunk_files), "ranked chunks to annotate with PPI.\n")

for i, chunk_file in enumerate(chunk_files):
    print("Annotating chunk", i + 1, "of", len(chunk_files))
    process_chunk(chunk_file, merged_ppi)

print("All chunks successfully annotated with PPI data.")