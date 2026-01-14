import os
import pandas as pd
import pyarrow.parquet as pq
import numpy as np
from pathlib import Path
from functools import reduce

# Define paths
go_folder = Path("/home/administrator/SLpred/GO")  # Folder with GO .parquet files
input_chunks_dir = Path("/home/administrator/SLpred/chunks_ppi")
output_chunks_dir = Path("/home/administrator/SLpred/chunks_go")

# Ensure output directory exists
output_chunks_dir.mkdir(parents=True, exist_ok=True)

# Load supporting files
cell_lines = pd.read_csv("/home/administrator/SLpred/DepMap_ID.csv").rename(columns={"0": "DepMap_ID"})

# Dictionary to map filenames to GO value columns
go_dict = {
    'go_BP_expression':'smallest_BP_GO_expression',
    'go_CC_expression':'smallest_CC_GO_expression',
    'go_BP_ranked_essentiality':'smallest_BP_GO_essentiality',
    'go_CC_ranked_essentiality':'smallest_CC_GO_essentiality'
}

# ---------- GO File Processing ----------
def process_go_files(go_files, cell_lines, go_dict):
    results = []

    for file in go_files:
        basename = os.path.basename(file).split(".")[0]
        value_col = go_dict.get(basename)

        if not value_col:
            print(f"Skipping file {file} as it does not match any known GO file names in go_dict.")
            continue

        if file.suffix == '.parquet':
            go_df = pd.read_parquet(file)
            go_df = go_df[go_df['cell_line'].isin(cell_lines['DepMap_ID'])]
            go_df = go_df[['cell_line', 'paralog_pair', value_col]]
            go_df = go_df.rename(columns={
                "cell_line": "DepMap_ID",
                "paralog_pair": "genepair",
                value_col: basename
            })

        results.append(go_df)
        print(f"Processed {file}")

    return results

# Get GO files and process
go_files = sorted(go_folder.glob("*.parquet"))
list_of_go_files = process_go_files(go_files, cell_lines, go_dict)

# Merge on DepMap_ID + genepair
go_merged_df = reduce(lambda left, right: pd.merge(left, right, on=["DepMap_ID", "genepair"], how="outer"), list_of_go_files)
print("GO feature merge complete. Shape:", go_merged_df.shape)
print(go_merged_df[:2])

# Optional renaming for final feature names
rename_dict = {
  'go_BP_expression':'smallest_BP_GO_expression',
  'go_BP_ranked_essentiality':'smallest_BP_GO_essentiality',
  'go_CC_expression':'go_CC_expression',
  'go_CC_ranked_essentiality': 'smallest_CC_GO_essentiality',
}

go_merged_df = go_merged_df.rename(columns=rename_dict)

go_feature_cols = list(rename_dict.values())

# ---------- Annotate Chunks ----------
def annotate_chunk(chunk_path, merged_go_df, output_dir, go_feature_cols):
    print(f"Annotating {chunk_path.name}...")
    df = pd.read_parquet(chunk_path)

    df = pd.merge(df, merged_go_df, on=["DepMap_ID", "genepair"], how="left")
    output_path = output_dir / chunk_path.name
    df.to_parquet(output_path)
    print(df[:2])
    print(f"Saved annotated chunk: {output_path.name}")

# Process all chunks
chunk_files = sorted(input_chunks_dir.glob("*.parquet"))
print(f"Annotating {len(chunk_files)} chunks...\n")

for chunk_file in chunk_files:
    annotate_chunk(chunk_file, go_merged_df, output_chunks_dir, go_feature_cols)

print("GO annotation complete.")