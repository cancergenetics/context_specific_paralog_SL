import pandas as pd
from pathlib import Path

# Path to directory with annotated Parquet chunks
chunk_dir = Path("/home/administrator/SLpred/chunks_processed")

# Output log file
log_path = Path("feature_report.txt")
log_file = open(log_path, "w")

# Feature definitions
feature_columns_1 = [
    'rMaxExp_A1A2', 'rMinExp_A1A2',
    'max_ranked_A1A2', 'min_ranked_A1A2',
    'max_cn', 'min_cn', 'Protein_Altering', 'Damaging',
    'min_sequence_identity', 'prediction_score',
    'weighted_PPI_expression', 'smallest_BP_GO_expression', 'go_CC_expression',
    'weighted_PPI_essentiality', 'smallest_BP_GO_essentiality', 'smallest_CC_GO_essentiality'
]

feature_columns_2 = feature_columns_1 + [
    'closest', 'WGD', 'family_size', 'cds_length_ratio', 'shared_domains',
    'has_pombe_ortholog', 'has_essential_pombe_ortholog',
    'has_cerevisiae_ortholog', 'has_essential_cerevisiae_ortholog',
    'conservation_score', 'mean_age', 'either_in_complex',
    'mean_complex_essentiality', 'colocalisation', 'interact', 'n_total_ppi',
    'fet_ppi_overlap', 'gtex_spearman_corr', 'gtex_min_mean_expr', 'gtex_max_mean_expr'
]

# Get chunk files
chunk_files = sorted(chunk_dir.glob("*.parquet"))
print(f"Scanning {len(chunk_files)} chunks...\n", file=log_file)

# Process each chunk
for chunk_path in chunk_files:
    print(f"Checking file: {chunk_path.name}", file=log_file)
    df = pd.read_parquet(chunk_path)

    # Check missing columns
    missing_1 = [col for col in feature_columns_1 if col not in df.columns]
    missing_2 = [col for col in feature_columns_2 if col not in df.columns]

    if not missing_1:
        print("All feature_columns_1 present.", file=log_file)
    else:
        print("Missing columns in feature_columns_1:", missing_1, file=log_file)

    if not missing_2:
        print("All feature_columns_2 present.", file=log_file)
    else:
        print("Missing columns in feature_columns_2:", missing_2, file=log_file)

    # NA counts
    print("  NA counts for present columns (feature_columns_2):", file=log_file)
    for col in feature_columns_2:
        if col in df.columns:
            na_count = df[col].isna().sum()
            print(f"    {col}: {na_count}", file=log_file)
    print("", file=log_file)

log_file.close()
print(f"Done. Report saved to: {log_path}")