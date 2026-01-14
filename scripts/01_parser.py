import os
import pyarrow.parquet as pq

paralog_pairs_path = "/home/administrator/SLpred/expanded_paralog_pairs.parquet"
output_dir = "/home/administrator/SLpred/chunks"

os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists

# Define the chunk size
chunk_size = 2_000_000 

# Open the Parquet file
parquet_file = pq.ParquetFile(paralog_pairs_path)
total_rows = parquet_file.metadata.num_rows
num_chunks = total_rows // chunk_size + (total_rows % chunk_size != 0)

print(f"Splitting into {num_chunks} chunks...")

for i, batch in enumerate(parquet_file.iter_batches(batch_size=chunk_size)):
    chunk = batch.to_pandas()  # Convert only this chunk to Pandas

    output_file = os.path.join(output_dir, f"all_paralogs_ann_{i}.parquet")
    chunk.to_parquet(output_file)

    print(f"Saved {output_file}")

print("Splitting completed successfully.")