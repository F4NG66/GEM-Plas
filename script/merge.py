import os
import pandas as pd

csv_file = "/root/sequence_clusters_with_reps.csv"
fasta_input_dir = "/root/blast_results"

fasta_output_dir = "/root/fasta_output"
os.makedirs(fasta_output_dir, exist_ok=True)

df = pd.read_csv(csv_file)

for file in os.listdir(fasta_input_dir):
    if file.endswith("_rep_hits.fasta"):
        src_path = os.path.join(fasta_input_dir, file)
        dst_path = os.path.join(fasta_output_dir, file)
        with open(src_path, "r") as src, open(dst_path, "w") as dst:
            dst.write(src.read())

for _, row in df.iterrows():
    cluster_id = row["cluster_id"]
    seq_id = row["sequence_id"]
    sequence = row["sequence"]

    fasta_file = os.path.join(fasta_output_dir, f"cluster{cluster_id}_rep_hits.fasta")

    with open(fasta_file, "a") as f:
        f.write(f">{seq_id}\n")
        f.write(f"{sequence}\n")

