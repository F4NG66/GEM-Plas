import os
import torch
import numpy as np
from Bio import AlignIO


fasta_dir = "/root/fasta_output"
msa_dir = "/root/msa_output1"
graph_dir = "/root/graphs"
os.makedirs(msa_dir, exist_ok=True)
os.makedirs(graph_dir, exist_ok=True)


def run_mafft(input_fasta, output_fasta):
    os.system(f"mafft --auto {input_fasta} > {output_fasta}")

for file in os.listdir(fasta_dir):
    if not file.endswith(".fasta"):
        continue

    cluster_name = file.replace(".fasta", "")
    input_fasta = os.path.join(fasta_dir, file)
    msa_fasta = os.path.join(msa_dir, f"{cluster_name}_aligned.fasta")

    print(f" Processing {cluster_name}...")


    run_mafft(input_fasta, msa_fasta)

