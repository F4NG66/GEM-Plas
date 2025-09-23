import os
import torch
import numpy as np
from Bio import AlignIO
from torch_geometric.data import Data
import random
import time
from tqdm import tqdm

fasta_dir = "/root/msa_output1"
graph_dir = "/root/graphs1"
os.makedirs(graph_dir, exist_ok=True)

def shannon_entropy(column):
    freqs = [column.count(aa)/len(column) for aa in set(column) if aa != "-"]
    return -sum(p*np.log2(p) for p in freqs if p > 0)

def mutual_information(col_i, col_j):
    aa_set = list(set(col_i + col_j) - {"-"})
    if not aa_set:
        return 0.0
    mi = 0.0
    for a in aa_set:
        for b in aa_set:
            p_ab = sum((x==a and y==b) for x,y in zip(col_i,col_j)) / len(col_i)
            p_a = col_i.count(a)/len(col_i)
            p_b = col_j.count(b)/len(col_j)
            if p_ab > 0:
                mi += p_ab * np.log2(p_ab/(p_a*p_b))
    return mi

files = [f for f in os.listdir(fasta_dir) if f.endswith(".fasta")]
total = len(files)
start_all = time.time()

for idx, file in enumerate(tqdm(files, desc="Processing clusters")):
    cluster_start = time.time()
    cluster_name = file.replace(".fasta", "")
    msa_fasta = os.path.join(fasta_dir, file)

    alignment = AlignIO.read(msa_fasta, "fasta")
    L = alignment.get_alignment_length()
    N = len(alignment)

    valid_cols = []
    for i in range(L):
        col = [rec.seq[i] for rec in alignment]
        gap_ratio = col.count("-") / N
        if gap_ratio < 0.5:
            valid_cols.append(i)

    node_features = []
    for i in valid_cols:
        col = [rec.seq[i] for rec in alignment]
        entropy = shannon_entropy(col)
        aa_counts = [col.count(aa)/N for aa in "ACDEFGHIKLMNPQRSTVWY"]
        node_features.append([entropy] + aa_counts)
    node_features = torch.tensor(node_features, dtype=torch.float)

    edge_index, edge_attr = [], []

    # backbone 边
    for k in range(len(valid_cols) - 1):
        edge_index.append([k, k+1])
        edge_index.append([k+1, k])
        edge_attr.append([1.0])
        edge_attr.append([1.0])

    for idx_i, i in enumerate(valid_cols):
        col_i = [rec.seq[i] for rec in alignment]
        scores = []
        for idx_j, j in enumerate(valid_cols):
            if idx_j <= idx_i:
                continue
            col_j = [rec.seq[j] for rec in alignment]
            mi = mutual_information(col_i, col_j)
            if mi > 0.3:
                scores.append((idx_j, mi))
        scores = sorted(scores, key=lambda x: x[1], reverse=True)[:30]
        for idx_j, mi in scores:
            edge_index.append([idx_i, idx_j])
            edge_index.append([idx_j, idx_i])
            edge_attr.append([mi])
            edge_attr.append([mi])

    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    edge_attr = torch.tensor(edge_attr, dtype=torch.float)

    graph = Data(x=node_features, edge_index=edge_index, edge_attr=edge_attr)
    torch.save(graph, os.path.join(graph_dir, f"{cluster_name}.pt"))

    cluster_time = time.time() - cluster_start
    elapsed = time.time() - start_all
    avg_time = elapsed / (idx + 1)
    remaining = avg_time * (total - idx - 1)

    print(f"\n {cluster_name}  {cluster_time:.1f}s | "
          f"{idx+1}/{total} | {remaining/60:.1f} min")


