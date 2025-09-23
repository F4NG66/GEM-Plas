#!/bin/bash
# 用法：bash batch_blast.sh

input="/root/representatives.fasta"
db="/root/databases/uniref/uniref90_blast"
outdir="/root/blast_results"
tmp="/root/blast_tmp"

mkdir -p "$outdir" "$tmp"


grep "^>" "$input" | sed 's/^>//' > "$tmp/ids.txt"

while read -r id; do
    echo "=== Processing $id ==="

 
    awk -v id="$id" '/^>/{p=($0==">"id)} p' "$input" > "$tmp/$id.fasta"

    
    blastp \
      -query "$tmp/$id.fasta" \
      -db "$db" \
      -outfmt "6 sseqid pident length evalue bitscore" \
      -max_target_seqs 500 \
      -evalue 1e-4 \
      -num_threads 4 \
      -out "$tmp/${id}_all.tsv"

    awk '$2 >= 40' "$tmp/${id}_all.tsv" | sort -k5,5nr | head -n 200 | cut -f1 > "$tmp/${id}_top100_ids.txt"


    blastdbcmd -db "$db" -entry_batch "$tmp/${id}_top100_ids.txt" -out "$outdir/${id}_hits.fasta"


    rm -f "$tmp/$id.fasta" "$tmp/${id}_all.tsv" "$tmp/${id}_top100_ids.txt"
done < "$tmp/ids.txt"


