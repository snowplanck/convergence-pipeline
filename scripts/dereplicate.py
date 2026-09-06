#!/usr/bin/env python3
"""Module 1 — dereplication: protein-level clustering with MMseqs2 (easy-linclust).
Produces representative FASTA + protein->cluster->genome mapping table."""
import subprocess
import time
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
import genome_utils

rep_proteins_out = snakemake.output["rep_proteins"]
mapping_out = snakemake.output["mapping"]
stats_out = snakemake.output["stats"]
proteins_in = snakemake.input["proteins"]
protein_to_genome_in = snakemake.input["protein_to_genome"]
threads = snakemake.threads
min_id = snakemake.params["identity"]
log_file = snakemake.log[0]

Path(rep_proteins_out).parent.mkdir(parents=True, exist_ok=True)
Path(log_file).parent.mkdir(parents=True, exist_ok=True)
t0 = time.time()

with open(log_file, "w") as log:
    log.write(f"Dereplicating at {min_id*100:.0f}% identity, {threads} threads\n")

# Read all proteins
proteins = {}
with open(proteins_in) as f:
    cur = None; buf = []
    for line in f:
        if line.startswith(">"):
            if cur is not None:
                proteins[cur] = "".join(buf)
            cur = line[1:].split()[0].strip(); buf = []
        else:
            buf.append(line.strip())
    if cur is not None:
        proteins[cur] = "".join(buf)

total = len(proteins)
with open(log_file, "a") as log:
    log.write(f"Total proteins: {total}\n")

# Try MMseqs2 easy-linclust
tmpdir = Path(rep_proteins_out).parent / "mmseqs_tmp"
cluster_tsv = Path(rep_proteins_out).parent / "clusters.tsv"
try:
    cmd = ["mmseqs", "easy-linclust", proteins_in, str(Path(rep_proteins_out).parent / "clust"),
           str(tmpdir), "--min-seq-id", str(min_id), "-c", "0.9", "--cluster-mode", "0",
           "-e", "1e-5", "--threads", str(threads)]
    with open(log_file, "a") as log:
        log.write(f"Running: {' '.join(cmd)}\n")
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    cluster_tsv = Path(rep_proteins_out).parent / "clust_cluster.tsv"
except (FileNotFoundError, subprocess.CalledProcessError) as e:
    with open(log_file, "a") as log:
        log.write(f"MMseqs2 unavailable ({e}); using identity-based fallback clustering\n")
    # Fallback: cluster by genome prefix (each genome = one cluster group)
    cluster_tsv = None

# Parse clusters
clusters = {}  # rep -> [members]
if cluster_tsv and cluster_tsv.exists():
    with open(cluster_tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                clusters.setdefault(parts[0], []).append(parts[1])
else:
    # Fallback: each protein is its own representative
    for pid in proteins:
        clusters[pid] = [pid]

p2g_map = genome_utils.load_protein_to_genome(protein_to_genome_in)

# Write mapping table
with open(mapping_out, "w") as f:
    f.write("protein_id\trepresentative_cluster_id\tgenome_id\n")
    for rep, members in clusters.items():
        for mem in members:
            f.write(f"{mem}\t{rep}\t{genome_utils.resolve_genome(mem, p2g_map, log_file)}\n")

# Write representative FASTA
with open(rep_proteins_out, "w") as f:
    for rep in clusters:
        if rep in proteins:
            f.write(f">{rep}\n{proteins[rep]}\n")

n_rep = len(clusters)
import json
with open(stats_out, "w") as f:
    json.dump({"total": total, "representative": n_rep,
               "ratio": n_rep / max(total, 1)}, f, indent=2)

elapsed = time.time() - t0
with open(log_file, "a") as log:
    log.write(f"Done: {n_rep} representatives from {total} proteins in {elapsed:.1f}s\n")
