#!/usr/bin/env python3
"""Homology triage: split representative proteins into resolved (Branch A
confident) vs unresolved (need Branch B). Reads homology_results.tsv and
produces three outputs: resolved table, unresolved FASTA, trait-candidate FASTA.
Proteins on the curated trait list are always sent to Branch B."""
import time
from pathlib import Path

homology_results = snakemake.input["homology_results"]
rep_proteins = snakemake.input["rep_proteins"]
resolved_out = snakemake.output["resolved"]
unresolved_out = snakemake.output["unresolved"]
trait_candidates_out = snakemake.output["trait_candidates"]
evalue_cutoff = float(snakemake.params["evalue"])
trait_list_path = snakemake.params.get("trait_list", "")
log_file = snakemake.log[0]

for p in [resolved_out, unresolved_out, trait_candidates_out]:
    Path(p).parent.mkdir(parents=True, exist_ok=True)
Path(log_file).parent.mkdir(parents=True, exist_ok=True)
t0 = time.time()

# Read trait list
trait_prots = set()
if trait_list_path and Path(trait_list_path).exists():
    with open(trait_list_path) as f:
        for line in f:
            line = line.strip()
            line and not line.startswith("#") and trait_prots.add(line)

# Parse homology results
results = {}
with open(homology_results) as f:
    header = f.readline()
    for line in f:
        cols = line.strip().split("\t")
        if not cols or cols[0].startswith("#"):
            continue
        results[cols[0]] = {
            "resolved": len(cols) > 1 and cols[1].lower() in ("true", "1"),
            "evidence": cols[2] if len(cols) > 2 else "",
            "ko": cols[3] if len(cols) > 3 else "",
            "ec": cols[4] if len(cols) > 4 else "",
            "go": cols[5] if len(cols) > 5 else "",
            "pfam": cols[6] if len(cols) > 6 else "",
            "cog": cols[7] if len(cols) > 7 else "",
        }

# Read protein sequences
proteins = {}
with open(rep_proteins) as f:
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

# Triage
resolved = []
unresolved_ids = []
trait_seqs = {}
for prot, ann in results.items():
    is_trait = prot in trait_prots
    if ann["resolved"] and not is_trait:
        resolved.append(prot)
    else:
        unresolved_ids.append(prot)
    if is_trait and prot in proteins:
        trait_seqs[prot] = proteins[prot]

with open(resolved_out, "w") as f:
    for p in resolved:
        f.write(p + "\n")
with open(unresolved_out, "w") as f:
    for p in unresolved_ids:
        if p in proteins:
            f.write(f">{p}\n{proteins[p]}\n")
with open(trait_candidates_out, "w") as f:
    for p, seq in trait_seqs.items():
        f.write(f">{p}\n{seq}\n")

elapsed = time.time() - t0
with open(log_file, "w") as log:
    log.write(f"Triage: {len(resolved)} resolved, {len(unresolved_ids)} unresolved, "
              f"{len(trait_seqs)} trait candidates in {elapsed:.1f}s\n")
