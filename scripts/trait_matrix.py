#!/usr/bin/env python3
"""Module 4 — build combined genome x trait matrix with provenance.
Reads orthogroups, homology results, structure results, derep map.
Writes trait_matrix_combined.tsv, trait_provenance.tsv, orthogroups.tsv."""
import csv
import json
import time
from collections import defaultdict
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
import genome_utils

orthogroups_in = snakemake.input["orthogroups"]
gene_count_in = snakemake.input["gene_count"]
homology_in = snakemake.input["homology_resolved"]
structure_in = snakemake.input["structure_results"]
mapping_in = snakemake.input["mapping"]
protein_to_genome_in = snakemake.input["protein_to_genome"]
matrix_out = snakemake.output["matrix"]
provenance_out = snakemake.output["provenance"]
orthogroups_out = snakemake.output["orthogroups_out"]
genome_dir = snakemake.params["genome_dir"]
log_file = snakemake.log[0]

for p in [matrix_out, provenance_out, orthogroups_out]:
    Path(p).parent.mkdir(parents=True, exist_ok=True)
Path(log_file).parent.mkdir(parents=True, exist_ok=True)
with open(log_file, "w") as _lf:
    pass  # truncate/create
t0 = time.time()

def read_fasta_ids(directory):
    ids = set()
    dpath = Path(directory)
    if dpath.exists():
        for p in dpath.iterdir():
            if p.suffix in (".fna", ".fa", ".fasta", ".gbk", ".gbff"):
                ids.add(p.stem)
    return ids

p2g_map = genome_utils.load_protein_to_genome(protein_to_genome_in)
def resolve_genome(prot):
    return genome_utils.resolve_genome(prot, p2g_map, log_file)

# Read homology resolved list
resolved_set = set()
if Path(homology_in).exists():
    with open(homology_in) as f:
        for line in f:
            resolved_set.add(line.strip())

# Read orthogroups
og = {}
with open(orthogroups_in) as f:
    header = f.readline()
    for line in f:
        cols = line.strip().split("\t")
        if len(cols) >= 2:
            og[cols[0]] = [g.strip() for c in cols[1:] if c.strip() for g in c.split(",") if g.strip()]

# Read structure results
struct = {}
with open(structure_in) as f:
    header = f.readline()
    for line in f:
        cols = line.strip().split("\t")
        if len(cols) >= 7:
            struct[cols[0]] = {"ec": cols[4] if len(cols) > 4 else "",
                                "fold": cols[6] if len(cols) > 6 else "",
                                "stage": cols[2]}

# Genome list: derive from annotations
all_genomes = read_fasta_ids(genome_dir)
for prot in list(struct.keys()) + list(resolved_set):
    all_genomes.add(resolve_genome(prot))
for members in og.values():
    for m in members:
        all_genomes.add(resolve_genome(m))
genomes = sorted(all_genomes)

# Build trait data
trait_data = defaultdict(lambda: defaultdict(lambda: {"count": 0, "provenance": set(), "confidence": 0}))

def add_trait(trait, genome, prov, conf):
    if not trait or trait in ("", "none"):
        return
    cell = trait_data[trait][genome]
    cell["count"] += 1
    cell["provenance"].add(prov)
    cell["confidence"] = max(cell["confidence"], conf)

# Homology
for prot in resolved_set:
    genome = resolve_genome(prot)
    add_trait(f"OG:resolved", genome, "homology", 1.0)

# Structure
for prot, ann in struct.items():
    genome = resolve_genome(prot)
    conf = 0.8 if ann["stage"] == "b1" else (0.7 if ann["stage"] == "b2" else 0.9)
    for ec in ann["ec"].split(";"):
        if ec:
            add_trait(f"EC:{ec}", genome, "structure", conf)
    if ann["fold"]:
        add_trait(f"CATH:{ann['fold']}", genome, "structure", conf)

# Orthogroups
for og_id, members in og.items():
    for mem in members:
        add_trait(f"OG:{og_id}", resolve_genome(mem), "homology", 1.0)

def provenance_label(s):
    h, st = "homology" in s, "structure" in s
    if h and st:
        return "consensus"
    if h:
        return "homology_only"
    if st:
        return "structure_only"
    return "unknown"

# Write matrix
traits = sorted(trait_data.keys())
with open(matrix_out, "w") as f:
    f.write("trait\t" + "\t".join(genomes) + "\n")
    for trait in traits:
        row = [trait]
        for g in genomes:
            cell = trait_data[trait].get(g, {"count": 0})
            row.append(str(cell["count"]))
        f.write("\t".join(row) + "\n")

# Write provenance
with open(provenance_out, "w") as f:
    f.write("trait\tgenome\tcount\tprovenance\tconfidence\n")
    for trait in traits:
        for g in genomes:
            cell = trait_data[trait].get(g)
            if cell and cell["count"] > 0:
                f.write(f"{trait}\t{g}\t{cell['count']}\t"
                        f"{provenance_label(cell['provenance'])}\t{cell['confidence']:.2f}\n")

# Write orthogroups presence/absence
with open(orthogroups_out, "w") as f:
    f.write("orthogroup\t" + "\t".join(genomes) + "\n")
    for og_id in sorted(og.keys()):
        row = [og_id]
        for g in genomes:
            present = any(resolve_genome(m) == g for m in og[og_id])
            row.append("1" if present else "0")
        f.write("\t".join(row) + "\n")

elapsed = time.time() - t0
with open(log_file, "a") as log:
    log.write(f"Trait matrix: {len(traits)} traits x {len(genomes)} genomes in {elapsed:.1f}s\n")
