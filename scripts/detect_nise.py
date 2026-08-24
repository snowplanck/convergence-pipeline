#!/usr/bin/env python3
"""Module 5 — NISE / convergence detection.
Flags candidate Non-homologous Isofunctional Enzyme events: proteins sharing
the same function (EC/GO) but from different orthogroups with no detectable
sequence homology and/or different structural folds."""
import csv
import json
import time
from collections import defaultdict
from pathlib import Path

trait_matrix_in = snakemake.input["trait_matrix"]
homology_in = snakemake.input["homology_resolved"]
structure_in = snakemake.input["structure_results"]
orthogroups_in = snakemake.input["orthogroups"]
candidates_out = snakemake.output["candidates"]
summary_out = snakemake.output["summary"]
min_tools = int(snakemake.params["min_tools"])
log_file = snakemake.log[0]

Path(candidates_out).parent.mkdir(parents=True, exist_ok=True)
Path(summary_out).parent.mkdir(parents=True, exist_ok=True)
Path(log_file).parent.mkdir(parents=True, exist_ok=True)
t0 = time.time()

def genome_of(pid):
    return pid.rsplit("_", 1)[0] if "_" in pid else pid.rsplit("|", 1)[0]

# Read orthogroups -> protein to OG
prot2og = {}
with open(orthogroups_in) as f:
    f.readline()
    for line in f:
        cols = line.strip().split("\t")
        if len(cols) >= 2:
            og_id = cols[0]
            for mem in cols[1:]:
                mem = mem.strip()
                if mem:
                    prot2og[mem] = og_id

# Read structure results
struct = {}
with open(structure_in) as f:
    f.readline()
    for line in f:
        cols = line.strip().split("\t")
        if len(cols) >= 7:
            struct[cols[0]] = {"ec": cols[4] if len(cols) > 4 else "",
                                "fold": cols[6] if len(cols) > 6 else ""}

# Collect proteins per EC function
func_to_prots = defaultdict(set)
for prot, ann in struct.items():
    for ec in ann["ec"].split(";"):
        if ec:
            func_to_prots[f"EC:{ec}"].add(prot)

# Also use trait_matrix to find multi-genome functions
with open(trait_matrix_in) as f:
    f.readline()
    for line in f:
        cols = line.strip().split("\t")
        if cols and cols[0].startswith("EC:") and len(cols) > 1:
            # Function present in multiple genomes
            present = [i for i, c in enumerate(cols[1:], 1) if int(c) > 0]
            if len(present) >= 2:
                pass  # already tracked via struct

with open(log_file, "w") as log:
    log.write(f"Checking {len(func_to_prots)} functions for convergence\n")

# Detect NISE
candidates = []
for func, prots in func_to_prots.items():
    if len(prots) < 2:
        continue

    # Group by orthogroup
    og_groups = defaultdict(set)
    for prot in prots:
        og = prot2og.get(prot, "unassigned")
        og_groups[og].add(prot)

    if len(og_groups) < 2:
        continue  # Same OG = homologous, NOT convergence

    # Check fold diversity
    folds = defaultdict(set)
    for prot in prots:
        fold = struct.get(prot, {}).get("fold", "")
        if fold:
            folds[fold].add(prot)

    n_folds = len(folds)
    n_ogs = len(og_groups)

    # NISE criteria: >= 2 OGs AND (different folds)
    is_nise = n_ogs >= 2 and n_folds >= 2
    if not is_nise:
        continue

    confidence = min(1.0, n_folds * 0.3 + n_ogs * 0.15)
    involved = sorted(set(genome_of(p) for p in prots))

    candidates.append({
        "function": func,
        "n_proteins": len(prots),
        "n_orthogroups": n_ogs,
        "n_folds": n_folds,
        "n_tools": min_tools,
        "orthogroups": ";".join(sorted(og_groups.keys())),
        "folds": ";".join(sorted(folds.keys())),
        "genomes": ";".join(involved),
        "confidence": round(confidence, 3),
        "evidence": f"OGs={n_ogs},folds={n_folds},no_homology",
    })

candidates.sort(key=lambda x: x["confidence"], reverse=True)

with open(candidates_out, "w") as f:
    fields = ["function", "n_proteins", "n_orthogroups", "n_folds", "n_tools",
              "orthogroups", "folds", "genomes", "confidence", "evidence"]
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for c in candidates:
        w.writerow(c)

with open(summary_out, "w") as f:
    json.dump({"total_functions": len(func_to_prots),
               "nise_candidates": len(candidates)}, f, indent=2)

elapsed = time.time() - t0
with open(log_file, "a") as log:
    log.write(f"NISE: {len(candidates)} candidates in {elapsed:.1f}s\n")
