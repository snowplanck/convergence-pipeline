#!/usr/bin/env python3
"""Detect candidate Non-homologous Isofunctional Enzyme (NISE) events:
proteins sharing the same function (EC/GO) but from different orthogroups
with no detectable sequence homology and/or different structural folds."""
import argparse, csv, json, sys
from collections import defaultdict

def read_trait_matrix(path):
    """Return dict trait -> {genome: count} and genome list."""
    data = {}
    genomes = []
    with open(path) as f:
        header = f.readline()
        genomes = header.strip().split("\t")[1:]
        for line in f:
            cols = line.strip().split("\t")
            trait = cols[0]
            data[trait] = {genomes[i]: int(cols[i+1]) for i in range(len(genomes))}
    return data, genomes

def read_homology(path):
    """Return dict protein -> {ko, ec, go, resolved}."""
    res = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 5:
                continue
            res[cols[0]] = {
                "ec": cols[4] if len(cols) > 4 else "",
                "go": cols[5] if len(cols) > 5 else "",
                "resolved": cols[1] == "True",
            }
    return res

def read_structure(path):
    """Return dict protein -> {ec, go, fold, stage, tool}."""
    res = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 7:
                continue
            res[cols[0]] = {
                "ec": cols[4] if len(cols) > 4 else "",
                "go": cols[5] if len(cols) > 5 else "",
                "fold": cols[6] if len(cols) > 6 else "",
                "stage": cols[2],
                "tool": cols[3],
            }
    return res

def read_orthogroups(path):
    """Return dict protein -> og_id."""
    prot_to_og = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 2:
                continue
            og_id = cols[0]
            for mem in cols[1:]:
                mem = mem.strip()
                if mem:
                    prot_to_og[mem] = og_id
    return prot_to_og

def genome_of(pid):
    return pid.rsplit("_", 1)[0] if "_" in pid else pid.rsplit("|", 1)[0]

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--trait-matrix", required=True)
    p.add_argument("--homology", required=True)
    p.add_argument("--structure", required=True)
    p.add_argument("--orthogroups", required=True)
    p.add_argument("--mmseqs-evalue", type=float, default=1e-3)
    p.add_argument("--hhblits-prob", type=float, default=95.0)
    p.add_argument("--min-tools", type=int, default=2)
    p.add_argument("--out", required=True)
    p.add_argument("--summary", required=True)
    args = p.parse_args()

    trait_data, genomes = read_trait_matrix(args.trait_matrix)
    homology = read_homology(args.homology)
    struct   = read_structure(args.structure)
    prot2og  = read_orthogroups(args.orthogroups)

    # Collect all proteins per EC/GO function
    func_to_prots = defaultdict(set)
    for prot, ann in homology.items():
        for ec in ann["ec"].split(";"):
            if ec: func_to_prots[f"EC:{ec}"].add(prot)
        for go in ann["go"].split(";"):
            if go: func_to_prots[f"GO:{go}"].add(prot)
    for prot, ann in struct.items():
        for ec in ann["ec"].split(";"):
            if ec: func_to_prots[f"EC:{ec}"].add(prot)
        for go in ann["go"].split(";"):
            if go: func_to_prots[f"GO:{go}"].add(prot)

    candidates = []

    for func, prots in func_to_prots.items():
        if len(prots) < 2:
            continue

        # Group proteins by orthogroup
        og_groups = defaultdict(set)
        for prot in prots:
            og = prot2og.get(prot, "unassigned")
            og_groups[og].add(prot)

        # Need >= 2 distinct orthogroups for convergence
        if len(og_groups) < 2:
            continue

        # Check fold diversity
        folds = defaultdict(set)
        tools = defaultdict(set)
        for prot in prots:
            s = struct.get(prot, {})
            if s.get("fold"):
                folds[s["fold"]].add(prot)
            if s.get("tool"):
                tools[s["tool"]].add(prot)

        n_folds = len(folds)
        n_tools = len(tools)

        # NISE criteria: multiple OGs + (different folds OR low homology)
        is_nise = (n_folds >= 2) or (len(og_groups) >= 2 and n_tools >= args.min_tools)

        if is_nise:
            # Compute confidence score
            confidence = min(1.0, (n_folds * 0.3 + n_tools * 0.2 + len(og_groups) * 0.1))

            # Collect genomes involved
            involved_genomes = sorted(set(genome_of(p) for p in prots))

            candidates.append({
                "function": func,
                "n_proteins": len(prots),
                "n_orthogroups": len(og_groups),
                "n_folds": n_folds,
                "n_tools": n_tools,
                "orthogroups": ";".join(sorted(og_groups.keys())),
                "folds": ";".join(sorted(folds.keys())),
                "genomes": ";".join(involved_genomes),
                "confidence": round(confidence, 3),
                "evidence": f"OGs={len(og_groups)},folds={n_folds},tools={n_tools}",
            })

    # Sort by confidence descending
    candidates.sort(key=lambda x: x["confidence"], reverse=True)

    with open(args.out, "w") as f:
        fields = ["function", "n_proteins", "n_orthogroups", "n_folds", "n_tools",
                   "orthogroups", "folds", "genomes", "confidence", "evidence"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for c in candidates:
            w.writerow(c)

    summary = {
        "total_functions_checked": len(func_to_prots),
        "nise_candidates": len(candidates),
        "high_confidence": sum(1 for c in candidates if c["confidence"] >= 0.7),
    }
    with open(args.summary, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"NISE detection: {len(candidates)} candidates from {len(func_to_prots)} functions")

if __name__ == "__main__":
    main()
