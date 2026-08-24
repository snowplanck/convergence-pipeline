#!/usr/bin/env python3
"""Build combined genome x trait matrix with provenance annotations.
Propagates representative-protein annotations back to all genomes via
the dereplication mapping table."""
import argparse, csv, json, sys
from collections import defaultdict

def genome_of(pid):
    return pid.rsplit("_", 1)[0] if "_" in pid else pid.rsplit("|", 1)[0]

def resolve_genome(prot, p2r):
    """Resolve genome_id for a protein, using derep map if available."""
    if prot in p2r:
        return p2r[prot]["genome"]
    return genome_of(prot)

def read_orthogroups(path):
    """Return dict og_id -> [protein_ids] and genome_count matrix."""
    og = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 2:
                continue
            og_id = cols[0]
            members = [c.strip() for c in cols[1:] if c.strip()]
            og[og_id] = members
    return og

def read_homology(path):
    """Return dict protein -> {ko, ec, go, pfam, resolved}."""
    res = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 3:
                continue
            res[cols[0]] = {
                "ko": cols[3] if len(cols) > 3 else "",
                "ec": cols[4] if len(cols) > 4 else "",
                "go": cols[5] if len(cols) > 5 else "",
                "pfam": cols[6] if len(cols) > 6 else "",
                "cog": cols[7] if len(cols) > 7 else "",
                "resolved": cols[1] == "True",
            }
    return res

def read_structure(path):
    """Return dict protein -> {ec, go, fold, stage}."""
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
            }
    return res

def read_derep_map(path):
    """Return dict protein_id -> {rep_id, genome_id} and rep->[proteins]."""
    prot_to_rep = {}
    rep_to_prots = defaultdict(list)
    if not path:
        return prot_to_rep, rep_to_prots
    with open(path) as f:
        header = f.readline()
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 3:
                continue
            prot, rep, genome = cols[0], cols[1], cols[2]
            prot_to_rep[prot] = {"rep": rep, "genome": genome}
            rep_to_prots[rep].append(prot)
    return prot_to_rep, rep_to_prots

def resolve_genome(prot, p2r):
    """Resolve genome_id for a protein, using derep map if available."""
    if prot in p2r:
        return p2r[prot]["genome"]
    return genome_of(prot)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--orthogroups", required=True)
    p.add_argument("--homology", required=True)
    p.add_argument("--structure", required=True)
    p.add_argument("--derep-map", default=None)
    p.add_argument("--genome-list", default=None)
    p.add_argument("--out-matrix", required=True)
    p.add_argument("--out-orthogroups", required=True)
    p.add_argument("--out-provenance", required=True)
    args = p.parse_args()

    og      = read_orthogroups(args.orthogroups)
    homology = read_homology(args.homology)
    struct   = read_structure(args.structure)
    p2r, r2p = read_derep_map(args.derep_map) if args.derep_map else ({}, {})

    # Genome list: derive from annotations if not provided
    genomes = []
    if args.genome_list:
        with open(args.genome_list) as f:
            for line in f:
                g = line.strip()
                if g: genomes.append(g)

    # Build trait set: KO, EC, GO, Pfam, CATH/SCOP fold, OG
    # trait -> {genome -> {count, provenance, confidence}}
    trait_data = defaultdict(lambda: defaultdict(lambda: {"count": 0,
                                                           "provenance": set(),
                                                           "confidence": 0}))

    def add_trait(trait, genome, prov, conf):
        if not trait or trait in ("", "none", "nan"):
            return
        cell = trait_data[trait][genome]
        cell["count"] += 1
        cell["provenance"].add(prov)
        cell["confidence"] = max(cell["confidence"], conf)

    # Process homology annotations (propagate via derep map)
    for prot, ann in homology.items():
        genome = resolve_genome(prot, p2r)
        if ann["resolved"]:
            for ko in ann["ko"].split(";"):
                add_trait(f"KO:{ko}", genome, "homology", 1.0)
            for ec in ann["ec"].split(";"):
                add_trait(f"EC:{ec}", genome, "homology", 1.0)
            for go in ann["go"].split(";"):
                add_trait(f"GO:{go}", genome, "homology", 1.0)
            for pf in ann["pfam"].split(";"):
                add_trait(f"Pfam:{pf}", genome, "homology", 1.0)

    # Process structure annotations
    for prot, ann in struct.items():
        genome = resolve_genome(prot, p2r)
        conf = 0.8 if ann["stage"] == "b1" else (0.7 if ann["stage"] == "b2" else 0.9)
        for ec in ann["ec"].split(";"):
            add_trait(f"EC:{ec}", genome, "structure", conf)
        for go in ann["go"].split(";"):
            add_trait(f"GO:{go}", genome, "structure", conf)
        if ann["fold"]:
            add_trait(f"CATH:{ann['fold']}", genome, "structure", conf)

    # Process orthogroups
    for og_id, members in og.items():
        for mem in members:
            genome = resolve_genome(mem, p2r)
            add_trait(f"OG:{og_id}", genome, "homology", 1.0)

    # Determine consensus provenance
    def provenance_label(prov_set):
        h = "homology" in prov_set
        s = "structure" in prov_set
        if h and s: return "consensus"
        if h:       return "homology_only"
        if s:       return "structure_only"
        return "unknown"

    # Auto-derive genome list from collected trait data if not provided
    if not genomes:
        all_genomes = set()
        for trait, gmap in trait_data.items():
            all_genomes.update(gmap.keys())
        genomes = sorted(all_genomes)

    # Write combined trait matrix
    traits = sorted(trait_data.keys())
    with open(args.out_matrix, "w") as f:
        f.write("trait\t" + "\t".join(genomes) + "\n")
        for trait in traits:
            row = [trait]
            for g in genomes:
                cell = trait_data[trait].get(g, {"count": 0})
                row.append(str(cell["count"]))
            f.write("\t".join(row) + "\n")

    # Write orthogroups presence/absence
    with open(args.out_orthogroups, "w") as f:
        f.write("orthogroup\t" + "\t".join(genomes) + "\n")
        for og_id in sorted(og.keys()):
            row = [og_id]
            for g in genomes:
                present = any(genome_of(m) == g for m in og[og_id])
                row.append("1" if present else "0")
            f.write("\t".join(row) + "\n")

    # Write provenance annotations
    with open(args.out_provenance, "w") as f:
        f.write("trait\tgenome\tcount\tprovenance\tconfidence\n")
        for trait in traits:
            for g in genomes:
                cell = trait_data[trait].get(g)
                if cell and cell["count"] > 0:
                    f.write(f"{trait}\t{g}\t{cell['count']}\t"
                            f"{provenance_label(cell['provenance'])}\t"
                            f"{cell['confidence']:.2f}\n")

    print(f"Trait matrix: {len(traits)} traits x {len(genomes)} genomes")

if __name__ == "__main__":
    main()
