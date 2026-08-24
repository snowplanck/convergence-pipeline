#!/usr/bin/env python3
"""Reconstruct complete KEGG modules from KofamScan hits using the
all-but-one reaction tolerance rule."""
import argparse, csv, json, sys
from collections import defaultdict

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--kofam", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--tolerance", type=int, default=1)
    p.add_argument("--modules-json", default=None)
    args = p.parse_args()

    # Parse KofamScan hits: protein_id -> set of KO
    ko_hits = defaultdict(set)
    with open(args.kofam) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if parts[0].startswith("#"):
                continue
            # KofamScan detail-tsv columns vary; use first meaningful hit per line
            # Format: protein, KO, score, evalue (threshold lines start with *)
            if parts[0].strip() == "*":
                continue
            protein = parts[0].strip()
            # Find the KO column (K<number>)
            for col in parts[1:]:
                col = col.strip()
                if col.startswith("K") and len(col) == 6:
                    ko_hits[protein].add(col)
                    break

    # If module definitions are unavailable, emit per-protein KO assignments
    if not args.modules_json:
        with open(args.output, "w") as f:
            f.write("genome_id\tmodule_id\tcompleteness\tkpresent\tkmissing\tkos\n")
            for prot, kos in ko_hits.items():
                genome = prot.rsplit("_", 1)[0] if "_" in prot else prot
                for ko in kos:
                    f.write(f"{genome}\t{ko}\t1.0\t1\t0\t{ko}\n")
        return

    with open(args.modules_json) as f:
        modules = json.load(f)

    with open(args.output, "w") as f:
        f.write("genome_id\tmodule_id\tcompleteness\tkpresent\tkmissing\tkos\n")
        for mod_id, req_kos in modules.items():
            # Aggregate KOs across all proteins in each genome
            genome_kos = defaultdict(set)
            for prot, kos in ko_hits.items():
                genome = prot.rsplit("_", 1)[0] if "_" in prot else prot
                genome_kos[genome] |= kos
            for genome, present in genome_kos.items():
                found = present & set(req_kos)
                missing = set(req_kos) - found
                n_found = len(found)
                n_total = len(req_kos)
                if n_total == 0:
                    continue
                completeness = n_found / n_total
                # Tolerance: count as complete if missing <= tolerance
                if len(missing) <= args.tolerance and n_found >= 2:
                    f.write(f"{genome}\t{mod_id}\t{completeness:.3f}\t"
                            f"{n_found}\t{len(missing)}\t{';'.join(sorted(found))}\n")

if __name__ == "__main__":
    main()
