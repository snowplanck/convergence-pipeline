#!/usr/bin/env python3
"""Parse Foldseek structural search results, extract fold/EC annotations,
and split proteins into resolved vs residual sets."""
import argparse, csv, sys
from collections import defaultdict

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--foldseek", required=True)
    p.add_argument("--input-fasta", required=True)
    p.add_argument("--tm-score-cutoff", type=float, default=0.5)
    p.add_argument("--fident-cutoff", type=float, default=0.3)
    p.add_argument("--out-predictions", required=True)
    p.add_argument("--out-resolved", required=True)
    p.add_argument("--out-unresolved", required=True)
    args = p.parse_args()

    # Foldseek output columns: query, target, evalue, bits, fident, alnlen, qcov, tcov, tmid
    best_hits = {}
    with open(args.foldseek) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            query = cols[0]
            target = cols[1]
            evalue = float(cols[2])
            fident = float(cols[4])
            tmid = float(cols[8]) if cols[8] != "N/A" else 0.0
            # TM-score proxy: use tmid (structural similarity)
            if query not in best_hits or evalue < best_hits[query]["evalue"]:
                best_hits[query] = {
                    "target": target, "evalue": evalue, "fident": fident,
                    "tmid": tmid
                }

    # Read sequences
    seqs = {}
    with open(args.input_fasta) as f:
        cur = None; buf = []
        for line in f:
            if line.startswith(">"):
                if cur: seqs[cur] = "".join(buf)
                cur = line[1:].split()[0].strip(); buf = []
            else: buf.append(line.strip())
        if cur: seqs[cur] = "".join(buf)

    resolved = {}
    unresolved = []
    for prot in seqs:
        hit = best_hits.get(prot)
        if hit and hit["tmid"] >= args.tm_score_cutoff:
            # Extract CATH/SCOP fold ID from target (e.g., 1.10.10.10_1abc)
            fold_id = hit["target"].split("_")[0] if "_" in hit["target"] else hit["target"]
            resolved[prot] = {
                "fold": fold_id, "evalue": hit["evalue"],
                "tmid": hit["tmid"], "fident": hit["fident"],
                "target": hit["target"]
            }
        else:
            unresolved.append(prot)

    with open(args.out_predictions, "w") as f:
        f.write("protein\tfold\ttarget\tevalue\ttmid\tfident\tresolved\n")
        for prot, info in resolved.items():
            f.write(f"{prot}\t{info['fold']}\t{info['target']}\t{info['evalue']}\t"
                    f"{info['tmid']:.3f}\t{info['fident']:.3f}\tTrue\n")
        for prot in unresolved:
            f.write(f"{prot}\tnone\t\t\t\t\tFalse\n")

    with open(args.out_resolved, "w") as f:
        for prot in resolved:
            f.write(f">{prot}\n{seqs[prot]}\n")

    with open(args.out_unresolved, "w") as f:
        for prot in unresolved:
            f.write(f">{prot}\n{seqs[prot]}\n")

    print(f"B2 Foldseek: {len(resolved)} resolved, {len(unresolved)} unresolved")

if __name__ == "__main__":
    main()
