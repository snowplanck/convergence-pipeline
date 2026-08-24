#!/usr/bin/env python3
"""Merge ProteInfer + CLEAN predictions, apply confidence threshold,
split into resolved vs unresolved sets."""
import argparse, csv, sys
from collections import defaultdict

def parse_proteinfer(path):
    """ProteInfer output: protein, label, probability."""
    res = defaultdict(list)
    if not path: return res
    with open(path) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 3 or cols[0].startswith("#"):
                continue
            prot, label, prob = cols[0], cols[1], float(cols[2])
            res[prot].append({"tool": "proteinfer", "label": label, "prob": prob})
    return res

def parse_clean(path):
    """CLEAN output: protein, ec, probability."""
    res = defaultdict(list)
    if not path: return res
    with open(path) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 3 or cols[0].startswith("#"):
                continue
            prot, ec, prob = cols[0], cols[1], float(cols[2])
            res[prot].append({"tool": "clean", "ec": ec, "prob": prob})
    return res

def read_fasta(path):
    seqs = {}
    with open(path) as f:
        cur = None; buf = []
        for line in f:
            if line.startswith(">"):
                if cur: seqs[cur] = "".join(buf)
                cur = line[1:].split()[0].strip(); buf = []
            else: buf.append(line.strip())
        if cur: seqs[cur] = "".join(buf)
    return seqs

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--proteinfer", default=None)
    p.add_argument("--clean", default=None)
    p.add_argument("--input-fasta", required=True)
    p.add_argument("--confidence", type=float, default=0.9)
    p.add_argument("--out-predictions", required=True)
    p.add_argument("--out-resolved", required=True)
    p.add_argument("--out-unresolved", required=True)
    args = p.parse_args()

    pf = parse_proteinfer(args.proteinfer)
    cl = parse_clean(args.clean)
    seqs = read_fasta(args.input_fasta)

    all_prots = list(seqs.keys())
    resolved = {}
    unresolved = []

    for prot in all_prots:
        best = None
        for hit in pf.get(prot, []):
            if best is None or hit["prob"] > best["prob"]:
                best = {"tool": "proteinfer", "label": hit["label"],
                        "prob": hit["prob"]}
        for hit in cl.get(prot, []):
            if best is None or hit["prob"] > best.get("prob", 0):
                best = {"tool": "clean", "ec": hit["ec"], "prob": hit["prob"]}
        if best and best["prob"] >= args.confidence:
            resolved[prot] = best
        else:
            unresolved.append(prot)

    with open(args.out_predictions, "w") as f:
        f.write("protein\ttool\tlabel\tec\tprobability\tresolved\n")
        for prot, info in resolved.items():
            f.write(f"{prot}\t{info['tool']}\t{info.get('label','')}\t"
                    f"{info.get('ec','')}\t{info['prob']:.4f}\tTrue\n")
        for prot in unresolved:
            f.write(f"{prot}\tnone\t\t\t0.0\tFalse\n")

    with open(args.out_resolved, "w") as f:
        for prot in resolved:
            f.write(f">{prot}\n{seqs[prot]}\n")

    with open(args.out_unresolved, "w") as f:
        for prot in unresolved:
            f.write(f">{prot}\n{seqs[prot]}\n")

    print(f"B1 ML: {len(resolved)} resolved, {len(unresolved)} unresolved")

if __name__ == "__main__":
    main()
