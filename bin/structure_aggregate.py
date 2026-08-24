#!/usr/bin/env python3
"""Aggregate B1/B2/B3 structure predictions into a single table with
provenance tags and stage resolution metadata."""
import argparse, csv, json, sys

def read_predictions(path, stage):
    recs = {}
    if not path: return recs
    with open(path) as f:
        header = f.readline()
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 6:
                continue
            recs[cols[0]] = {
                "stage": stage,
                "tool": cols[1] if len(cols) > 1 else "",
                "label": cols[2] if len(cols) > 2 else "",
                "ec": cols[3] if len(cols) > 3 else "",
                "go": cols[7] if len(cols) > 7 else "",
                "fold": cols[2] if stage == "b2" else (cols[3] if len(cols) > 3 else ""),
                "tmid": float(cols[4]) if len(cols) > 4 and cols[4] else 0,
                "evalue": float(cols[4]) if len(cols) > 4 and cols[4] else 999,
                "prob": float(cols[4]) if len(cols) > 4 and cols[4] else 0,
                "resolved": cols[-1] == "True" if cols else False,
            }
    return recs

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--b1", default=None)
    p.add_argument("--b2", default=None)
    p.add_argument("--b3", default=None)
    p.add_argument("--out", required=True)
    p.add_argument("--triage-stats", required=True)
    args = p.parse_args()

    b1 = read_predictions(args.b1, "b1")
    b2 = read_predictions(args.b2, "b2")
    b3 = read_predictions(args.b3, "b3")

    all_prots = set(b1) | set(b2) | set(b3)

    with open(args.out, "w") as f:
        f.write("protein\tresolved\tstage\ttool\tec\tgo\tfold\ttmid\tevalue\tprob\n")
        n_b1 = n_b2 = n_b3 = n_unres = 0
        for prot in sorted(all_prots):
            if prot in b1 and b1[prot]["resolved"]:
                info = b1[prot]; stage = "b1"; n_b1 += 1
            elif prot in b2 and b2[prot]["resolved"]:
                info = b2[prot]; stage = "b2"; n_b2 += 1
            elif prot in b3:
                info = b3[prot]; stage = "b3"; n_b3 += 1
            else:
                info = {"stage": "unresolved", "tool": "", "ec": "", "go": "",
                        "fold": "", "tmid": 0, "evalue": 999, "prob": 0}
                stage = "unresolved"; n_unres += 1

            f.write(f"{prot}\t{stage != 'unresolved'}\t{stage}\t{info['tool']}\t"
                    f"{info['ec']}\t{info['go']}\t{info['fold']}\t"
                    f"{info['tmid']:.3f}\t{info['evalue']}\t{info['prob']:.4f}\n")

    stats = {
        "structure_total": len(all_prots),
        "resolved_b1": n_b1,
        "resolved_b2": n_b2,
        "resolved_b3": n_b3,
        "unresolved_final": n_unres,
        "pct_cheap": (n_b1 + n_b2) / max(len(all_prots), 1) * 100,
    }
    with open(args.triage_stats, "w") as f:
        json.dump(stats, f, indent=2)

    print(f"Structure aggregate: B1={n_b1}, B2={n_b2}, B3={n_b3}, unresolved={n_unres}")

if __name__ == "__main__":
    main()
