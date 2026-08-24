#!/usr/bin/env python3
"""Generate triage funnel statistics report documenting cost savings
achieved by the staged architecture."""
import argparse, json, sys

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--structure-stats", default="")
    p.add_argument("--homology-stats", default="")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    # Parse structure triage stats — handle JSONL (one JSON object per line)
    struct = {}
    if args.structure_stats:
        try:
            with open(args.structure_stats) as f:
                for line in f:
                    line = line.strip()
                    if line:
                        struct.update(json.loads(line))
        except Exception:
            pass

    homology = {}
    if args.homology_stats:
        try:
            with open(args.homology_stats) as f:
                for line in f:
                    line = line.strip()
                    if line:
                        homology.update(json.loads(line))
        except Exception:
            pass

    total = homology.get("total_proteins", 0)
    res_homology = homology.get("resolved_homology", 0)
    unresolved = homology.get("unresolved", 0)
    res_b1 = struct.get("resolved_b1", 0)
    res_b2 = struct.get("resolved_b2", 0)
    res_b3 = struct.get("resolved_b3", 0)
    unres_final = struct.get("unresolved_final", 0)

    pct_homology = res_homology / max(total, 1) * 100
    pct_b1 = res_b1 / max(total, 1) * 100
    pct_b2 = res_b2 / max(total, 1) * 100
    pct_b3 = res_b3 / max(total, 1) * 100
    pct_cheap = (res_homology + res_b1 + res_b2) / max(total, 1) * 100
    pct_expensive = res_b3 / max(total, 1) * 100

    with open(args.out, "w") as f:
        f.write("stage\tprotein_count\tpct_of_total\ttier\n")
        f.write(f"total_input\t{total}\t100.00\tall\n")
        f.write(f"resolved_homology\t{res_homology}\t{pct_homology:.2f}\tcheap\n")
        f.write(f"unresolved_after_homology\t{unresolved}\t{unresolved/max(total,1)*100:.2f}\t-\n")
        f.write(f"resolved_b1_ml\t{res_b1}\t{pct_b1:.2f}\tcheap\n")
        f.write(f"resolved_b2_prostt5\t{res_b2}\t{pct_b2:.2f}\tcheap\n")
        f.write(f"resolved_b3_3d\t{res_b3}\t{pct_b3:.2f}\texpensive\n")
        f.write(f"unresolved_final\t{unres_final}\t{unres_final/max(total,1)*100:.2f}\t-\n")
        f.write(f"total_cheap\t{res_homology + res_b1 + res_b2}\t{pct_cheap:.2f}\tcheap\n")
        f.write(f"total_expensive\t{res_b3}\t{pct_expensive:.2f}\texpensive\n")

    print(f"Triage funnel: {pct_cheap:.1f}% cheap, {pct_expensive:.1f}% expensive (target <5%)")

if __name__ == "__main__":
    main()
