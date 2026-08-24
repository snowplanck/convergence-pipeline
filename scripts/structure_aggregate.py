#!/usr/bin/env python3
"""Aggregate B1/B2/B3 structure predictions into a single table."""
import time
from pathlib import Path

b1_in = snakemake.input["b1"]
b2_in = snakemake.input["b2"]
b3_in = snakemake.input["b3"]
combined_out = snakemake.output["combined"]
log_file = snakemake.log[0]
Path(combined_out).parent.mkdir(parents=True, exist_ok=True)
Path(log_file).parent.mkdir(parents=True, exist_ok=True)
t0 = time.time()


def read_preds(path):
    """Read a predictions TSV using its own header to map columns by name,
    since B1/B2/B3 each use a different column order/schema."""
    recs = {}
    if not Path(path).exists():
        return recs
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}

        def get(cols, name, default=""):
            i = idx.get(name)
            if i is None or i >= len(cols):
                return default
            return cols[i]

        for line in f:
            cols = line.rstrip("\n").split("\t")
            if not cols or not cols[0]:
                continue
            prot = cols[0]

            resolved_raw = get(cols, "resolved", None)
            # B3 has no "resolved" column: reaching B3 counts as resolved.
            resolved = (resolved_raw == "True") if resolved_raw not in (None, "") else True

            prob_raw = get(cols, "probability", get(cols, "tmid", "0"))
            try:
                prob = float(prob_raw) if prob_raw else 0.0
            except ValueError:
                prob = 0.0

            recs[prot] = {
                "tool": get(cols, "tool", ""),
                "ec": get(cols, "ec", ""),
                "fold": get(cols, "fold", ""),
                "go": get(cols, "go", ""),
                "tmid": get(cols, "tmid", ""),
                "prob": prob,
                "resolved": resolved,
            }
    return recs


b1 = read_preds(b1_in)
b2 = read_preds(b2_in)
b3 = read_preds(b3_in)
all_prots = set(b1) | set(b2) | set(b3)

with open(log_file, "w") as log:
    log.write(f"Aggregating: B1={len(b1)}, B2={len(b2)}, B3={len(b3)}, total={len(all_prots)}\n")

n_b1 = n_b2 = n_b3 = 0
with open(combined_out, "w") as f:
    f.write("protein\tresolved\tstage\ttool\tec\tgo\tfold\ttmid\tprobability\n")
    for prot in sorted(all_prots):
        if prot in b1 and b1[prot]["resolved"]:
            info = b1[prot]; stage = "b1"; n_b1 += 1
        elif prot in b2 and b2[prot]["resolved"]:
            info = b2[prot]; stage = "b2"; n_b2 += 1
        elif prot in b3:
            info = b3[prot]; stage = "b3"; n_b3 += 1
        else:
            info = {"tool": "", "ec": "", "fold": "", "go": "", "tmid": "", "prob": 0}
            stage = "unresolved"
        f.write(f"{prot}\t{stage != 'unresolved'}\t{stage}\t{info['tool']}\t"
                f"{info['ec']}\t{info.get('go','')}\t{info['fold']}\t"
                f"{info.get('tmid','')}\t{info['prob']:.4f}\n")

elapsed = time.time() - t0
with open(log_file, "a") as log:
    log.write(f"Done: B1={n_b1}, B2={n_b2}, B3={n_b3} in {elapsed:.1f}s\n")
