#!/usr/bin/env python3
"""Branch B1 — sequence-only ML (ProteInfer / CLEAN).
Input: unresolved proteins from homology triage + trait candidates (residual set).
Output: predictions table, resolved FASTA, still-unresolved FASTA.
Only processes the residual, NOT the full protein set."""
import os
import subprocess
import tempfile
import time
from pathlib import Path

import tool_utils as tu

unresolved_in = snakemake.input["unresolved"]
trait_candidates_in = snakemake.input["trait_candidates"]
predictions_out = snakemake.output["predictions"]
resolved_out = snakemake.output["resolved"]
still_unresolved_out = snakemake.output["still_unresolved"]
confidence = float(snakemake.params["confidence"])
test_mode = tu.parse_bool(snakemake.params["test_mode"])
external_tools_dir = snakemake.params.get("external_tools_dir", "")
b1_method = str(snakemake.params.get("b1_method", "auto")).lower()
outdir = snakemake.params.get("outdir", "results")
log_file = snakemake.log[0]

tu.prepend_tools_dir(external_tools_dir)
tool_report = Path(outdir) / "tool_availability.tsv"

for p in [predictions_out, resolved_out, still_unresolved_out]:
    Path(p).parent.mkdir(parents=True, exist_ok=True)
Path(log_file).parent.mkdir(parents=True, exist_ok=True)
t0 = time.time()

# Merge unresolved + trait candidates
def read_fasta(path):
    seqs = {}
    if not Path(path).exists():
        return seqs
    with open(path) as f:
        cur = None; buf = []
        for line in f:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                cur = line[1:].split()[0].strip(); buf = []
            else:
                buf.append(line.strip())
        if cur is not None:
            seqs[cur] = "".join(buf)
    return seqs

seqs = read_fasta(unresolved_in)
seqs.update(read_fasta(trait_candidates_in))
proteins = list(seqs.keys())

with open(log_file, "w") as log:
    log.write(f"B1: {len(proteins)} proteins in residual set (test_mode={test_mode})\n")

# In test mode, simulate: resolve ~70% with high confidence, leave rest unresolved
resolved = {}
unresolved = []
if test_mode:
    import random
    random.seed(42)
    for prot in proteins:
        if "NISE" in prot:
            # Synthetic non-homologous controls must always fall through to
            # the structural branch (B2), never get resolved here by chance.
            unresolved.append(prot)
            continue
        # Simulate confidence score
        score = random.uniform(0.7, 1.0) if random.random() < 0.7 else random.uniform(0.3, 0.89)
        if score >= confidence:
            # Assign a plausible EC based on protein length to ensure variety
            ec = random.choice(["3.4.21.4", "3.2.1.1", "1.1.1.1", "3.1.3.1", "2.7.7.7"])
            resolved[prot] = {"tool": "ProteInfer", "ec": ec, "prob": score}
        else:
            unresolved.append(prot)
else:
    # ---- Real mode: ProteInfer / CLEAN (sequence-only EC/GO prediction) ----
    device = tu.detect_device(snakemake.params.get("force_device", "auto"))
    # Resolve which tool to call.
    # Priority: explicit b1_method -> locate binary on PATH (or external_tools_dir).
    candidate_cmds = []
    if b1_method in ("proteinfer", "proteinfier", "auto"):
        candidate_cmds.append("proteinfer")
    if b1_method in ("clean", "clean_predict", "auto"):
        candidate_cmds.append("clean_predict")
    tool_bin = None
    for cmd in candidate_cmds:
        if tu.which(cmd):
            tool_bin = cmd
            break

    ok = bool(tool_bin)
    detail = (f"binary '{tool_bin}' found" if tool_bin
              else "no ProteInfer/CLEAN binary found on PATH "
                   "(set config b1_method + install the tool, or point "
                   "external_tools_dir at its location)")
    tu.record_tool(tool_report, "proteinfer/clean", "B1", ok, detail, device)
    with open(log_file, "a") as log:
        log.write(f"B1 real-mode tool: {tool_bin or 'NONE'} (device={device}); {detail}\n")

    if tool_bin:
        # Write residual FASTA and invoke the tool. The tool is expected to read
        # the FASTA and emit a TSV (protein<TAB>ec<TAB>probability[<TAB>go])
        # to the path given by --out / stdout.
        with tempfile.TemporaryDirectory() as td:
            fasta_path = Path(td) / "b1_input.faa"
            with open(fasta_path, "w") as f:
                for prot, seq in seqs.items():
                    f.write(f">{prot}\n{seq}\n")
            out_path = Path(td) / "b1_output.tsv"
            cmd = [tool_bin, "--input", str(fasta_path), "--out", str(out_path),
                   "--device", device, "--confidence", str(confidence)]
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                # Parse the tool's TSV output.
                if out_path.exists():
                    with open(out_path) as f:
                        for line in f:
                            cols = line.rstrip("\n").split("\t")
                            if not cols or not cols[0] or cols[0].startswith("#"):
                                continue
                            prot = cols[0]
                            ec = cols[1] if len(cols) > 1 else ""
                            prob = float(cols[2]) if len(cols) > 2 and cols[2] else 0.0
                            if prot in seqs and ec and prob >= confidence:
                                resolved[prot] = {"tool": tool_bin,
                                                  "ec": ec, "prob": prob}
            except Exception as e:
                with open(log_file, "a") as log:
                    log.write(f"B1 tool '{tool_bin}' failed ({e}); "
                              f"marking ALL residual proteins unresolved\n")
                tu.log_failure(log_file, "ProteInfer/CLEAN",
                               f"tool errored: {e}")
    else:
        with open(log_file, "a") as log:
            log.write("B1: no sequence-only ML tool available; "
                      "marking residual proteins unresolved (honest no-hit)\n")
        tu.log_failure(log_file, "ProteInfer/CLEAN", detail)

    unresolved = [p for p in proteins if p not in resolved]

with open(predictions_out, "w") as f:
    f.write("protein\ttool\tec\tprobability\tresolved\n")
    for prot, info in resolved.items():
        f.write(f"{prot}\t{info['tool']}\t{info.get('ec','')}\t{info['prob']:.4f}\tTrue\n")
    for prot in unresolved:
        f.write(f"{prot}\tnone\t\t0.0\tFalse\n")

with open(resolved_out, "w") as f:
    for prot in resolved:
        f.write(f">{prot}\n{seqs[prot]}\n")
with open(still_unresolved_out, "w") as f:
    for prot in unresolved:
        if prot in seqs:
            f.write(f">{prot}\n{seqs[prot]}\n")

elapsed = time.time() - t0
with open(log_file, "a") as log:
    log.write(f"B1 done: {len(resolved)} resolved, {len(unresolved)} unresolved in {elapsed:.1f}s\n")
