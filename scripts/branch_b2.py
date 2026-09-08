#!/usr/bin/env python3
"""Branch B2 — ProstT5 3Di prediction + Foldseek structural search.
Input: proteins still UNRESOLVED after B1 (residual set only).
Output: fold/EC predictions, resolved FASTA, still-unresolved FASTA.
This is the key cost optimization: ProstT5 -> 3Di is ~1000x faster than
full 3D folding. ProstT5 uses GPU when available (auto-detected); the
Foldseek search step is CPU-based regardless."""
import os
import subprocess
import tempfile
import time
from pathlib import Path

import tool_utils as tu

unresolved_in = snakemake.input["unresolved"]
predictions_out = snakemake.output["predictions"]
resolved_out = snakemake.output["resolved"]
still_unresolved_out = snakemake.output["still_unresolved"]
tmscore_cutoff = float(snakemake.params["tmscore"])
test_mode = tu.parse_bool(snakemake.params["test_mode"])
external_tools_dir = snakemake.params.get("external_tools_dir", "")
foldseek_db = snakemake.params.get("foldseek_db", "")
outdir = snakemake.params.get("outdir", "results")
log_file = snakemake.log[0]

tu.prepend_tools_dir(external_tools_dir)
tool_report = Path(outdir) / "tool_availability.tsv"

for p in [predictions_out, resolved_out, still_unresolved_out]:
    Path(p).parent.mkdir(parents=True, exist_ok=True)
Path(log_file).parent.mkdir(parents=True, exist_ok=True)
t0 = time.time()

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
proteins = list(seqs.keys())

with open(log_file, "w") as log:
    log.write(f"B2: {len(proteins)} proteins in B1 residual (test_mode={test_mode})\n")

resolved = {}
unresolved = []

if test_mode:
    # Simulate: resolve ~80% of this smaller residual, leave ~20% for B3
    import random
    random.seed(43)
    # Simulate different folds — crucial for NISE detection
    folds = ["3.40.50", "2.40.10", "3.20.20", "3.90.180", "1.10.600"]
    for prot in proteins:
        score = random.uniform(0.4, 1.0)
        fold = random.choice(folds)
        # NISE positive-control proteins get different folds but same EC
        if "NISE" in prot:
            ec = "3.4.21.1"
            fold = "2.40.10" if "NISE_A" in prot else "3.20.20"
        else:
            ec = random.choice(["3.4.21.4", "3.2.1.1", "1.1.1.1"])
        if score >= tmscore_cutoff:
            resolved[prot] = {"fold": fold, "ec": ec, "tmid": score, "evalue": 1e-8}
        else:
            unresolved.append(prot)
else:
    # ---- Real mode: ProstT5 -> 3Di -> Foldseek ----
    device = tu.detect_device(snakemake.params.get("force_device", "auto"))

    # 1) Tool/binary availability
    foldseek_bin = tu.which("foldseek")
    fs_ok = bool(foldseek_bin)
    fs_detail = (f"binary '{foldseek_bin}' found"
                 if foldseek_bin else "foldseek binary not found on PATH")
    tu.record_tool(tool_report, "foldseek", "B2", fs_ok, fs_detail, "cpu")

    # 2) Foldseek reference database availability (must be downloaded locally)
    db_ok, db_detail = tu.check_db([foldseek_db])
    if foldseek_db and db_ok:
        # A foldseek DB directory contains files named <db>.* (e.g. .dbtype/.index).
        db_files = list(Path(foldseek_db).glob("*")) if Path(foldseek_db).is_dir() else [Path(foldseek_db)]
        db_ok = any(Path(foldseek_db).with_suffix(s) for s in ("",)) or bool(db_files)
        if not db_files:
            db_ok = False
            db_detail = f"foldseek DB path exists but is empty: {foldseek_db}"
    tu.record_tool(tool_report, "foldseek_db", "B2", db_ok, db_detail, "cpu")

    # 3) ProstT5 (GPU when available) availability
    prostt5_ok = False
    prostt5_detail = ""
    model = None
    if tu.which("prostt5"):
        prostt5_ok = True
        prostt5_detail = "prostt5 CLI found"
    else:
        try:
            from transformers import AutoTokenizer, AutoModelForSeq2SeqLM  # noqa
            prostt5_ok = True
            prostt5_detail = "transformers available; will load Rostlab/prostt5"
        except Exception as e:
            prostt5_detail = f"transformers/prostt5 unavailable: {e}"
    tu.record_tool(tool_report, "prostt5", "B2", prostt5_ok, prostt5_detail, device)

    with open(log_file, "a") as log:
        log.write(f"B2 device={device}; prostt5={prostt5_ok}; "
                  f"foldseek={fs_ok}; foldseek_db={db_ok}\n")

    if not fs_ok or not db_ok:
        # Honestly fail: missing binary or (critically) missing reference DB.
        reason = []
        if not fs_ok:
            reason.append("foldseek binary missing")
        if not db_ok:
            reason.append(
                f"Foldseek reference database missing at '{foldseek_db}'. "
                f"Download it with, e.g.: "
                f"foldseek databases PDB databases/foldseek/pdb tmp  (or CATH/SCOP/AlphaFold)")
        reason = "; ".join(reason)
        with open(log_file, "a") as log:
            log.write(f"B2: cannot run structural search — {reason}. "
                      f"Marking all residual proteins unresolved (honest no-hit).\n")
        tu.log_failure(log_file, "Foldseek/ProstT5", reason)
        unresolved = list(proteins)
    else:
        # Genuine run: translate AA -> 3Di with ProstT5, then Foldseek search.
        try:
            with tempfile.TemporaryDirectory() as td:
                td = Path(td)
                # Translate each protein to 3Di using ProstT5 (GPU if present).
                three_di = {}
                if tu.which("prostt5"):
                    faa = td / "query.faa"
                    with open(faa, "w") as f:
                        for p, s in seqs.items():
                            f.write(f">{p}\n{s}\n")
                    out3di = td / "query.3di.faa"
                    subprocess.run([tu.which("prostt5"), "--input", str(faa),
                                    "--output", str(out3di), "--device", device],
                                   check=True, capture_output=True, text=True)
                    for p, s in read_fasta(out3di).items():
                        three_di[p] = s
                else:
                    # transformers-based ProstT5
                    import torch
                    from transformers import AutoTokenizer, AutoModelForSeq2SeqLM
                    tok = AutoTokenizer.from_pretrained("Rostlab/prostt5", use_fast=False)
                    # use_fast=False: the fast (Rust) tokenizer implementation
                    # mis-detects this SentencePiece model's algorithm on some
                    # transformers/tokenizers version combinations, raising
                    # "trying to run a Unigram model but trained with a
                    # different algorithm". The slow (pure Python) tokenizer
                    # reads the same file correctly.
                    model = AutoModelForSeq2SeqLM.from_pretrained("Rostlab/prostt5")
                    dev = torch.device("cuda" if device == "cuda" else "cpu")
                    model = model.to(dev)
                    for p, s in seqs.items():
                        inp = tok(f">>3Di</>{s}", return_tensors="pt",
                                  add_special_tokens=False).to(dev)
                        out = model.generate(**inp, max_length=2000)
                        three_di[p] = tok.batch_decode(out, skip_special_tokens=True)[0]

                # Write 3Di FASTA and search with Foldseek.
                three_di_faa = td / "query.3di.faa"
                with open(three_di_faa, "w") as f:
                    for p, s in three_di.items():
                        f.write(f">{p}\n{s}\n")
                qdb = td / "query_3di"
                res = td / "result.m8"
                subprocess.run([foldseek_bin, "createdb", str(three_di_faa),
                                str(qdb)], check=True, capture_output=True, text=True)
                subprocess.run([foldseek_bin, "search", str(qdb),
                                str(foldseek_db), str(res), str(td / "tmp"),
                                "--format-mode", "0"],
                               check=True, capture_output=True, text=True)
                # Parse best hit per query from the tabular result.
                best = {}
                if res.with_suffix("").exists() or res.exists():
                    # foldseek writes <res> (no suffix) as the m8 when format-mode 0
                    result_file = res if res.exists() else res.with_suffix("")
                    with open(result_file) as f:
                        for line in f:
                            c = line.rstrip("\n").split("\t")
                            if len(c) < 12:
                                continue
                            q, t, fident, evalue = c[0], c[1], c[2], c[10]
                            key = float(evalue)
                            if q not in best or key < best[q][0]:
                                best[q] = (key, t, float(fident), float(evalue))
                for p in proteins:
                    if p in best:
                        _, target, fident, evalue = best[p]
                        # Use the matched target's domain id as the fold;
                        # a matching structure is itself strong evidence.
                        fold = target
                        tmid = fident
                        # EC annotation from the matched structure is DB-specific;
                        # left empty unless an annotation side-car is supplied.
                        ec = ""
                        if tmid >= tmscore_cutoff or float(evalue) <= 1e-3:
                            resolved[p] = {"fold": fold, "ec": ec,
                                           "tmid": tmid, "evalue": evalue}
            with open(log_file, "a") as log:
                log.write(f"B2: ProstT5+Foldseek resolved {len(resolved)} proteins\n")
        except Exception as e:
            with open(log_file, "a") as log:
                log.write(f"B2: ProstT5/Foldseek run failed ({e}); "
                          f"marking all residual proteins unresolved\n")
            tu.log_failure(log_file, "ProstT5/Foldseek", f"tool errored: {e}")
            unresolved = list(proteins)

    unresolved = [p for p in proteins if p not in resolved]

with open(predictions_out, "w") as f:
    f.write("protein\tfold\tevalue\ttmid\tec\tresolved\n")
    for prot, info in resolved.items():
        f.write(f"{prot}\t{info['fold']}\t{info['evalue']}\t{info['tmid']:.3f}\t"
                f"{info.get('ec','')}\tTrue\n")
    for prot in unresolved:
        f.write(f"{prot}\tnone\t\t\t\tFalse\n")

with open(resolved_out, "w") as f:
    for prot in resolved:
        f.write(f">{prot}\n{seqs[prot]}\n")
with open(still_unresolved_out, "w") as f:
    for prot in unresolved:
        if prot in seqs:
            f.write(f">{prot}\n{seqs[prot]}\n")

elapsed = time.time() - t0
with open(log_file, "a") as log:
    log.write(f"B2 done: {len(resolved)} resolved, {len(unresolved)} residual for B3 in {elapsed:.1f}s\n")
