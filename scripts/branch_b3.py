#!/usr/bin/env python3
"""Branch B3 — full 3D structure prediction (ESMFold / ColabFold) + DeepFRI.
Input: proteins still UNRESOLVED after B2 (the smallest residual set).
In test_mode this is a fast stub that returns plausible dummy output.
This is the EXPENSIVE stage — gated by resources: gpu=1 and run with a
per-protein timeout so one pathological protein cannot hang the whole run."""
import importlib.util
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import tool_utils as tu

unresolved_in = snakemake.input["unresolved"]
predictions_out = snakemake.output["predictions"]
tmscore_cutoff = float(snakemake.params["tmscore"])
test_mode = tu.parse_bool(snakemake.params["test_mode"])
external_tools_dir = snakemake.params.get("external_tools_dir", "")
alphafold_db = snakemake.params.get("alphafold_db", "")
foldseek_db = snakemake.params.get("foldseek_db", "")
b3_method = str(snakemake.params.get("b3_method", "esmfold")).lower()
per_protein_timeout = float(snakemake.params.get("b3_per_protein_timeout", 1800.0))
outdir = snakemake.params.get("outdir", "results")
force_device = snakemake.params.get("force_device", "auto")
log_file = snakemake.log[0]

tu.prepend_tools_dir(external_tools_dir)
tool_report = Path(outdir) / "tool_availability.tsv"

Path(predictions_out).parent.mkdir(parents=True, exist_ok=True)
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
    log.write(f"B3: {len(proteins)} proteins in B2 residual — EXPENSIVE stage "
              f"(test_mode={test_mode})\n")

results = []
if test_mode:
    # Stub: assign plausible fold/EC annotations
    import random
    random.seed(44)
    folds = ["3.40.50", "2.40.10", "3.20.20"]
    for prot in proteins:
        if "NISE" in prot:
            # Same convergence signature used in B2: same EC, different fold
            # depending on which member of the synthetic pair this is.
            ec = "3.4.21.1"
            fold = "2.40.10" if "NISE_A" in prot else "3.20.20"
        else:
            fold = random.choice(folds)
            ec = random.choice(["3.4.21.4", "1.1.1.1", "3.1.3.1"])
        results.append({"protein": prot, "stage": "b3",
                        "fold": fold, "ec": ec, "go": "", "tmid": 0.75})
else:
    # ---- Real mode ----
    device = tu.detect_device(force_device)
    worker = Path(snakemake.params["worker_script"])

    # Tool availability (honest: report exactly what is and isn't present).
    af_ok, af_detail = tu.check_db([alphafold_db])
    esmfold_ok = False
    esmfold_detail = ""
    try:
        if (importlib.util.find_spec("transformers") is not None
                and importlib.util.find_spec("fair_esm") is not None):
            esmfold_ok = True
            esmfold_detail = "transformers + fair-esm available"
        else:
            esmfold_detail = "transformers and/or fair-esm not importable"
    except Exception as e:
        esmfold_detail = f"import check failed: {e}"
    colabfold_bin = tu.which("colabfold_predict")
    colabfold_ok = (b3_method == "colabfold") and bool(colabfold_bin)
    deepfri_ok = bool(tu.which("deepfri")) or (importlib.util.find_spec("clean") is not None)

    tu.record_tool(tool_report, "alphafold_db", "B3", af_ok, af_detail, "")
    tu.record_tool(tool_report, "esmfold", "B3", esmfold_ok, esmfold_detail, device)
    tu.record_tool(tool_report, "colabfold", "B3", colabfold_ok,
                   (f"binary '{colabfold_bin}' found" if colabfold_bin
                    else "colabfold_predict not found (only used if b3_method=colabfold)"),
                   device)
    tu.record_tool(tool_report, "deepfri/clean-contact", "B3", deepfri_ok,
                   ("available" if deepfri_ok else "not available (EC/GO omitted, not fabricated)"),
                   device)

    can_fold = af_ok or esmfold_ok or colabfold_ok
    with open(log_file, "a") as log:
        log.write(f"B3 device={device}; alphafold_db={af_ok}; esmfold={esmfold_ok}; "
                  f"colabfold={colabfold_ok}; deepfri={deepfri_ok}; "
                  f"method={b3_method}; per_protein_timeout={per_protein_timeout}s\n")

    if not can_fold:
        with open(log_file, "a") as log:
            log.write("B3: no folding route available (no AlphaFold DB, ESMFold, "
                      "or ColabFold). Marking all residual proteins unresolved "
                      "(honest no-hit). To enable: download AlphaFold DB, or "
                      "ensure transformers+fair-esm (ESMFold) are installed, or "
                      "set b3_method=colabfold with colabfold_predict on PATH.\n")
        tu.log_failure(log_file, "B3 folding",
                       "no AlphaFold DB, ESMFold, or ColabFold available")
    else:
        if b3_method == "colabfold" and not colabfold_ok:
            with open(log_file, "a") as log:
                log.write("B3: b3_method=colabfold but colabfold_predict missing; "
                          "falling back to ESMFold if available.\n")
        effective_method = "colabfold" if (b3_method == "colabfold" and colabfold_ok) else "esmfold"
        # Resolve a usable folding method in priority order.
        if effective_method == "colabfold" and not colabfold_ok:
            effective_method = "esmfold"

        n_resolved = 0
        n_timed_out = 0
        for prot, seq in seqs.items():
            out_json = Path(tempfile.gettempdir()) / f"b3_{prot}.json"
            cmd = [sys.executable, str(worker),
                   "--protein", prot, "--sequence", seq,
                   "--method", effective_method,
                   "--alphafold-db", str(alphafold_db),
                   "--foldseek-db", str(foldseek_db),
                   "--device", device,
                   "--out-json", str(out_json)]
            timed_out = False
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True,
                               timeout=per_protein_timeout)
            except subprocess.TimeoutExpired:
                timed_out = True
                with open(log_file, "a") as log:
                    log.write(f"B3: protein {prot} TIMED OUT after "
                              f"{per_protein_timeout}s (marked unresolved, continuing)\n")
                n_timed_out += 1
                if out_json.exists():
                    out_json.unlink()
                continue
            except Exception as e:
                with open(log_file, "a") as log:
                    log.write(f"B3: protein {prot} worker error: {e}\n")
            if out_json.exists():
                try:
                    res = json.loads(out_json.read_text())
                except Exception:
                    res = {"status": "error"}
                out_json.unlink()
                if res.get("status") == "resolved":
                    n_resolved += 1
                    results.append({"protein": prot, "stage": "b3",
                                    "fold": res.get("fold", ""),
                                    "ec": res.get("ec", ""),
                                    "go": res.get("go", ""),
                                    "tmid": res.get("tmid", "")})
                else:
                    with open(log_file, "a") as log:
                        log.write(f"B3: protein {prot} unresolved "
                                  f"({res.get('message','')})\n")
        with open(log_file, "a") as log:
            log.write(f"B3: resolved {n_resolved} proteins, timed out {n_timed_out}\n")

with open(predictions_out, "w") as f:
    f.write("protein\tstage\tfold\tec\ttmid\tgo\n")
    for r in results:
        f.write(f"{r['protein']}\t{r['stage']}\t{r['fold']}\t{r['ec']}\t"
                f"{r['tmid']}\t{r.get('go','')}\n")

elapsed = time.time() - t0
with open(log_file, "a") as log:
    log.write(f"B3 done: {len(results)} proteins processed in {elapsed:.1f}s\n")
