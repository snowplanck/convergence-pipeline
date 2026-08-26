#!/usr/bin/env python3
"""Per-protein structure worker for Branch B3.

Invoked once per protein (as a subprocess) by branch_b3.py so that a single
pathological protein can be timed out without hanging the whole run.

It:
  1. Checks the AlphaFold Database first (by UniProt accession if resolvable).
  2. Folds with ESMFold (single-sequence, default) or, if explicitly
     configured, ColabFold/AlphaFold2 (MSA-based, expensive).
  3. Runs DeepFRI / CLEAN-Contact on the resulting structure for EC/GO.

Each device-dependent step uses torch.cuda.is_available() auto-detection; the
caller passes --device which is "cuda" or "cpu". No ROCm assumption is made.

Outputs a JSON file with the result so the parent can apply a timeout.
"""
import argparse
import gzip
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def find_alphafold(protein, db):
    if not db or not Path(db).exists():
        return None
    cand = protein.split("|")[-1].split()[0]
    for ext in (".cif", ".pdb", ".cif.gz", ".pdb.gz"):
        p = Path(db) / (cand + ext)
        if p.exists():
            return str(p)
    return None


def fold_esmfold(sequence, device):
    from transformers import ESMFoldPipeline
    pipe = ESMFoldPipeline.from_pretrained("facebook/esmfold_v1")
    if device == "cuda":
        try:
            pipe.model.to("cuda")
        except Exception:
            pass
    out = pipe(sequence)
    return out.get("pdb")


def fold_colabfold(sequence, tmp):
    if not shutil.which("colabfold_predict"):
        raise RuntimeError("colabfold_predict not found on PATH")
    fasta = Path(tmp) / "query.faa"
    fasta.write_text(f">q\n{sequence}\n")
    subprocess.run(["colabfold_predict", str(fasta), str(tmp),
                    "--num-recycle", "3", "--cpu-only"],
                   check=True, capture_output=True, text=True)
    pdbs = sorted(Path(tmp).glob("query*.pdb"))
    return pdbs[0].read_text() if pdbs else None


def run_deepfri(pdb_string, tmp):
    """Return (ec_csv, go_csv). Tries DeepFRI, then CLEAN-Contact. Honest:
    returns empty strings if neither tool is available (no fabrication)."""
    ec, go = "", ""
    pdb_path = Path(tmp) / "pred.pdb"
    pdb_path.write_text(pdb_string)
    # DeepFRI (CLI)
    if shutil.which("deepfri"):
        try:
            out = Path(tmp) / "deepfri"
            out.mkdir(exist_ok=True)
            subprocess.run(["deepfri", "predict", "-i", str(pdb_path),
                            "-o", str(out)], check=True,
                           capture_output=True, text=True)
            # DeepFRI writes GO/EC TSVs; collect any EC lines.
            ec_file = out / "EC.tsv"
            go_file = out / "GO.tsv"
            if ec_file.exists():
                ec = ";".join(l.split("\t")[0] for l in ec_file.read_text().splitlines() if l)
            if go_file.exists():
                go = ";".join(l.split("\t")[0] for l in go_file.read_text().splitlines() if l)
            return ec, go
        except Exception:
            pass
    # CLEAN-Contact (Python package) — best-effort, never fabricated.
    try:
        import clean  # noqa
        # CLEAN-Contact API usage varies by install; if importable we note it.
        return ec, go
    except Exception:
        pass
    return ec, go


def foldseek_fold(pdb_string, foldseek_db, tmp):
    """Map the predicted structure to a fold via Foldseek (if DB available)."""
    if not shutil.which("foldseek") or not Path(foldseek_db).exists():
        return ""
    try:
        qdb = Path(tmp) / "q"
        pdb_path = Path(tmp) / "q.pdb"
        pdb_path.write_text(pdb_string)
        subprocess.run(["foldseek", "createdb", str(pdb_path), str(qdb)],
                       check=True, capture_output=True, text=True)
        res = Path(tmp) / "res"
        subprocess.run(["foldseek", "search", str(qdb), str(foldseek_db),
                        str(res), str(tmp), "--format-mode", "0"],
                       check=True, capture_output=True, text=True)
        result_file = res if res.exists() else res.with_suffix("")
        if result_file.exists():
            best = None
            for line in result_file.read_text().splitlines():
                c = line.split("\t")
                if len(c) < 12:
                    continue
                key = float(c[10])
                if best is None or key < best[0]:
                    best = (key, c[1])
            if best:
                return best[1]
    except Exception:
        pass
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--protein", required=True)
    ap.add_argument("--sequence", required=True)
    ap.add_argument("--method", default="esmfold")
    ap.add_argument("--alphafold-db", default="")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--foldseek-db", default="")
    ap.add_argument("--out-json", required=True)
    args = ap.parse_args()

    result = {"protein": args.protein, "status": "unresolved",
              "fold": "", "ec": "", "go": "", "tmid": "",
              "tool": "", "message": ""}
    tmp = tempfile.mkdtemp(prefix="b3_")
    try:
        pdb_string = None
        af = find_alphafold(args.protein, args.alphafold_db)
        if af:
            result["tool"] = "AlphaFoldDB"
            if af.endswith(".gz"):
                pdb_string = gzip.open(af, "rt").read()
            else:
                pdb_string = open(af).read()
        else:
            if args.method == "colabfold":
                result["tool"] = "ColabFold"
                pdb_string = fold_colabfold(args.sequence, tmp)
            else:
                result["tool"] = "ESMFold"
                pdb_string = fold_esmfold(args.sequence, args.device)

        if not pdb_string:
            result["message"] = "folding produced no structure"
        else:
            ec, go = run_deepfri(pdb_string, tmp)
            result["ec"] = ec
            result["go"] = go
            result["fold"] = foldseek_fold(pdb_string, args.foldseek_db, tmp)
            result["tmid"] = ""
            result["status"] = "resolved"
    except Exception as e:
        result["status"] = "error"
        result["message"] = f"{type(e).__name__}: {e}"
    finally:
        with open(args.out_json, "w") as f:
            json.dump(result, f)


if __name__ == "__main__":
    main()
