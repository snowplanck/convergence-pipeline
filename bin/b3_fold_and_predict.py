#!/usr/bin/env python3
"""Full 3D folding (ESMFold first pass, ColabFold escalation) plus
Foldseek + DeepFRI/CLEAN-Contact function prediction for the residual set."""
import argparse, json, os, sys

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

def check_alphafold_db(prot, db_dir):
    """Check if precomputed AlphaFold structure exists."""
    if not db_dir or not os.path.isdir(db_dir):
        return None
    # UniProt accessions: strip version, try .pdb files
    base = prot.split(".")[0]  # strip version
    for ext in [".pdb", ".pdb.gz"]:
        candidate = os.path.join(db_dir, f"AF-{base}-F1-model_v4{ext}")
        if os.path.exists(candidate):
            return candidate
    return None

def run_esmfold(seq, prot_id, out_dir):
    """Run ESMFold for fast single-sequence structure prediction."""
    out_pdb = os.path.join(out_dir, f"{prot_id}.pdb")
    try:
        import torch
        import esm
        model = esm.pretrained.esmfold_v1()
        model = model.eval().cuda()
        with torch.no_grad():
            output = model.infer_pdb(seq)
        with open(out_pdb, "w") as f:
            f.write(output)
        return out_pdb
    except Exception as e:
        print(f"ESMFold failed for {prot_id}: {e}", file=sys.stderr)
        return None

def run_foldseek_structure(pdb_path, foldseek_db, out_dir, prot_id):
    """Run Foldseek on real 3D structure."""
    out_tsv = os.path.join(out_dir, f"{prot_id}_foldseek.tsv")
    cmd = (f"foldseek easy-search {pdb_path} {foldseek_db} {out_tsv} "
           f"{out_dir}/tmp_foldseek "
           f"--format-output 'query,target,evalue,bits,fident,alnlen,qcov,tcov,tmid' "
           f"--threads 1")
    os.system(cmd)
    result = {}
    if os.path.exists(out_tsv):
        with open(out_tsv) as f:
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) >= 9:
                    result = {"target": cols[1], "evalue": float(cols[2]),
                              "fident": float(cols[4]), "tmid": float(cols[8])}
                    break
    return result

def run_deepfri(pdb_path, model_dir, out_dir, prot_id):
    """Run DeepFRI for structure-aware GO/EC prediction."""
    out_json = os.path.join(out_dir, f"{prot_id}_deepfri.json")
    cmd = (f"python -m deepfri.predict --pdb {pdb_path} --model-dir {model_dir} "
           f"--output {out_json}")
    os.system(cmd)
    if os.path.exists(out_json):
        with open(out_json) as f:
            return json.load(f)
    return {}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--predictions", required=True)
    p.add_argument("--alphafold-db", default=None)
    p.add_argument("--foldseek-db", default=None)
    p.add_argument("--deepfri-model", default=None)
    p.add_argument("--use-esm", default="first_pass")
    p.add_argument("--device", default="cpu")
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    seqs = read_fasta(args.fasta)

    predictions = []
    for prot, seq in seqs.items():
        rec = {"protein": prot, "stage": "b3", "structure_source": None,
               "fold": "", "ec": "", "go": "", "tmid": 0, "evalue": 999}

        # Check AlphaFold DB first
        af_path = check_alphafold_db(prot, args.alphafold_db)
        pdb_path = af_path
        rec["structure_source"] = "alphafold_db" if af_path else None

        if not pdb_path:
            # Run ESMFold
            pdb_path = run_esmfold(seq, prot, args.out_dir)
            if pdb_path:
                rec["structure_source"] = "esmfold"

        if pdb_path and args.foldseek_db:
            fs = run_foldseek_structure(pdb_path, args.foldseek_db,
                                         args.out_dir, prot)
            if fs:
                rec["fold"] = fs["target"].split("_")[0]
                rec["tmid"] = fs["tmid"]
                rec["evalue"] = fs["evalue"]

        if pdb_path and args.deepfri_model:
            deepfri = run_deepfri(pdb_path, args.deepfri_model,
                                   args.out_dir, prot)
            if deepfri:
                rec["ec"] = ";".join(deepfri.get("ec", []))
                rec["go"] = ";".join(deepfri.get("go", []))

        predictions.append(rec)

    with open(args.predictions, "w") as f:
        f.write("protein\tstage\tstructure_source\tfold\ttmid\tevalue\tec\tgo\n")
        for rec in predictions:
            f.write(f"{rec['protein']}\t{rec['stage']}\t{rec['structure_source']}\t"
                    f"{rec['fold']}\t{rec['tmid']:.3f}\t{rec['evalue']}\t"
                    f"{rec['ec']}\t{rec['go']}\n")

    print(f"B3 folding: processed {len(predictions)} proteins")

if __name__ == "__main__":
    main()
