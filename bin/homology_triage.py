#!/usr/bin/env python3
"""Merge homology-based annotations, apply confidence thresholds, and split
proteins into resolved (Branch A confident) vs unresolved (needs Branch B).
Proteins on the curated trait list are always routed to Branch B for cross-
validation regardless of homology confidence."""
import argparse, csv, json, sys
from collections import defaultdict

def parse_eggnog(path):
    """Return dict protein -> list of {cog, go, ec, evalue, bitscore}."""
    res = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 7:
                continue
            prot = cols[0]
            # eggNOG columns: query, seed_ortholog, evalue, score, ...
            evalue = float(cols[2]) if cols[2] != "*" else 999
            bitscore = float(cols[3]) if cols[3] != "*" else 0
            gos = cols[5] if len(cols) > 5 and cols[5] != "*" else ""
            kos = cols[6] if len(cols) > 6 and cols[6] != "*" else ""
            cog = cols[7] if len(cols) > 7 and cols[7] != "*" else ""
            ec = cols[8] if len(cols) > 8 and cols[8] != "*" else ""
            res[prot].append({
                "source": "eggnog", "evalue": evalue, "bitscore": bitscore,
                "go": gos, "ko": kos, "cog": cog, "ec": ec
            })
    return res

def parse_kofam(path, threshold_mode):
    """Return dict protein -> list of {ko, ec, score, evalue, passed}."""
    res = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 5:
                continue
            # KofamScan: * protein  KO  score  evalue  (threshold-dependent)
            if cols[0].strip() == "*":
                # Marked as above threshold
                prot = cols[1].strip()
                ko = cols[2].strip()
                score = float(cols[3]) if cols[3] != "*" else 0
                res[prot].append({"source": "kofam", "ko": ko, "score": score,
                                  "passed": True, "ec": ""})
            else:
                prot = cols[0].strip()
                ko = cols[1].strip()
                score = float(cols[2]) if cols[2] != "*" else 0
                res[prot].append({"source": "kofam", "ko": ko, "score": score,
                                  "passed": False, "ec": ""})
    return res

def parse_interpro(path):
    """Return dict protein -> list of {pfam, go, ipr}."""
    res = defaultdict(list)
    with open(path) as f:
        for line in f:
            cols = line.split("\t")
            if len(cols) < 12:
                continue
            prot = cols[0]
            pfam = cols[4] if cols[4] != "-" else ""
            ipr = cols[11] if len(cols) > 11 and cols[11] != "-" else ""
            go = cols[13] if len(cols) > 13 and cols[13] != "-" else ""
            res[prot].append({"source": "interpro", "pfam": pfam,
                              "interpro": ipr, "go": go})
    return res

def read_fasta_ids(path):
    ids = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].split()[0].strip())
    return ids

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--eggnog", required=True)
    p.add_argument("--kofam", required=True)
    p.add_argument("--interpro", required=True)
    p.add_argument("--proteins", required=True)
    p.add_argument("--evalue-cutoff", type=float, default=1e-5)
    p.add_argument("--kofam-threshold", default="score")
    p.add_argument("--trait-list", default=None)
    p.add_argument("--out-results", required=True)
    p.add_argument("--out-unresolved", required=True)
    p.add_argument("--out-trait-candidates", required=True)
    p.add_argument("--out-stats", required=True)
    args = p.parse_args()

    eggnog  = parse_eggnog(args.eggnog)
    kofam   = parse_kofam(args.kofam, args.kofam_threshold)
    interpro = parse_interpro(args.interpro)

    trait_prots = set()
    if args.trait_list and os.path.exists(args.trait_list):
        with open(args.trait_list) as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    trait_prots.add(line.strip())

    all_prots = read_fasta_ids(args.proteins)
    resolved = []
    unresolved_ids = []
    trait_faa_seqs = {}

    n_resolved_egg = 0
    n_resolved_kofam = 0
    n_resolved_interpro = 0

    # Read sequences for trait-candidates output
    seqs = {}
    with open(args.proteins) as f:
        cur = None; buf = []
        for line in f:
            if line.startswith(">"):
                if cur: seqs[cur] = "".join(buf)
                cur = line[1:].split()[0].strip(); buf = []
            else: buf.append(line.strip())
        if cur: seqs[cur] = "".join(buf)

    for prot in all_prots:
        conf = False
        sources = []
        annotation = {"protein": prot, "evidence": [], "ko": "", "ec": "",
                       "go": "", "pfam": "", "cog": "", "resolved": False}

        for hit in eggnog.get(prot, []):
            sources.append("eggnog")
            if hit["evalue"] < args.evalue_cutoff:
                conf = True
                n_resolved_egg += 1
                annotation["evalue"] = hit["evalue"]
            if hit["ko"]:  annotation["ko"]  = annotation["ko"]  + ";" + hit["ko"]
            if hit["ec"]:  annotation["ec"]  = annotation["ec"]  + ";" + hit["ec"]
            if hit["go"]:  annotation["go"]  = annotation["go"]  + ";" + hit["go"]
            if hit["cog"]: annotation["cog"] = annotation["cog"] + ";" + hit["cog"]

        for hit in kofam.get(prot, []):
            sources.append("kofam")
            if hit.get("passed") or args.kofam_threshold == "off":
                conf = True
                n_resolved_kofam += 1
            if hit["ko"]: annotation["ko"] = annotation["ko"] + ";" + hit["ko"]

        for hit in interpro.get(prot, []):
            sources.append("interpro")
            if hit["pfam"]:
                conf = True
                n_resolved_interpro += 1
            if hit["pfam"]: annotation["pfam"] = annotation["pfam"] + ";" + hit["pfam"]
            if hit["go"]:   annotation["go"]    = annotation["go"]   + ";" + hit["go"]

        annotation["evidence"] = ";".join(set(sources))

        # Trait-list proteins always go to Branch B (cross-validation)
        if prot in trait_prots:
            conf = False

        annotation["resolved"] = conf
        if conf:
            resolved.append(annotation)
        else:
            unresolved_ids.append(prot)

        if prot in trait_prots and prot in seqs:
            trait_faa_seqs[prot] = seqs[prot]

    # Write outputs
    fields = ["protein", "resolved", "evidence", "ko", "ec", "go", "pfam", "cog", "evalue"]
    with open(args.out_results, "w") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in resolved:
            w.writerow(r)
        # Include unresolved with resolved=False for completeness
        for pid in unresolved_ids:
            w.writerow({"protein": pid, "resolved": False, "evidence": "none"})

    with open(args.out_unresolved, "w") as f:
        for pid in unresolved_ids:
            if pid in seqs:
                f.write(f">{pid}\n{seqs[pid]}\n")

    with open(args.out_trait_candidates, "w") as f:
        for pid, seq in trait_faa_seqs.items():
            f.write(f">{pid}\n{seq}\n")

    stats = {
        "total_proteins": len(all_prots),
        "resolved_homology": len(resolved),
        "unresolved": len(unresolved_ids),
        "resolved_by_eggnog": n_resolved_egg,
        "resolved_by_kofam": n_resolved_kofam,
        "resolved_by_interpro": n_resolved_interpro,
    }
    with open(args.out_stats, "w") as f:
        json.dump(stats, f, indent=2)

    print(f"Homology triage: {len(resolved)} resolved, {len(unresolved_ids)} unresolved")

if __name__ == "__main__":
    main()
