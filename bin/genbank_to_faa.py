#!/usr/bin/env python3
"""Extract protein FASTA + metadata from annotated GenBank files."""
import sys, csv
from Bio import SeqIO

def main():
    gb_path = sys.argv[1]
    faa_out = sys.argv[2]
    meta_out = sys.argv[3]

    proteins = []
    meta = {"genome_id": "", "source": gb_path, "taxonomy": "",
            "n_contigs": 0, "n_proteins": 0}
    seen = set()
    n_contigs = 0

    for rec in SeqIO.parse(gb_path, "genbank"):
        n_contigs += 1
        organism = rec.annotations.get("organism", "")
        if not meta["taxonomy"] and organism:
            meta["taxonomy"] = organism
        if not meta["genome_id"]:
            meta["genome_id"] = rec.id
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            if "translation" not in feat.qualifiers:
                continue
            loc = str(feat.location)
            if loc in seen:
                continue
            seen.add(loc)
            prot_id = feat.qualifiers.get("protein_id", [""])[0]
            if not prot_id:
                prot_id = feat.qualifiers.get("locus_tag", [""])[0]
            product = feat.qualifiers.get("product", ["unknown"])[0]
            seq = feat.qualifiers["translation"][0]
            proteins.append((prot_id, product, seq))

    meta["n_contigs"] = n_contigs
    meta["n_proteins"] = len(proteins)

    with open(faa_out, "w") as f:
        for pid, _, seq in proteins:
            f.write(f">{pid}\n{seq}\n")

    with open(meta_out, "w") as f:
        w = csv.DictWriter(f, fieldnames=["genome_id", "source", "taxonomy",
                                           "n_contigs", "n_proteins"], delimiter="\t")
        w.writeheader()
        w.writerow(meta)

if __name__ == "__main__":
    main()
