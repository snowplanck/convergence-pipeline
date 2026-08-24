#!/usr/bin/env python3
"""Generate minimal synthetic test genomes for pipeline dry-runs.
Creates N small FASTA files with short contigs to test workflow wiring
without real genomic data."""
import argparse, os, random

random.seed(42)

def random_seq(length):
    return "".join(random.choices("ACGT", k=length))

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n-genomes", type=int, default=5)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--out-csv", required=True)
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    rows = ["id,path,type"]

    for i in range(args.n_genomes):
        name = f"test_genome_{i:02d}"
        fasta_path = os.path.join(args.out_dir, f"{name}.fna")
        with open(fasta_path, "w") as f:
            n_contigs = random.randint(1, 3)
            for c in range(n_contigs):
                f.write(f">{name}_contig_{c}\n")
                f.write(random_seq(5000) + "\n")
        rows.append(f"{name},{fasta_path},fasta")

    with open(args.out_csv, "w") as f:
        f.write("\n".join(rows) + "\n")

    print(f"Generated {args.n_genomes} test genomes in {args.out_dir}")

if __name__ == "__main__":
    main()
