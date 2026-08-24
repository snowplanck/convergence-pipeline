#!/usr/bin/env python3
"""Generate toy genome FASTA files and positive-control proteins.
Run once before testing the pipeline."""
import random
from pathlib import Path

random.seed(42)
TOY_DIR = Path(__file__).parent.parent / "data" / "toy_genomes"
TOY_DIR.mkdir(parents=True, exist_ok=True)
BASES = "ACGT"


def random_seq(length):
    return "".join(random.choices(BASES, k=length))


def make_protein_seq(length_aa):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(random.choices(amino_acids, k=length_aa))


# NISE pair: different sequences, same EC:3.4.21.1, different folds
nise_a = make_protein_seq(300)
nise_b = make_protein_seq(310)

# Homologous pair: high identity, same OG, EC:3.2.1.1
hom_core = make_protein_seq(280)
hom_b = list(hom_core)
for i in range(0, len(hom_b), 20):
    hom_b[i] = random.choice("ACDEFGHIKLMNPQRSTVWY")
hom_b = "".join(hom_b)

# Build genomes
genomes = {
    "toy_genome_A": [("contig_1", 4000), ("contig_2", 3000)],
    "toy_genome_B": [("contig_1", 3500)],
    "toy_genome_C": [("contig_1", 4500), ("contig_2", 2500), ("contig_3", 3000)],
}

for gname, contigs in genomes.items():
    fpath = TOY_DIR / f"{gname}.fna"
    with open(fpath, "w") as f:
        for cname, clen in contigs:
            f.write(f">{gname}_{cname}\n")
            f.write(random_seq(clen) + "\n")
    print(f"Created {fpath}")

# Positive controls
control_dir = Path(__file__).parent.parent / "data" / "toy_positive_controls"
control_dir.mkdir(parents=True, exist_ok=True)
control_fasta = control_dir / "positive_controls.faa"

with open(control_fasta, "w") as f:
    f.write(">NISE_A_convergent_kinase_A\n" + nise_a + "\n")
    f.write(">NISE_B_convergent_kinase_B\n" + nise_b + "\n")
    f.write(">HOM_amylase_1\n" + hom_core + "\n")
    f.write(">HOM_amylase_2\n" + hom_b + "\n")

print(f"Created {control_fasta}")
identity = sum(a == b for a, b in zip(hom_core, hom_b)) / len(hom_core)
print(f"Homologous pair identity: {identity:.1%}")
