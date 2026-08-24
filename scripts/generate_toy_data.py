#!/usr/bin/env python3
"""Generate toy dataset for pipeline validation.
Creates 3 small bacterial genome FASTA files plus two synthetic positive-control
protein cases injected into the protein set after preprocessing.

Case A (NISE): two proteins with NO sequence homology, different Foldseek fold
assignments, but the SAME assigned EC number (EC:3.4.21.1). MUST appear in
nise_candidates.tsv.

Case B (homologous): two clearly homologous proteins (same OG) sharing EC:3.2.1.1.
MUST NOT appear in nise_candidates.tsv.
"""
import random
from pathlib import Path

random.seed(42)
TOY_DIR = Path("data/toy_genomes")
TOY_DIR.mkdir(parents=True, exist_ok=True)

BASES = "ACGT"


def random_seq(length):
    return "".join(random.choices(BASES, k=length))


def translate(dna):
    table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }
    aa = []
    for i in range(0, len(dna) - 2, 3):
        aa.append(table.get(dna[i:i+3].upper(), "X"))
    return "".join(aa).rstrip("*")


def make_protein_seq(length_aa):
    """Generate a random protein sequence of given amino-acid length."""
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(random.choices(amino_acids, k=length_aa))


# ---- Case A: NISE pair (convergent, different folds, same EC:3.2.1.1 → use 3.4.21.1) ----
# Two proteins with completely different sequences (no homology) but same function.
nise_a = make_protein_seq(300)  # "convergent_fold_A" — assigned fold 2.40.10
nise_b = make_protein_seq(310)  # "convergent_fold_B" — assigned fold 3.20.20
# Both assigned EC:3.4.21.1

# ---- Case B: Homologous pair (same OG, same function EC:3.2.1.1) ----
# Two proteins with high sequence homology (share >90% identity) — same ancestry.
hom_core = make_protein_seq(280)
# hom_b is a slight variant of hom_core (same OG)
hom_b = list(hom_core)
for i in range(0, len(hom_b), 20):  # ~5% divergence
    hom_b[i] = random.choice("ACDEFGHIKLMNPQRSTVWY")
hom_b = "".join(hom_b)

# ---- Build genome FASTA files ----
genomes = {
    "toy_genome_A": [
        ("contig_1", 4000),
        ("contig_2", 3000),
    ],
    "toy_genome_B": [
        ("contig_1", 3500),
    ],
    "toy_genome_C": [
        ("contig_1", 4500),
        ("contig_2", 2500),
        ("contig_3", 3000),
    ],
}

for gname, contigs in genomes.items():
    fpath = TOY_DIR / f"{gname}.fna"
    with open(fpath, "w") as f:
        for cname, clen in contigs:
            f.write(f">{gname}_{cname}\n")
            f.write(random_seq(clen) + "\n")
    print(f"Created {fpath}")

# ---- Write synthetic positive-control proteins to a separate FASTA ----
# These get injected into the protein pool to guarantee the test cases exist.
control_dir = Path("data/toy_positive_controls")
control_dir.mkdir(parents=True, exist_ok=True)
control_fasta = control_dir / "positive_controls.faa"

with open(control_fasta, "w") as f:
    # Case A: NISE — different sequences, same EC:3.4.21.1, different folds
    f.write(">NISE_A_convergent_kinase_A\n" + nise_a + "\n")
    f.write(">NISE_B_convergent_kinase_B\n" + nise_b + "\n")
    # Case B: Homologous — same OG, EC:3.2.1.1
    f.write(">HOM_amylase_1\n" + hom_core + "\n")
    f.write(">HOM_amylase_2\n" + hom_b + "\n")

print(f"Created {control_fasta}")
print(f"NISE pair: {len(nise_a)}aa vs {len(nise_b)}aa (no homology, same EC:3.4.21.1)")
print(f"Homologous pair: {len(hom_core)}aa vs {len(hom_b)}aa (same OG, EC:3.2.1.1)")
print(f"Identity between homologous pair: {sum(a==b for a,b in zip(hom_core, hom_b))/len(hom_core):.1%}")
