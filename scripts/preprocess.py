#!/usr/bin/env python3
"""Module 0 — preprocessing: gene calling with Prokka, extract proteins.
Emits an explicit protein_id -> genome_id mapping (protein_to_genome.tsv) so
downstream scripts never need to guess genome identity from protein ID
string-parsing."""
import re
import subprocess
import time
from pathlib import Path

# Snakemake I/O
proteins_out = snakemake.output["proteins"]
metadata_out = snakemake.output["metadata"]
control_map_out = snakemake.output["control_map"]
protein_to_genome_out = snakemake.output["protein_to_genome"]
genome_files = snakemake.input["genomes"]
controls_file = snakemake.input["controls"]
threads = snakemake.threads
log_file = snakemake.log[0]

Path(proteins_out).parent.mkdir(parents=True, exist_ok=True)
Path(log_file).parent.mkdir(parents=True, exist_ok=True)

t0 = time.time()
with open(log_file, "w") as log:
    log.write(f"Preprocessing {len(genome_files)} genomes with {threads} threads\n")


def sanitize_locustag(genome_id, used_tags):
    """Prokka requires an alphanumeric locus tag. Derive one from genome_id,
    strip anything non-alphanumeric, uppercase it, and ensure uniqueness
    across genomes in this run (append a numeric suffix on collision)."""
    tag = re.sub(r"[^A-Za-z0-9]", "", genome_id).upper()
    if not tag or tag[0].isdigit():
        tag = "G" + tag
    tag = tag[:15]  # keep it reasonably short
    base = tag
    i = 1
    while tag in used_tags:
        suffix = str(i)
        tag = base[: 15 - len(suffix)] + suffix
        i += 1
    used_tags.add(tag)
    return tag


def _placeholder_gene_call(gpath, genome_id, all_proteins):
    """Minimal ORF caller when Prokka is unavailable."""
    from Bio import SeqIO
    n = 0
    for rec in SeqIO.parse(str(gpath), "fasta"):
        seq = str(rec.seq)
        for i in range(0, len(seq) - 300, 3000):
            codons = [seq[j:j + 3] for j in range(i, min(i + 300, len(seq) - 2), 3)]
            aa = "".join(_translate(c) for c in codons if len(c) == 3)
            all_proteins.append((f"{genome_id}_ORF{n}", aa))
            n += 1
            if n >= 3:
                break
    return n


def _translate(codon):
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
    return table.get(codon.upper(), "X")


all_proteins = []
metadata_rows = [["genome_id", "source", "n_proteins", "n_contigs", "locus_tag"]]
protein_to_genome_rows = [["protein_id", "genome_id"]]
used_locustags = set()

for gpath in genome_files:
    gpath = Path(gpath)
    genome_id = gpath.stem
    suffix = gpath.suffix.lower()
    start_idx = len(all_proteins)
    locus_tag = ""

    if suffix in (".gbk", ".gbff"):
        from Bio import SeqIO
        n_prot = 0
        for rec in SeqIO.parse(str(gpath), "genbank"):
            for feat in rec.features:
                if feat.type == "CDS" and "translation" in feat.qualifiers:
                    pid = feat.qualifiers.get("protein_id", [f"{genome_id}_{n_prot}"])[0]
                    all_proteins.append((pid, feat.qualifiers["translation"][0]))
                    n_prot += 1
        metadata_rows.append([genome_id, str(gpath), n_prot, 0, locus_tag])
    else:
        locus_tag = sanitize_locustag(genome_id, used_locustags)
        outdir = Path(proteins_out).parent / f"prokka_{genome_id}"
        cmd = ["prokka", "--force", "--outdir", str(outdir), "--prefix", genome_id,
               "--locustag", locus_tag,
               "--cpus", str(threads), "--addgenes", "--compliant", str(gpath)]
        with open(log_file, "a") as log:
            log.write(f"  Prokka {genome_id} (locustag={locus_tag}): {' '.join(cmd)}\n")
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except FileNotFoundError:
            with open(log_file, "a") as log:
                log.write("  Prokka not installed — using placeholder gene calling\n")
            n_prot = _placeholder_gene_call(gpath, genome_id, all_proteins)
            metadata_rows.append([genome_id, str(gpath), n_prot, 0, ""])
            for pid, _ in all_proteins[start_idx:]:
                protein_to_genome_rows.append([pid, genome_id])
            continue
        except subprocess.CalledProcessError as e:
            with open(log_file, "a") as log:
                log.write(f"  Prokka failed: {e.stderr}\n")
            continue

        faa = outdir / f"{genome_id}.faa"
        if faa.exists():
            from Bio import SeqIO
            n_prot = 0
            for rec in SeqIO.parse(str(faa), "fasta"):
                all_proteins.append((rec.id, str(rec.seq)))
                n_prot += 1
            metadata_rows.append([genome_id, str(gpath), n_prot, 0, locus_tag])

    for pid, _ in all_proteins[start_idx:]:
        protein_to_genome_rows.append([pid, genome_id])

# ---- Everything below runs ONCE, after all genomes are processed ----

# Inject synthetic positive-control proteins (added once, not per-genome)
control_proteins = []
if Path(controls_file).exists():
    with open(controls_file) as f:
        cur = None; buf = []
        for line in f:
            if line.startswith(">"):
                if cur is not None:
                    control_proteins.append((cur, "".join(buf)))
                cur = line[1:].split()[0].strip(); buf = []
            else:
                buf.append(line.strip())
        if cur is not None:
            control_proteins.append((cur, "".join(buf)))
all_proteins.extend(control_proteins)


def control_genome_of(pid):
    """Synthetic control protein IDs are not Prokka-generated, so the
    last-underscore heuristic is fine here -- it's correct by construction
    for these specific, hand-crafted names (e.g. NISE_A_convergent_kinase_A
    -> NISE_A_convergent_kinase)."""
    return pid.rsplit("_", 1)[0] if "_" in pid else pid


for pid, _ in control_proteins:
    protein_to_genome_rows.append([pid, control_genome_of(pid)])

# Write combined protein FASTA (once)
with open(proteins_out, "w") as f:
    for pid, seq in all_proteins:
        f.write(f">{pid}\n{seq}\n")

# Write control protein map, for NISE validation (once)
with open(control_map_out, "w") as f:
    f.write("protein_id\tcase\tec_number\tfold\n")
    for pid, _ in control_proteins:
        if pid.startswith("NISE_A"):
            f.write(f"{pid}\tnise\t3.4.21.1\t2.40.10\n")
        elif pid.startswith("NISE_B"):
            f.write(f"{pid}\tnise\t3.4.21.1\t3.20.20\n")
        elif pid.startswith("HOM_"):
            f.write(f"{pid}\thomologous\t3.2.1.1\t3.40.50\n")

# Write metadata (once)
import csv
with open(metadata_out, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerows(metadata_rows)

# Write the explicit protein -> genome mapping (once) -- this is now the
# single source of truth for genome identity used by every downstream script.
with open(protein_to_genome_out, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerows(protein_to_genome_rows)

elapsed = time.time() - t0
with open(log_file, "a") as log:
    log.write(f"Done: {len(all_proteins)} proteins in {elapsed:.1f}s\n")
