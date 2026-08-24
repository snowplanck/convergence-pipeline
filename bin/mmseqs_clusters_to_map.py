#!/usr/bin/env python3
"""Convert MMseqs2 cluster output to a protein->cluster->genome mapping table
and extract representative sequences."""
import sys
from collections import defaultdict

def parse_fasta_ids(faa):
    """Return dict of seqid -> header index (assume ids are genome|prot form)."""
    ids = {}
    with open(faa) as f:
        for i, line in enumerate(f):
            if line.startswith(">"):
                sid = line[1:].split()[0].strip()
                ids[sid] = i
    return ids

def main():
    cluster_tsv = sys.argv[1]      # MMseqs2 _cluster.tsv: rep\\tmember
    out_map     = sys.argv[2]      # protein_id, cluster_rep, genome_id
    out_faa     = sys.argv[3]      # representative sequences FASTA
    input_faa   = sys.argv[4] if len(sys.argv) > 4 else "all_proteins.faa"

    # cluster_tsv: representative_id \\t member_id
    clusters = defaultdict(list)
    with open(cluster_tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            rep, mem = parts[0], parts[1]
            clusters[rep].append(mem)

    def genome_of(pid):
        return pid.rsplit("|", 1)[0] if "|" in pid else pid.rsplit("_", 1)[0]

    with open(out_map, "w") as m:
        m.write("protein_id\trepresentative_cluster_id\tgenome_id\n")
        for rep, members in clusters.items():
            for mem in members:
                m.write(f"{mem}\t{rep}\t{genome_of(mem)}\n")

    # Extract representative sequences
    # Build rep->sequence mapping by reading input FASTA
    rep_seqs = {}
    current = None
    buf = []
    with open(input_faa) as f:
        for line in f:
            if line.startswith(">"):
                if current and current in clusters:
                    rep_seqs[current] = "".join(buf)
                sid = line[1:].split()[0].strip()
                current = sid
                buf = []
            else:
                buf.append(line.strip())
        if current and current in clusters:
            rep_seqs[current] = "".join(buf)

    with open(out_faa, "w") as f:
        for rep, seq in rep_seqs.items():
            f.write(f">{rep}\n{seq}\n")

if __name__ == "__main__":
    main()
