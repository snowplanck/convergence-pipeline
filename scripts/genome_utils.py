#!/usr/bin/env python3
"""Shared helper for resolving protein -> genome identity.

Genome identity must ALWAYS be resolved via the explicit
protein_to_genome.tsv mapping (produced once, during preprocessing),
never guessed from protein ID string-parsing. The only sanctioned
exception is the last-underscore heuristic used for the synthetic
positive-control proteins at the point they are first assigned a genome
(inside scripts/preprocess.py) -- everywhere else, this module's
resolve_genome() is the single source of truth.
"""
from pathlib import Path


def load_protein_to_genome(mapping_path):
    """Load protein_to_genome.tsv into a dict: protein_id -> genome_id."""
    mapping = {}
    path = Path(mapping_path)
    if not path.exists():
        return mapping
    with open(path) as f:
        header = f.readline()
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 2 and cols[0]:
                mapping[cols[0]] = cols[1]
    return mapping


def resolve_genome(protein_id, mapping, log_file=None):
    """Resolve protein_id to its genome_id using the explicit mapping.
    Falls back to the last-underscore heuristic ONLY if the protein is not
    found in the mapping, and logs a warning when that happens so any gap
    in protein_to_genome.tsv coverage is visible rather than silent."""
    genome = mapping.get(protein_id)
    if genome is not None:
        return genome

    fallback = protein_id.rsplit("_", 1)[0] if "_" in protein_id else protein_id
    msg = (f"WARNING: protein {protein_id!r} not found in protein_to_genome.tsv, "
           f"falling back to heuristic -> {fallback!r}\n")
    if log_file:
        with open(log_file, "a") as log:
            log.write(msg)
    return fallback
