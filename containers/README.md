# =============================================================================
# Container / environment definitions
# Each tool stack gets its own environment to avoid dependency conflicts
# between ML-heavy packages (ProstT5, DeepFRI, CLEAN, ESMFold).
# =============================================================================

# --- 1. Cheap annotation stack (Modules 0, 1, 2) ---
# File: containers/annotation.Dockerfile
# Base: biocontainers/ubuntu + prokka, mmseqs2, eggnog-mapper,
#       kofamscan, interproscan, orthofinder

FROM quay.io/biocontainers/base:2.18

RUN apt-get update && apt-get install -y \
    openjdk-11-jre python3 perl libdb-perl && \
    rm -rf /var/lib/apt/lists/*

# MMseqs2
RUN wget -q https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz && \
    tar xzf mmseqs-linux-avx2.tar.gz -C /opt && \
    ln -s /opt/mmseqs/bin/mmseqs /usr/local/bin/mmseqs

# eggNOG-mapper
RUN pip install eggnog-mapper

# KofamScan (requires ko_list + profiles from KEGG)
# Volumes: KOFAM_KO_LIST, KOFAM_PROFILE

# OrthoFinder
RUN pip install orthofinder

ENTRYPOINT ["/bin/bash"]
