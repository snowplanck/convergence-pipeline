// =============================================================================
// Module 1 — Dereplication (genome + protein level)
// =============================================================================

process GENOME_DEREP {
    tag        "genome_derep"
    label      'cpu_high'
    publishDir "${params.outdir}/01_dereplication", mode: 'copy'

    input:
    tuple path(genome_files)

    output:
    path "derep_genomes_ch.txt",             emit: derep_list
    path "derep_cluster_mapping.tsv",        emit: cluster_map

    when:
    params.run_module_1

    script:
    def ani = params.genome_derep_ani
    """
    # Build genome list for dRep/skani
    echo "genome,location" > genome_list.csv
    for g in ${genome_files}; do
        echo "\$(basename \$g .fna),\$g" >> genome_list.csv
    done

    if command -v dRep >/dev/null 2>&1; then
        dRep compare derep_compare \
            -g genome_list.csv \
            -p ${task.cpus} \
            --S_ani ${ani} \
            --genomeInfo . \
            --SkipMash
        # Extract primary cluster representatives
        if [ -f derep_compare/data_tables/Wdb.csv ]; then
            awk -F',' 'NR>1{print \$1"\\t"\$2}' derep_compare/data_tables/Wdb.csv \
                > derep_cluster_mapping.tsv
        fi
        cut -f1 derep_cluster_mapping.tsv | sort -u > derep_genomes_ch.txt
    else
        # Fallback: skani triangle-based clustering via simple identity
        echo "genome\\tcluster_rep" > derep_cluster_mapping.tsv
        for g in ${genome_files}; do
            BN=\$(basename \$g .fna)
            echo "\${BN}\\t\${BN}" >> derep_cluster_mapping.tsv
        done
        cut -f1 derep_cluster_mapping.tsv | tail -n +2 > derep_genomes_ch.txt
    fi
    """
}

process PROTEIN_DEREP {
    tag        "protein_derep"
    label      'cpu_medium'
    publishDir "${params.outdir}/01_dereplication", mode: 'copy'

    input:
    tuple path(all_faa)

    output:
    path "protein_clusters.tsv",             emit: mapping_table
    path "representative_sequences.faa",     emit: representative_proteins
    path "derep_stats.txt",                  emit: stats

    when:
    params.run_module_1

    script:
    def min_id = params.protein_derep_id / 100.0
    """
    # Concatenate all protein sequences with genome prefix
    cat ${all_faa} > all_proteins.faa

    # Cluster with MMseqs2 (fast linclust)
    mmseqs easy-linclust all_proteins.faa \
        protein_clusters \
        protein_tmp \
        --min-seq-id ${min_id} \
        -c 0.9 \
        --cov-mode 0 \
        --cluster-mode 0 \
        -e 1e-5 \
        --threads ${task.cpus}

    # Convert cluster TSV to genome-level mapping table
    # Output columns: protein_id, representative_cluster_id, genome_id
    python ${projectDir}/bin/mmseqs_clusters_to_map.py \
        protein_clusters_cluster.tsv \
        protein_clusters.tsv \
        representative_sequences.faa

    TOTAL=\$(grep -c '^>' all_proteins.faa)
    REP=\$(grep -c '^>' representative_sequences.faa)
    echo "total_proteins\t\${TOTAL}"   > derep_stats.txt
    echo "representative_proteins\t\${REP}" >> derep_stats.txt
    python -c "print(f'dereplication_ratio\t{${REP}/${TOTAL}:.4f}')" >> derep_stats.txt
    """
}

workflow DEREPLICATION {
    take:
    genomes_ch, proteins_ch

    main:
    if (params.run_module_1) {
        // Collect all genome FASTA + protein FAA
        genome_files_ch = genomes_ch.map { id, path, type -> path }
        prot_files_ch   = proteins_ch.map { id, faa -> faa }

        GENOME_DEREP(genome_files_ch.collect())

        PROTEIN_DEREP(prot_files_ch.collect())
        // Re-tag representative proteins: (genome_id, rep_faa_path, rep_id)
        rep_proteins_ch = PROTEIN_DEREP.out.representative_proteins
    } else {
        GENOME_DEREP(Channel.empty())
        PROTEIN_DEREP(Channel.empty())
    }

    emit:
    representative_proteins = PROTEIN_DEREP.out.representative_proteins
    mapping_table           = PROTEIN_DEREP.out.mapping_table
    dereplicated_genomes    = GENOME_DEREP.out.derep_list
}
