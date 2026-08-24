// =============================================================================
// Module 0 — Preprocessing: gene calling, protein extraction, taxonomy
// =============================================================================

process PROKKA_ANNOTATE {
    tag        "$genome_id"
    label      'cpu_medium'
    publishDir "${params.outdir}/00_preprocessing/${genome_id}", mode: 'copy'

    input:
    tuple val(genome_id), path(genome), val(ftype)

    output:
    tuple val(genome_id), path("${genome_id}.faa"), emit: proteins
    tuple val(genome_id), path("${genome_id}.metadata.tsv"), emit: metadata
    path "${genome_id}.gff",  emit: gff,  optional: true
    path "${genome_id}.gbk",  emit: gbk,  optional: true

    script:
    def prot_prefix = "${genome_id}"
    if (ftype == 'genbank' || genome.name.endsWith('.gbk') || genome.name.endsWith('.gbff')) {
        // GenBank already annotated — extract proteins directly
        """
        # Extract protein sequences from annotated GenBank
        python ${projectDir}/bin/genbank_to_faa.py ${genome} ${genome_id}.faa ${genome_id}.metadata.tsv
        touch ${genome_id}.gbk
        """
    } else {
        // Raw FASTA → de novo annotation with Prokka
        """
        prokka --force \
               --outdir prokka_out \
               --prefix ${genome_id} \
               --cpus ${task.cpus} \
               --addgenes \
               --compliant \
               ${genome}

        cp prokka_out/${genome_id}.faa ${genome_id}.faa
        cp prokka_out/${genome_id}.gff ${genome_id}.gff
        cp prokka_out/${genome_id}.gbk ${genome_id}.gbk 2>/dev/null || true

        # Build metadata record
        printf "genome_id\\tsource\\ttaxonomy\\tn_contigs\\tn_proteins\\n" > ${genome_id}.metadata.tsv
        NCONTAGS=\$(grep -c '^>' ${genome} | head -1 || echo 0)
        NPROT=\$(grep -c '^>' ${genome_id}.faa)
        printf "${genome_id}\\t${genome}\\tunclassified\\t\${NCONTAGS}\\t\${NPROT}\\n" >> ${genome_id}.metadata.tsv
        """
    }
}

process GTDB_TAXONOMY {
    tag        "$genome_id"
    label      'cpu_high'
    publishDir "${params.outdir}/00_preprocessing/taxonomy", mode: 'copy'

    input:
    tuple path(fasta_files)

    output:
    path "gtdbtk.bac120.summary.tsv", emit: taxonomy, optional: true

    when:
    params.run_module_0

    script:
    """
    # GTDB-Tk classification (wrapper — runs only if classify step desired)
    if command -v gtdbtk >/dev/null 2>&1; then
        gtdbtk classify_wf --cpus ${task.cpus} --out_dir gtdbtk_out --genome_dir genomes_dir
        cp gtdbtk_out/gtdbtk.bac120.summary.tsv . 2>/dev/null || true
    else
        echo "GTDB-Tk not available; skipping taxonomy assignment"
    fi
    """
}

workflow PREPROCESSING {
    take:
    genomes_ch

    main:
    ann_ch = PROKKA_ANNOTATE(genomes_ch)

    emit:
    proteins = ann_ch.proteins
    metadata = ann_ch.metadata
}
