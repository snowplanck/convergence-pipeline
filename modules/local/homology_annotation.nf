// =============================================================================
// Module 2 — Branch A: Homology-based annotation (runs on ALL representatives)
// Includes triage decision point that splits resolved vs unresolved proteins
// =============================================================================

process ORTHOFINDER {
    tag        "orthology"
    label      'cpu_high'
    publishDir "${params.outdir}/02_homology/orthofinder", mode: 'copy'

    input:
    path representative_faa

    output:
    path "OrthoFinder/Results_*/Orthogroups.tsv", emit: orthogroups
    path "OrthoFinder/Results_*/Orthogroups.GeneCountMatrix.csv", emit: gene_count

    when:
    params.run_module_2

    script:
    """
    mkdir -p of_input && cp ${representative_faa} of_input/
    orthofinder -f of_input -t ${task.cpus} -a ${task.cpus} -og
    # Standardize output name
    RESULT_DIR=\$(ls -d of_input/OrthoFinder/Results_* | head -1)
    cp \$RESULT_DIR/Orthogroups.tsv .
    cp \$RESULT_DIR/Orthogroups.GeneCountMatrix.csv . 2>/dev/null || true
    """
}

process EGGNOG_MAPPER {
    tag        "eggnog"
    label      'cpu_medium'
    publishDir "${params.outdir}/02_homology/eggnog", mode: 'copy'

    input:
    path representative_faa

    output:
    path "eggnog_annotations.tsv", emit: annotations

    when:
    params.run_module_2

    script:
    """
    mkdir -p eggnog_work
    emapper.py -i ${representative_faa} \
        --output eggnog_annotations \
        --cpu ${task.cpus} \
        --data_dir \${EGGNOG_DATA:-eggnog_work} \
        -o eggnog_out
    cp eggnog_out/emapper.annotations eggnog_annotations.tsv
    """
}

process KOFAMSCAN {
    tag        "kofam"
    label      'cpu_medium'
    publishDir "${params.outdir}/02_homology/kofam", mode: 'copy'

    input:
    path representative_faa

    output:
    path "kofamscan_output.tsv", emit: ko_hits
    path "kegg_modules_reconstructed.tsv", emit: modules

    when:
    params.run_module_2

    script:
    """
    exec_annotation -f detail-tsv \
        --cpu ${task.cpus} \
        -p \${KOFAM_PROFILE:-/opt/kofam} \
        -k \${KOFAM_KO_LIST:-/opt/kofam/ko_list} \
        --tmp-dir kofam_tmp \
        -o kofamscan_output.tsv \
        ${representative_faa}

    # Reconstruct KEGG modules (all-but-one tolerance)
    python ${projectDir}/bin/reconstruct_kegg_modules.py \
        kofamscan_output.tsv \
        kegg_modules_reconstructed.tsv \
        --tolerance 1
    """
}

process INTERPROSCAN {
    tag        "interpro"
    label      'cpu_medium'
    publishDir "${params.outdir}/02_homology/interpro", mode: 'copy'

    input:
    path representative_faa

    output:
    path "interproscan.tsv", emit: domains

    when:
    params.run_module_2

    script:
    """
    interproscan.sh -i ${representative_faa} \
        -f TSV,GFF3 \
        -goterms \
        -pa \
        -cpu ${task.cpus} \
        -dp \
        -o interproscan.tsv
    """
}

process HOMOLOGY_TRIAGE {
    tag        "triage"
    label      'cpu_low'
    publishDir "${params.outdir}/02_homology", mode: 'copy'

    input:
    path eggnog_tsv
    path kofam_tsv
    path interpro_tsv
    path representative_faa
    path trait_list

    output:
    path "homology_results.tsv",         emit: homology_results
    path "unresolved_proteins.faa",      emit: unresolved_proteins
    path "trait_candidates.faa",         emit: trait_candidates
    path "triage_stats.json",            emit: homology_stats

    when:
    params.run_module_2

    script:
    def evalue = params.homology_evalue
    def trait_opt = trait_list ? "--trait-list ${trait_list}" : ""
    """
    python ${projectDir}/bin/homology_triage.py \
        --eggnog ${eggnog_tsv} \
        --kofam ${kofam_tsv} \
        --interpro ${interpro_tsv} \
        --proteins ${representative_faa} \
        --evalue-cutoff ${evalue} \
        --kofam-threshold ${params.kofam_threshold} \
        ${trait_opt} \
        --out-results homology_results.tsv \
        --out-unresolved unresolved_proteins.faa \
        --out-trait-candidates trait_candidates.faa \
        --out-stats triage_stats.json
    """
}

workflow HOMOLOGY_ANNOTATION {
    take:
    rep_proteins_ch

    main:
    if (params.run_module_2) {
        // rep_proteins_ch is a single path (representative FASTA)
        ORTHOFINDER(rep_proteins_ch)
        EGGNOG_MAPPER(rep_proteins_ch)
        KOFAMSCAN(rep_proteins_ch)
        INTERPROSCAN(rep_proteins_ch)

        trait_ch = params.trait_list && file(params.trait_list).exists()
            ? Channel.fromPath(params.trait_list, type: 'file')
            : Channel.value('NO_TRAIT_FILE')

        HOMOLOGY_TRIAGE(
            EGGNOG_MAPPER.out.annotations,
            KOFAMSCAN.out.ko_hits,
            INTERPROSCAN.out.domains,
            rep_proteins_ch,
            trait_ch
        )
    } else {
        ORTHOFINDER(Channel.empty())
        EGGNOG_MAPPER(Channel.empty())
        KOFAMSCAN(Channel.empty())
        INTERPROSCAN(Channel.empty())
        HOMOLOGY_TRIAGE(Channel.empty(), Channel.empty(), Channel.empty(), Channel.empty(), Channel.empty())
    }

    emit:
    homology_results    = HOMOLOGY_TRIAGE.out.homology_results
    orthogroups         = ORTHOFINDER.out.orthogroups
    unresolved_proteins = HOMOLOGY_TRIAGE.out.unresolved_proteins
    trait_candidates    = HOMOLOGY_TRIAGE.out.trait_candidates
    homology_stats      = HOMOLOGY_TRIAGE.out.homology_stats
}
