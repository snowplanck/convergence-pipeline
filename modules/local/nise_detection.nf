// =============================================================================
// Module 5 — NISE / convergence detection
// =============================================================================

process NISE_DETECT {
    tag        "nise"
    label      'cpu_medium'
    publishDir "${params.outdir}/05_nise", mode: 'copy'

    input:
    path trait_matrix_tsv
    path homology_results_tsv
    path structure_results_tsv
    path orthogroups_tsv

    output:
    path "nise_candidates.tsv", emit: nise_candidates
    path "nise_summary.json",   emit: summary

    when:
    params.run_module_5

    script:
    """
    python ${projectDir}/bin/detect_nise.py \
        --trait-matrix ${trait_matrix_tsv} \
        --homology ${homology_results_tsv} \
        --structure ${structure_results_tsv} \
        --orthogroups ${orthogroups_tsv} \
        --mmseqs-evalue ${params.nise_mmseqs_evalue} \
        --hhblits-prob ${params.nise_hhblits_prob} \
        --min-tools ${params.nise_min_tools} \
        --out nise_candidates.tsv \
        --summary nise_summary.json
    """
}

workflow NISE_DETECTION {
    take:
    trait_matrix_ch, homology_results_ch, structure_results_ch, orthogroups_ch

    main:
    if (params.run_module_5) {
        NISE_DETECT(
            trait_matrix_ch.collect(),
            homology_results_ch.collect(),
            structure_results_ch.collect(),
            orthogroups_ch.collect()
        )
    } else {
        NISE_DETECT(Channel.empty(), Channel.empty(), Channel.empty(), Channel.empty())
    }

    emit:
    nise_candidates = NISE_DETECT.out.nise_candidates
}
