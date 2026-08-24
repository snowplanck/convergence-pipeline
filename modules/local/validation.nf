// =============================================================================
// Module 7 — Validation (positive controls + curated DB cross-check)
// =============================================================================

process VALIDATE {
    tag        "validate"
    label      'cpu_low'
    publishDir "${params.outdir}/07_validation", mode: 'copy'

    input:
    path nise_candidates_tsv
    path structure_results_tsv
    path trait_matrix_tsv

    output:
    path "validation_report.md",     emit: report
    path "control_recovery.tsv",     emit: control_recovery

    when:
    params.run_module_7

    script:
    """
    python ${projectDir}/bin/validate.py \
        --nise ${nise_candidates_tsv} \
        --structure ${structure_results_tsv} \
        --trait-matrix ${trait_matrix_tsv} \
        --positive-controls ${projectDir}/assets/positive_controls.json \
        --swissprot \${SWISSPROT_FASTA:-/db/swissprot/uniprot_sprot.fasta} \
        --out-report validation_report.md \
        --out-recovery control_recovery.tsv
    """
}

workflow VALIDATION {
    take:
    nise_ch, struct_ch, trait_matrix_ch

    main:
    if (params.run_module_7) {
        VALIDATE(
            nise_ch.collect(),
            struct_ch.collect(),
            trait_matrix_ch.collect()
        )
    } else {
        VALIDATE(Channel.empty(), Channel.empty(), Channel.empty())
    }

    emit:
    report = VALIDATE.out.report
}
