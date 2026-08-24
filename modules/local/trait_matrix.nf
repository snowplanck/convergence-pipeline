// =============================================================================
// Module 4 — Combined trait matrix construction
// =============================================================================

process BUILD_TRAIT_MATRIX {
    tag        "trait_matrix"
    label      'cpu_low'
    publishDir "${params.outdir}/04_trait_matrix", mode: 'copy'

    input:
    path orthogroups_tsv
    path homology_results_tsv
    path structure_results_tsv
    path derep_map_tsv

    output:
    path "trait_matrix_combined.tsv", emit: trait_matrix
    path "orthogroups.tsv",           emit: orthogroups_presence
    path "trait_provenance.tsv",       emit: provenance

    when:
    params.run_module_4

    script:
    def map_file = file(derep_map_tsv).exists() ? derep_map_tsv : ''
    """
    python ${projectDir}/bin/build_trait_matrix.py \
        --orthogroups ${orthogroups_tsv} \
        --homology ${homology_results_tsv} \
        --structure ${structure_results_tsv} \
        --derep-map ${map_file} \
        --out-matrix trait_matrix_combined.tsv \
        --out-orthogroups orthogroups.tsv \
        --out-provenance trait_provenance.tsv
    """
}

workflow TRAIT_MATRIX {
    take:
    orthogroups_ch, homology_results_ch, structure_results_ch,
    derep_map_ch

    main:
    if (params.run_module_4) {
        // Derep map: if dereplication was skipped, pass empty-channel
        map_input = derep_map_ch.ifEmpty(Channel.fromPath('NO_FILE', type: 'any'))
        BUILD_TRAIT_MATRIX(
            orthogroups_ch.collect(),
            homology_results_ch.collect(),
            structure_results_ch.collect(),
            map_input.collect()
        )
    } else {
        BUILD_TRAIT_MATRIX(Channel.empty(), Channel.empty(), Channel.empty(), Channel.empty())
    }

    emit:
    trait_matrix = BUILD_TRAIT_MATRIX.out.trait_matrix
}
