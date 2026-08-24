// =============================================================================
// Module 6 — Clustering comparison: function-based vs identity-based
// =============================================================================

process CLUSTER_COMPARE {
    tag        "cluster_compare"
    label      'cpu_low'
    publishDir "${params.outdir}/06_clustering", mode: 'copy'

    input:
    path trait_matrix_tsv
    path orthogroups_tsv

    output:
    path "functional_clusters.tsv",     emit: functional_clusters
    path "identity_clusters.tsv",       emit: identity_clusters
    path "clustering_comparison.md",    emit: comparison_report
    path "divergence_metrics.json",     emit: metrics

    when:
    params.run_module_6

    script:
    """
    Rscript ${projectDir}/bin/compare_clusterings.R \
        --trait-matrix ${trait_matrix_tsv} \
        --orthogroups ${orthogroups_tsv} \
        --method ${params.cluster_method} \
        --mantel-perms ${params.mantel_permutations} \
        --out-functional functional_clusters.tsv \
        --out-identity identity_clusters.tsv \
        --out-report clustering_comparison.md \
        --out-metrics divergence_metrics.json
    """
}

workflow CLUSTERING_COMPARISON {
    take:
    trait_matrix_ch, orthogroups_ch

    main:
    if (params.run_module_6) {
        CLUSTER_COMPARE(
            trait_matrix_ch.collect(),
            orthogroups_ch.collect()
        )
    } else {
        CLUSTER_COMPARE(Channel.empty(), Channel.empty())
    }

    emit:
    functional_clusters = CLUSTER_COMPARE.out.functional_clusters
    identity_clusters   = CLUSTER_COMPARE.out.identity_clusters
}
