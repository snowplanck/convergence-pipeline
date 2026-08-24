// =============================================================================
// Triage funnel statistics report
// =============================================================================

process TRIAGE_SUMMARY {
    label 'cpu_low'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path structure_stats
    path homology_stats

    output:
    path "pipeline_run_stats.tsv", emit: stats

    script:
    """
    python ${projectDir}/bin/triage_report.py \
        --structure-stats ${structure_stats} \
        --homology-stats ${homology_stats} \
        --out pipeline_run_stats.tsv
    """
}

workflow TRIAGE_REPORT {
    take:
    triage_stats_ch, homology_stats_ch

    main:
    // Collect all triage JSON blobs from structure module into one file
    collected_stats = triage_stats_ch.collectFile(name: 'structure_triage.jsonl', newLine: true)
    collected_homology = homology_stats_ch.collectFile(name: 'homology_stats.jsonl', newLine: true)
    TRIAGE_SUMMARY(collected_stats, collected_homology)

    emit:
    stats = TRIAGE_SUMMARY.out.stats
}
