// =============================================================================
// Main orchestration — Comparative Genomics Convergence Pipeline
// =============================================================================
include { PREPROCESSING            } from './modules/local/preprocessing'
include { DEREPLICATION            } from './modules/local/dereplication'
include { HOMOLOGY_ANNOTATION      } from './modules/local/homology_annotation'
include { STRUCTURE_ANNOTATION     } from './modules/local/structure_annotation'
include { TRAIT_MATRIX             } from './modules/local/trait_matrix'
include { NISE_DETECTION           } from './modules/local/nise_detection'
include { CLUSTERING_COMPARISON    } from './modules/local/clustering_comparison'
include { VALIDATION               } from './modules/local/validation'
include { TRIAGE_REPORT            } from './modules/local/triage_report'

workflow CONVERGENCE_PIPELINE {

    // ---- Channel: input genomes ----
    genomes_ch = Channel
        .fromPath(params.input_genomes)
        .splitCsv(header: true)
        .map { row -> tuple(row.id, file(row.path), row.type) }
        .ifEmpty { error "No genomes found in ${params.input_genomes}" }

    // ---- Module 0: preprocessing (gene calling) ----
    if (params.run_module_0) {
        PREPROCESSING(genomes_ch)
        proteins_ch   = PREPROCESSING.out.proteins
        metadata_ch   = PREPROCESSING.out.metadata
    } else {
        proteins_ch = Channel
            .fromPath(params.input_genomes)
            .splitCsv(header: true)
            .map { row -> tuple(row.id, file(row.proteins_faa)) }
        metadata_ch = Channel.empty()
    }

    // ---- Module 1: dereplication (genome + protein level) ----
    if (params.run_module_1) {
        DEREPLICATION(genomes_ch, proteins_ch)
        rep_proteins_ch = DEREPLICATION.out.representative_proteins
        derep_map_ch    = DEREPLICATION.out.mapping_table
        // Keep original genomes_ch (tuple) for downstream; derep list not needed
        // since trait matrix auto-derives genomes from annotations
    } else {
        // Dereplication skipped: concatenate all proteins into one representative FASTA
        rep_proteins_ch = proteins_ch.map { id, faa -> faa }
            .collectFile(name: 'representative_sequences.faa') { f -> f.readLines().join('\n') }
            .map { file(it) }
        derep_map_ch = Channel.value('')
    }

    // ---- Module 2: Branch A — homology annotation (ALL representatives) ----
    HOMOLOGY_ANNOTATION(rep_proteins_ch)
    homology_results_ch   = HOMOLOGY_ANNOTATION.out.homology_results
    orthogroups_ch        = HOMOLOGY_ANNOTATION.out.orthogroups
    unresolved_ch         = HOMOLOGY_ANNOTATION.out.unresolved_proteins
    trait_candidates_ch   = HOMOLOGY_ANNOTATION.out.trait_candidates

    // ---- Module 3: Branch B — staged structure annotation ----
    if (params.run_module_3) {
        STRUCTURE_ANNOTATION(unresolved_ch, trait_candidates_ch)
        struct_results_ch = STRUCTURE_ANNOTATION.out.structure_results
        triage_stats_ch   = STRUCTURE_ANNOTATION.out.triage_stats
    } else {
        struct_results_ch = Channel.empty()
        triage_stats_ch   = Channel.empty()
    }

    // ---- Module 4: combined trait matrix ----
    if (params.run_module_4) {
        TRAIT_MATRIX(
            orthogroups_ch,
            homology_results_ch,
            struct_results_ch,
            derep_map_ch
        )
        trait_matrix_ch = TRAIT_MATRIX.out.trait_matrix
    }

    // ---- Module 5: NISE / convergence detection ----
    if (params.run_module_5 && trait_matrix_ch) {
        NISE_DETECTION(
            trait_matrix_ch,
            homology_results_ch,
            struct_results_ch,
            orthogroups_ch
        )
        nise_ch = NISE_DETECTION.out.nise_candidates
    }

    // ---- Module 6: functional vs identity clustering ----
    if (params.run_module_6 && trait_matrix_ch) {
        CLUSTERING_COMPARISON(trait_matrix_ch, orthogroups_ch)
    }

    // ---- Module 7: validation ----
    if (params.run_module_7) {
        VALIDATION(
            nise_ch.ifEmpty(Channel.empty()),
            struct_results_ch.ifEmpty(Channel.empty()),
            trait_matrix_ch.ifEmpty(Channel.empty())
        )
    }

    // ---- Triage funnel statistics ----
    TRIAGE_REPORT(
        triage_stats_ch.ifEmpty(Channel.empty()),
        HOMOLOGY_ANNOTATION.out.homology_stats
    )
}

workflow.onComplete {
    log.info """
        =========================================================
        Convergence Genomics Pipeline — COMPLETE
        Output directory: ${params.outdir}
        =========================================================
        """.stripIndent()
}
