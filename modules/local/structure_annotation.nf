// =============================================================================
// Module 3 — Branch B: Structure-based, homology-independent annotation
// Staged triage: B1 (cheap ML) → B2 (ProstT5+3Di) → B3 (full 3D fold)
// =============================================================================

process B1_SEQUENCE_ML {
    tag        "b1_ml"
    label      'gpu_cheap'
    publishDir "${params.outdir}/03_structure/b1_sequence_ml", mode: 'copy'

    input:
    path unresolved_faa
    path trait_faa

    output:
    path "b1_predictions.tsv",       emit: predictions
    path "b1_resolved.faa",          emit: resolved
    path "b1_unresolved.faa",        emit: still_unresolved

    when:
    params.run_module_3

    script:
    def conf = params.b1_confidence
    """
    # Merge unresolved set with trait candidates (cross-validation set)
    cat ${unresolved_faa} ${trait_faa} 2>/dev/null | seqkit rmdup -n -o b1_input.faa -D b1_duplicates.txt

    # Run ProteInfer (sequence-only EC/GO prediction) — GPU but cheap
    if command -v proteinfer >/dev/null 2>&1; then
        proteinfer predict \
            --output b1_proteinfer.tsv \
            --batch-size 64 \
            b1_input.faa
    fi

    # Run CLEAN (contrastive EC predictor)
    if command -v CLEAN >/dev/null 2>&1; then
        python -m CLEAN.predict \
            --fasta b1_input.faa \
            --output b1_clean.tsv \
            --batch_size 64 \
            -T \${CLEAN_TABLE:-/opt/CLEAN/enzyme.tsv}
    fi

    # Merge predictions, apply confidence threshold, split resolved/unresolved
    python ${projectDir}/bin/b1_merge.py \
        --proteinfer b1_proteinfer.tsv \
        --clean b1_clean.tsv \
        --input-fasta b1_input.faa \
        --confidence ${conf} \
        --out-predictions b1_predictions.tsv \
        --out-resolved b1_resolved.faa \
        --out-unresolved b1_unresolved.faa
    """
}

process B2_PROSTT5_FOLDSEEK {
    tag        "b2_fold"
    label      'gpu_cheap'
    publishDir "${params.outdir}/03_structure/b2_prostt5_foldseek", mode: 'copy'

    input:
    path b1_unresolved_faa

    output:
    path "b2_predictions.tsv",           emit: predictions
    path "b2_resolved.faa",              emit: resolved
    path "b2_unresolved.faa",            emit: still_unresolved

    when:
    params.run_module_3

    script:
    def tm  = params.b2_tmscore
    def fid = params.b2_fident
    """
    # Predict 3Di structural tokens with ProstT5 (~1000x faster than 3D fold)
    python ${projectDir}/bin/prostt5_predict.py \
        --fasta ${b1_unresolved_faa} \
        --output b2_3di.fasta \
        --device cuda

    # Foldseek search of 3Di tokens against structural databases
    foldseek easy-search \
        b2_3di.fasta \
        \${FOLDSEEK_DB:-/db/foldseek/v4} \
        b2_foldseek.tsv \
        b2_foldseek_tmp \
        --format-output "query,target,evalue,bits,fident,alnlen,qcov,tcov,tmid" \
        --threads ${task.cpus}

    # Parse hits → fold/EC annotation + split resolved vs residual
    python ${projectDir}/bin/b2_parse_foldseek.py \
        --foldseek b2_foldseek.tsv \
        --input-fasta ${b1_unresolved_faa} \
        --tm-score-cutoff ${tm} \
        --fident-cutoff ${fid} \
        --out-predictions b2_predictions.tsv \
        --out-resolved b2_resolved.faa \
        --out-unresolved b2_unresolved.faa
    """
}

process B3_3D_FOLDING {
    tag        "b3_fold"
    label      'gpu_expensive'
    publishDir "${params.outdir}/03_structure/b3_3d_folding", mode: 'copy'

    input:
    path b2_unresolved_faa

    output:
    path "b3_predictions.tsv",           emit: predictions
    path "b3_structures/**",             emit: structures, optional: true

    when:
    params.run_module_3

    script:
    """
    mkdir -p b3_structures

    # Check AlphaFold DB first to skip precomputed entries
    python ${projectDir}/bin/b3_fold_and_predict.py \
        --fasta ${b2_unresolved_faa} \
        --out-dir b3_structures \
        --predictions b3_predictions.tsv \
        --alphafold-db \${ALPHAFOLD_DB:-/db/alphafold} \
        --foldseek-db \${FOLDSEEK_DB:-/db/foldseek/v4} \
        --use-esm first_pass \
        --deepfri-model \${DEEPFRI_MODEL:-/opt/deepfri} \
        --device cuda
    """
}

process STRUCTURE_AGGREGATE {
    tag        "struct_agg"
    label      'cpu_low'
    publishDir "${params.outdir}/03_structure", mode: 'copy'

    input:
    path b1_predictions
    path b2_predictions
    path b3_predictions

    output:
    path "structure_results.tsv",    emit: structure_results
    path "triage_stats.json",        emit: triage_stats

    when:
    params.run_module_3

    script:
    """
    python ${projectDir}/bin/structure_aggregate.py \
        --b1 ${b1_predictions} \
        --b2 ${b2_predictions} \
        --b3 ${b3_predictions} \
        --out structure_results.tsv \
        --triage-stats triage_stats.json
    """
}

workflow STRUCTURE_ANNOTATION {
    take:
    unresolved_ch, trait_candidates_ch

    main:
    if (params.run_module_3) {
        // unresolved_ch and trait_candidates_ch are single-path channels
        B1 = B1_SEQUENCE_ML(unresolved_ch, trait_candidates_ch)
        B2 = B2_PROSTT5_FOLDSEEK(B1.still_unresolved)
        B3 = B3_3D_FOLDING(B2.still_unresolved)

        STRUCTURE_AGGREGATE(
            B1.predictions,
            B2.predictions,
            B3.predictions
        )
    } else {
        B1 = B1_SEQUENCE_ML(Channel.empty(), Channel.empty())
        B2 = B2_PROSTT5_FOLDSEEK(Channel.empty())
        B3 = B3_3D_FOLDING(Channel.empty())
        STRUCTURE_AGGREGATE(Channel.empty(), Channel.empty(), Channel.empty())
    }

    emit:
    structure_results = STRUCTURE_AGGREGATE.out.structure_results
    triage_stats      = STRUCTURE_AGGREGATE.out.triage_stats
}
