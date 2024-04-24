#!/usr/bin/env nextflow

// tools: skesa for assembly, mlst for genotyping, quast for quality assessment
// mamba: create -n nextflow -c bioconda -c conda-forge 
// mamba: skesa mlst quast nextflow

params.reads = "${PWD}/reads"


process run_SKESA {
    publishDir "${PWD}/asm", mode: 'copy'
    
    input:
    tuple val(sample_id), path(read_1), path(read_2)
    
    output:
    tuple val(sample_id), path("${sample_id}.fa")
    
    script:
    """
    mkdir -p ${PWD}/asm
    skesa --reads ${read_1},${read_2} > ${sample_id}.fa
    """
}

process run_QUAST {
    publishDir "${PWD}/quast", mode: 'copy'
    
    input:
    tuple val(sample_id), path(assembly)
    
    output:
    path 'quast_results'
    
    script:
    """
    mkdir -p ${PWD}/quast
    quast.py ${assembly}
    """
}

process run_MLST {
    publishDir "${PWD}/mlst", mode: 'copy'
    
    input:
    tuple val(sample_id), path(assembly)
    
    output:
    path("${sample_id}_MLST_Summary.tsv")
    
    script:
    """
    mkdir -p ${PWD}/mlst
    mlst ${assembly} > ${sample_id}_MLST_Summary.tsv
    """
}

workflow {
    read_pairs_channel = channel.fromFilePairs("${params.reads}/E*_{R1,R2}_phix.fq.gz", flat: true, checkIfExists: true)
    
    run_SKESA(read_pairs_channel)
    run_QUAST(ASSEMBLE.out)
    run_MLST(ASSEMBLE.out)
}