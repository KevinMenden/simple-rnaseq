#!/usr/bin/env nextflow
/*
========================================================================================
                         simple-rnaseq
========================================================================================

----------------------------------------------------------------------------------------
*/



/*
 * SET UP CONFIGURATION VARIABLES
 */



// Validate inputs
if( params.star_index ){
    star_index = Channel
        .fromPath(params.star_index)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
}
else if ( params.fasta ){
    ch_fasta_for_star_index = Channel.fromPath(params.fasta)
        .ifEmpty { exit 1, "Fasta file not found: ${params.fasta}" }
}
else {
    exit 1, "No reference genome specified!"
}

if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .into { gtf_makeSTARindex; gtf_star}
} else {
    exit 1, "No GTF annotation specified!"
}


/*
 * Create a channel for input read files
 */

Channel
        .fromFilePairs( params.reads, size: 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { read_files_fastqc; read_files_star }





/*
 * STEP 2  - Build STAR index
 */
if(!params.star_index){

    process makeSTARindex {
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from ch_fasta_for_star_index
        file gtf from gtf_makeSTARindex

        output:
        file "star" into star_index

        script:
        """
        mkdir star
        STAR \\
            --runMode genomeGenerate \\
            --runThreadN ${task.cpus} \\
            --sjdbGTFfile $gtf \\
            --genomeDir star/ \\
            --genomeFastaFiles $fasta
        """
    }
}


/**
 * STEP 3 - STAR alignment
 */


process star {
    tag "$name"
    publishDir "${params.outdir}/STAR", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".bam") == -1) "logs/$filename"
                else  filename }

    input:
    set val(name), file(reads) from read_files_star
    file index from star_index.collect()
    file gtf from gtf_star.collect()

    output:
    file '*.bam' into star_aligned
    file "*.out" into alignment_logs
    file "*SJ.out.tab"
    file "*Log.out" into star_log

    script:

    """
    STAR --genomeDir $index \\
        --sjdbGTFfile $gtf \\
        --readFilesIn $reads \\
        --runThreadN ${task.cpus} \\
        --outSAMtype BAM SortedByCoordinate \\
        --readFilesCommand zcat \\
        --runDirPerm All_RWX \\
        --outFileNamePrefix $name \\
        --outFilterMatchNmin ${params.min_aln_length}
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {
 log.info "Pipeline Complete"
}
