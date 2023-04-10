#!/bin/bash nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2

include { split_if_contains } from './modules/functions.nf'
include { STOREFASTACHUNKS } from './modules/utils.nf'


if (params.target_db == "PDB") {
    include { GET_PDB_STRUCTURES } from './workflows/get_structures'
} else if (params.target_db == "UniProtKB") {
    include { GET_UNIPROT_STRUCTURES } from './workflows/get_structures'
}


// Prepare input channels
seqs = Channel.fromPath( "$params.seqs" )
target_db = Channel.fromPath( "${params.dbdir}/${params.target_db}",checkIfExists: true ).map { item -> [ item.baseName, item] }


seqs.view()
log.info """\
         S T R U C T U R E S  F E T C H E R  ~  version 1.0.0    
         ======================================="""
log.info "seqs                    : ${params.seqs}"
log.info "target_db               : ${params.target_db}"
log.info "min_id_mmseqs           : ${params.min_id_mmseqs}"
log.info "min_cov_mmseqs          : ${params.min_cov_mmseqs}"
log.info "covmode_mmseqs          : ${params.covmode_mmseqs}"
log.info "min_id_filter           : ${params.min_id_filter}"
log.info "min_cov_filter          : ${params.min_cov_filter}"
log.info "foldseek_convert        : ${params.foldseek_convert}"


workflow FETCH_STRUCTURES {

    // Prepare input channels
    if(params.splitfa == "true"){
        seqs = seqs.splitFasta( by: params.fasta_chunk_size, file: true ).map { item -> [item.simpleName, item.baseName,item] }
        STOREFASTACHUNKS(seqs)
        fastas = STOREFASTACHUNKS.out.fasta_chunk
    }else{
        fastas = seqs.map { item -> [ item.baseName,item] }
    }

    // Fetch the structures from the target database
    if(params.target_db == "PDB") {
      structures = GET_PDB_STRUCTURES (fastas, target_db).out.structures
    }
    else if (params.target_db == "UniProtKB") {
      structures = GET_UNIPROT_STRUCTURES (fastas, target_db, params.min_id_filter, params.min_cov_filter).out.structures
    }

    emit: 
        structures
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    structures = FETCH_STRUCTURES()

}

/*
 * completion handler
 */
 workflow.onComplete {

     def msg = """\
         Pipeline execution summary
         ---------------------------
         Completed at: ${workflow.complete}
         Duration    : ${workflow.duration}
         Success     : ${workflow.success}
         workDir     : ${workflow.workDir}
         exit status : ${workflow.exitStatus}
         """
         .stripIndent()

}
