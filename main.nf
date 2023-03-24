#!/bin/bash nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2



if (params.target_db == "PDB") {
    include { GET_PDB_STRUCTURES } from './workflows/get_structures'
} else if (params.target_db == "UniProtKB") {
    include { GET_UNIPROT_STRUCTURES } from './workflows/get_structures'
}


// Prepare input channels
seqs = Channel.fromPath( "$params.seqs" ).map { item -> [ item.baseName,item] }
target_db = Channel.fromPath( "${params.dbdir}/${params.target_db}",checkIfExists: true ).map { item -> [ item.baseName, item] }


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



workflow FETCH_STRUCTURES {

    if(params.target_db == "PDB") {
      GET_PDB_STRUCTURES (seqs, target_db)
    }
    else if (params.target_db == "UniProtKB") {
      GET_UNIPROT_STRUCTURES (seqs, target_db, params.min_id_filter, params.min_cov_filter)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    FETCH_STRUCTURES ()
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
