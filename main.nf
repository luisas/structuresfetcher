#!/bin/bash nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2

include { split_if_contains } from './modules/functions.nf'
include { STOREFASTACHUNKS } from './modules/utils.nf'
include { FOLDSEEK_CONVERT } from './workflows/foldseek.nf'
include { SAVE_MERGED; SAVE_MERGED_DIR } from './modules/utils.nf'
include { GET_PDB_STRUCTURES; GET_UNIPROT_STRUCTURES } from './workflows/get_structures'



// Prepare input channels
seqs = Channel.fromPath( "${params.seqs}" )
if(params.target_db !in ['AF2_PRED','PDB_fetch_test'] ){
    target_db = Channel.fromPath( "${params.dbdir}/${params.target_db}",checkIfExists: true ).map { item -> [ item.baseName, item] }
}

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


seqs.view()
print(params.seqs)

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
        structures = GET_PDB_STRUCTURES (fastas, target_db, params.min_id_filter, params.min_cov_filter).unique()
        structures_db = structures.unique().transpose().map{  it -> [it[0].split("\\.")[0], it[1], it[2], it[3], it[4] ] }.groupTuple(by: [0,1,2,3])
        SAVE_MERGED(structures_db, "structures")
    }
    else if (params.target_db == "UniProtKB") {
        structures = GET_UNIPROT_STRUCTURES (fastas, target_db, params.min_id_filter, params.min_cov_filter)
        structures_db = structures.transpose().map{  it -> [it[0].split("\\.")[0], it[1], it[2], it[3], it[4] ] }.groupTuple(by: [0,1,2,3])
        structures_db.view()
        SAVE_MERGED(structures_db, "structures")
    }
    else if (params.target_db == "AF2_PRED") {
        structures = Channel.fromPath(params.structures_path)
        structures = structures.map{  it -> [it.getParent().getParent().getBaseName(),"AF2_PRED", it] }.groupTuple(by: [0,1])
                               .map{ it -> [it[0], it[1], "1.0", "1.0", it[2]]}
    }
    if (params.foldseek_convert) {
        structures_db.view()
        foldseek_db = FOLDSEEK_CONVERT(structures_db)
        foldseek_db = foldseek_db.map{  it -> [it[0].split("\\.")[0], it[1], it[2], it[3], it[4] ] }.groupTuple(by: [0,1,2,3])
        SAVE_MERGED_DIR(foldseek_db, "foldseek")
    }

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
