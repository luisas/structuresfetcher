#!/bin/bash nextflow
//params.outdir = 'results'

include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()


process MMSEQS_SEARCH {
    container 'luisas/mmseqs2test'
    storeDir "${params.outdir}/mmseqs/${db_id}/$id/"
    label 'process_medium_high'
    tag "$id in $db_id"

    input:
    tuple val(id), file(seqs)
    tuple val(db_id), file(db)

    output:
    tuple val(id), val(db_id), path("hits.m8"), emit: hits
    path ".command.trace", emit: metricFile

    script:
    """
    mmseqs easy-search --min-seq-id ${params.min_id_mmseqs} -c ${params.min_cov_mmseqs} --cov-mode ${params.covmode_mmseqs} ${seqs} ${db}/${db_id} hits.m8 tmp --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov
    """

}


process FILTER_HITS {
    container 'luisas/python:bio3'
    storeDir "${params.outdir}/filtered_ids/${db_id}/id_${min_id_filter}_cov_${min_cov_filter}/$id/"
    label 'process_low'
    tag "$id in $db_id"

    input:
    tuple val(id), val(db_id), file(hits)
    val(min_id_filter)
    val(min_cov_filter)

    output:
    tuple val(id), val(db_id), file("${id}_filtered_hits.m8"), file("${id}_template.txt"), file("${id}_ids_to_download.txt"), val(min_id_filter), val(min_cov_filter), emit: filtered_hits
    path ".command.trace", emit: metricFile


    script:
    """
    filter_hits.py ${hits} "${id}_filtered_hits.m8" "${id}_template.txt" "${id}_ids_to_download.txt" "${id}_chains.txt" ${min_id_filter} ${min_cov_filter} 
    """
}


process FETCH_STRUCTURES_AF2DB {
    container 'luisas/python:bio3'
    storeDir "${params.outdir}/structures/fetched/${db_id}/$id/"
    label 'process_low'
    tag "$id in $db_id"

    input:
    tuple val(id), val(db_id), file(hits), file(template), file(ids_to_download), val(min_id_filter), val(min_cov_filter)

    output:
    tuple val(id), val(db_id), file("*.pdb"), emit: fetched_structures
    path ".command.trace", emit: metricFile

    script:
    """

    # ----------------------------------------------------
    # Download the structures from the AlphaFold2 database
    # ----------------------------------------------------

    function validate_url(){
      if [[ `wget -S --spider \$1  2>&1 | grep 'HTTP/1.1 200 OK'` ]]; then echo "true"; else echo "false";  fi
    }

    for id in \$(cat $ids_to_download); do url="https://alphafold.ebi.ac.uk/files/AF-\$id-F1-model_v4.pdb"; if `validate_url \$url == "true"`; then wget \$url; else echo "does not exist"; fi ; done


    # ----------------------------------------------------
    # Rechange the protein names according to what appears in the template file
    # ----------------------------------------------------
    getlinks_uniprot.py ${template} "make_links_tmp.sh" "pdb"
    [ -f ./make_links_tmp.sh ] && tr ', ' ' ' < make_links_tmp.sh > make_links.sh
    [ -f ./make_links.sh ] && bash ./make_links.sh
    rm ./AF-*

    # ----------------------------------------------------
    # Here cut them according to hits
    # ----------------------------------------------------
    cut_structures.py ${hits}
    [ -f ./cut_structures_tmp.sh ] && tr ', ' ' ' < cut_structures_tmp.sh > cut_structures.sh
    [ -f ./cut_structures.sh ] && bash ./cut_structures.sh
    rm temp.pdb
    """
}





