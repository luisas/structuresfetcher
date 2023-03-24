#!/bin/bash nextflow
//params.outdir = 'results'

include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()


process MMSEQS_PREP_DB {
    container 'soedinglab/mmseqs2'
    storeDir "${params.dbdir}/dbs/${params.target_db}"
    label 'process_small'
    tag "$db_id"

    input:
    tuple val(db_id), file(db)

    output:
    path db

    script:
    """
    mmseqs createindex $db/$db_id tmp
    """
}


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

    script:
    """
    mmseqs easy-search --min-seq-id ${params.min_id_mmseqs} -c ${params.min_cov_mmseqs} --cov-mode ${params.covmode_mmseqs} ${seqs} ${db}/${db_id} hits.m8 tmp --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov
    """

}


process FILTER_HITS {
    container 'luisas/python:bio3'
    storeDir "${params.outdir}/filtered_ids/${db_id}/$id/"
    label 'process_low'
    tag "$id in $db_id"

    input:
    tuple val(id), val(db_id), file(hits)

    output:
    tuple val(id), val(db_id), file("filtered_hits.m8"), file("template.txt"), file("ids_to_download.txt"), emit: filtered_hits

    script:
    """
    filter_hits.py ${hits} "filtered_hits.m8"
    """
}


process FETCH_STRUCTURES_UNIPROT {
    container 'luisas/python:bio3'
    storeDir "${params.outdir}/structures/fetched/${db_id}/$id/"
    label 'process_small'
    tag "$id in $db_id"

    input:
    tuple val(id), val(db_id), file(hits), file(template), file(ids_to_download)

    output:
    tuple val(id), val(db_id), file(hits), file(template), file("*_ref.pdb"), emit: fetched_structures

    script:
    """
    #bash "${path_templates}/scripts/fetch_${db_id}.sh"

    function validate_url(){
      if [[ `wget -S --spider \$1  2>&1 | grep 'HTTP/1.1 200 OK'` ]]; then echo "true"; else echo "false";  fi
    }

    sort $ids_to_download | uniq > ids_to_download_uniq.txt

    for id in \$(cat ids_to_download_uniq.txt); do url="https://alphafold.ebi.ac.uk/files/AF-\$id-F1-model_v4.pdb"; if `validate_url \$url == "true"`; then wget \$url; else echo "does not exist"; fi ; done


    # Horrible coding - please change if keeping
    # Generate links from name of the protein to the one matched_ref.pdb
    python3 "${path_templates}/scripts/getlinks_uniprot.py" ${template} "make_links_tmp.sh" "pdb"
    [ -f ./make_links_tmp.sh ] && tr ', ' ' ' < make_links_tmp.sh > make_links.sh
    [ -f ./make_links.sh ] && bash ./make_links.sh
    [ -f ./make_links.sh ] && cat ./make_links.sh
    """
}


