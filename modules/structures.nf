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
    storeDir "${params.outdir}/structures/fetched/${db_id}/id_${min_id_filter}_cov_${min_cov_filter}/$id/"
    label 'process_medium'
    tag "$id in $db_id"

    input:
    tuple val(id), val(db_id), file(hits), file(template), file(ids_to_download), val(min_id_filter), val(min_cov_filter)

    output:
    tuple val(id), val(db_id),  val(min_id_filter), val(min_cov_filter), file("*.pdb"), emit: fetched_structures
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

    echo "LUISA" > LUISA.txt

    # ----------------------------------------------------
    # Rechange the protein names according to what appears in the template file
    # ----------------------------------------------------
    getlinks_uniprot.py ${template} "make_links_tmp.sh" "pdb"
    [ -f ./make_links_tmp.sh ] && tr ', ' ' ' < make_links_tmp.sh > make_links.sh
    [ -f ./make_links.sh ] && bash ./make_links.sh

    echo "LUISA" > LUISA2.txt
    
    # ----------------------------------------------------
    # Here cut them according to hits
    # ----------------------------------------------------
    cut_structures.py ${hits}
    [ -f ./cut_structures_tmp.sh ] && tr ', ' ' ' < cut_structures_tmp.sh > cut_structures.sh
    [ -f ./cut_structures.sh ] && bash ./cut_structures.sh
    rm temp.pdb
    """
}





// Currently bad, no caching but just or trying out things TODO --> fix
process FETCH_STRUCTURES_PDB {
    container 'luisas/python:bio3'
    storeDir "${params.outdir}/structures/fetched/${db_id}/id_${min_id_filter}_cov_${min_cov_filter}/$id/"
    label 'process_small'
    tag "$id in $db_id"

    input:
    tuple val(id), val(db_id), file(hits), file(template), file(ids_to_download), val(min_id_filter), val(min_cov_filter)

    output:
    tuple val(id), val(db_id),val(min_id_filter), val(min_cov_filter), file(hits), file(template), file("*_ref.pdb"), emit: fetched_structures

    script:
    """
    for id in \$(cat $ids_to_download); do wget https://files.rcsb.org/download/\$id.pdb; done

    # Horrible coding - please change if keeping
    # Generate links from name of the protein to the one matched_ref.pdb
    getlinks_pdb.py ${template} "make_links_tmp.sh" "pdb"
    [ -f ./make_links_tmp.sh ] && tr ', ' ' ' < make_links_tmp.sh > make_links.sh
    [ -f ./make_links.sh ] && bash ./make_links.sh
    [ -f ./make_links.sh ] && cat ./make_links.sh
    """
}



process FETCH_FASTA{

  container 'luisas/python:bio3'
  storeDir "${params.outdir}/structures/fetched/${db_id}/id_${min_id_filter}_cov_${min_cov_filter}/$id/"
  label 'process_small'
  tag "$id in $db_id"

  input:
  tuple val(id), val(db_id), file(hits), file(template), file(ids_to_download), val(min_id_filter), val(min_cov_filter)

  output:
  tuple val(id), val(db_id),val(min_id_filter), val(min_cov_filter), file(hits), file(template), file("*_ref.fa"), emit: fastas

  script:
  """
  sort $ids_to_download | uniq > ids_to_download_uniq.txt

  for id in \$(cat ids_to_download_uniq.txt); do wget https://www.rcsb.org/fasta/entry/\$id -O \$id".fa"; done

  getlinks_pdb.py ${template} "make_links_tmp.sh" "fa"
  [ -f ./make_links_tmp.sh ] && tr ', ' ' ' < make_links_tmp.sh > make_links.sh
  [ -f ./make_links.sh ] && bash ./make_links.sh
  """

}


process PREP_STRUCTURES {
    container 'luisas/tcoffee_python'
    storeDir "${params.outdir}/structures/fetched_preprocessed/${db_id}/id_${min_id_filter}_cov_${min_cov_filter}/$id/"
    label 'process_low'
    tag "$seq_id in $id"

    input:
    tuple val(id), val(seq_id), val(db_id),val(min_id_filter), val(min_cov_filter), file(hits), file(template), file(pdb), file(fasta_ref)

    output:
    tuple val(id), val(db_id), val(seq_id),val(min_id_filter), val(min_cov_filter), file(hits), file(template), file("${seq_id}_ref_ready.pdb"), file("${seq_id}_correct_shift"), emit: structures
    tuple val(id), file("${seq_id}.fa"), emit: fastas

    script:
    """
    # Extract the specific protein hit from the hits file, where all are stored
    tr "/" "_" < $hits > hitsfile
    awk '/^${seq_id}/' hitsfile > hits.txt

    # Extract chain from hits
    CHAIN=\$(awk '{ print \$2 }' hits.txt |  awk -F '[_]' '{ print \$2 }')

    # Extract chain from PDB 
    t_coffee -other_pg extract_from_pdb -chain \$CHAIN -force -infile $pdb > ${seq_id}_fixed.pdb
    # Get Fasta from PDB
    pdb-seq.py ${seq_id}_fixed.pdb > ${seq_id}.fa

    # Extract the right chain from the fasta
    extract_fasta_chain.py ${fasta_ref} \$CHAIN ${fasta_ref.baseName}_chain.fa

    # Then find-match only finds target sequence if 100% match
    # We want it like this also because its 100% btw the real fasta and the PDB
    SHIFT=\$(find-match.py ${fasta_ref.baseName}_chain.fa ${seq_id}.fa)
    echo \$SHIFT

    # TODO: I need to find a solution for when some mismatch is happening - report

    # Extract the chunk from the template
    awk -v shift=\$SHIFT -v chain=\$CHAIN '{print "extract_structure_uniprot.py", \$9, \$10, \$1"_fixed.pdb", \$1"_ref_ready " chain, shift}' hits.txt > extract_structures.sh
    bash ./extract_structures.sh

    if [[ "\$SHIFT" != "-1" ]]; then
        echo 'passed' > ${seq_id}_correct_shift
    else
        touch ${seq_id}_correct_shift
    fi
    """
}


process ADD_HEADER{
  container 'luisas/structural_regression:17'
  tag "${id}"
  storeDir "${params.outdir}/structures/fetched_preprocessed/${db_id}/id_${min_id_filter}_cov_${min_cov_filter}/$id/"

  label "process_medium"

  input:
  tuple val (id), val(db_id), val(min_id_filter), val(min_cov_filter), path (structures)

  output:
  tuple val(id), val(db_id), val(min_id_filter), val(min_cov_filter), path ("pdbs/*.pdb"), emit: pdb


  script:
  """
  mkdir pdbs

  for seq in ${structures}; do
    seq_id=\$(basename \$seq .pdb)
    t_coffee -other_pg extract_from_pdb -infile \$seq > pdbs/\$seq_id.pdb
  done
  """
}



process RENAME_PDB{
  container 'luisas/structural_regression:17'
  tag "${id}"
  storeDir "${params.outdir}/structures/ready/${db_id}/id_${min_id_filter}_cov_${min_cov_filter}/$id/"

  label "process_low"

  input:
  tuple val(id),val(db_id), val(seq_id), val(min_id_filter), val(min_cov_filter), path (pdb_header)

  output:
	tuple val(id), val(db_id), val(min_id_filter), val(min_cov_filter), path ("${seq_id}.pdb"), emit: pdb


  script:
  """
  # Add the headers
  cp $pdb_header ${seq_id}.pdb
  """
}



process FILTER_STRUCTURES{

  container 'luisas/tcoffee_python'
  tag "${id}"
  storeDir "${params.outdir}/structures/fetched_preprocessed/${db_id}/id_${min_id_filter}_cov_${min_cov_filter}/$id/"

  label "process_medium"

  input:
  tuple val (id), val(db_id), val(min_id_filter), val(min_cov_filter), path (structures)

  output:
  tuple val(id), val(db_id), val(min_id_filter), val(min_cov_filter), path ("pdbs/*.pdb"), emit: pdb


  script:
  """
  # Add the headers
  cp $pdb_header ${seq_id}.pdb
  """


}