

process STOREFASTACHUNKS{

    storeDir "${params.outdir}/fasta_chunks/$id/$chunk/"
    label "process_low"

    input:
    tuple val(id),val(chunk), file(fasta)

    output:
    tuple val(chunk), file("fasta/*.fa"), emit:  fasta_chunk
    
    script:
    """
    echo "Just storing the file"
    mkdir -p fasta
    cp $fasta fasta/
    """
}

process SAVE_MERGED{
    
    storeDir "${params.outdir}/DB/${folder_name}/${db}/id_${min_id_filter}_cov_${min_cov_filter}/"
    label "process_low"

    input:
    tuple val(id), val(db),val(min_id_filter), val(min_cov_filter),file(files)
    val(folder_name)

    output:
    tuple val(id), val(db),val(min_id_filter), val(min_cov_filter), file("$id"), emit: all_files

    script:
    """
    mkdir -p $id
    for file in $files; do cp \$file $id/ ; done
    """
}

process SAVE_MERGED_DIR{
    
    storeDir "${params.outdir}/DB/${folder_name}/${db}/id_${min_id_filter}_cov_${min_cov_filter}/"
    label "process_low"

    input:
    tuple val(id), val(db),val(min_id_filter), val(min_cov_filter), file(dirs)
    val(folder_name)

    output:
    tuple val(id), val(db) ,val(min_id_filter), val(min_cov_filter), file("$id"), emit: all_files

    script:
    """
    mkdir -p $id
    for dir in $dirs; do cp \$dir/* $id/ ; done
    """
}