

process STOREFASTACHUNKS{

    storeDir "${params.outdir}/fasta_chunks/$id/$chunk/"

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