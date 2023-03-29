


workflow FOLDSEEK_CONVERT{
    take:
        structures

    main:
        echo "FOLDSEEK_CONVERT"
        echo "structures: ${structures}"
        echo "params.foldseek_convert: ${params.foldseek_convert}"
        echo "params.outdir: ${params.outdir}"
        echo "params.outdir_foldseek: ${params.outdir_foldseek}"
    
    emit: 
        foldseek_db

}