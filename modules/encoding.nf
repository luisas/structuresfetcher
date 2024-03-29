process STRUCTURE_TO_3DI{
    container 'luisas/foldseek_tcoffee:2'
    tag "$id"
    storeDir "${params.outdir}/foldseek/mapping/${db}/id_${min_id_filter}_cov_${min_cov_filter}/$id/"
    label 'process_low'
    
    input:
    tuple val(id), val(db),val(min_id_filter), val(min_cov_filter), path (structures)
    
    output:
    tuple val(id), val(db),val(min_id_filter), val(min_cov_filter), path ("*3di.out"), emit: mapping
    path ".command.trace", emit: metricFile

    
    script:
    """
    # Convert structures to 3di
    
    for structure in *.pdb; do
        st_id=\$(echo \$structure | cut -d'.' -f1)
        foldseek structureto3didescriptor \$structure \${st_id}_3di
        cut -f1,2,3 \${st_id}_3di > \${st_id}_3di_temp.out
        echo -n `cut -f1 \${st_id}_3di_temp.out | cut -f1 -d' '` > \${st_id}_3di.out
        echo -e -n ' \t ' >> \${st_id}_3di.out
        echo -n `cut -f2 \${st_id}_3di_temp.out` >> \${st_id}_3di.out
        echo -e -n ' \t ' >> \${st_id}_3di.out
        echo `cut -f3 \${st_id}_3di_temp.out` >> \${st_id}_3di.out
    done
    """
}


process MERGE_MAPPINGS {

  container 'luisas/python:bio3'
  storeDir "${params.outdir}/foldseek/mapping_merged/${db}/id_${min_id_filter}_cov_${min_cov_filter}/$id/"
  label 'process_low'

  input:
  tuple val(id), val(db),val(min_id_filter), val(min_cov_filter), file(files)

  output:
  tuple val(id), val(db), val(min_id_filter), val(min_cov_filter),file("${id}.mapping"), emit: mapping

  script:
  """
  for file in $files; do cat \$file >>  "${id}.mapping"; done
  """
}

process PREP_FS_SEQS{
    container 'luisas/foldseek_tcoffee:2'
    tag "$id"
    storeDir "${params.outdir}/foldseek/tcoffee_templates/${db}/id_${min_id_filter}_cov_${min_cov_filter}/$id/"
    label 'process_low'
    
    input:
    tuple val(id), val(db),val(min_id_filter), val(min_cov_filter), path(mapping)
    
    output:
    tuple val (id), val(db),val(min_id_filter), val(min_cov_filter), path("${id}_fs"), emit: foldseek_db
    
    script:
    """
    fs_prep.py $mapping ${id}_fs
    """
}