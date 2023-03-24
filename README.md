# S T R U C T E R E S  _  F E T C H E R   ~ v1.0.0

```
nextflow run main.nf -profile homfam_uniprot
```

default filtering parameters: 
  - min_id_mmseqs = 0.7
  - min_cov_mmseqs = 0.7
  - covmode_mmseqs = 2 ( coverage is computed on the query sequence ) 
  - min_id_filter = 0.9
  - min_cov_filter= 1.0
  
The mmseqs filters (min_id_mmseqs, min_cov_mmseqs, covmode_mmseqs) define the parameters with which the search is performed. 
The filter parameters define the the filters applied to the mmseqs results and the results of this filter define which structures will be eventually downloaded. The filtering parameters cannot be more permissive than the mmseqs ones. 
The reason behind this double filtering is that we may want to run the search once with more permissive filters first and experiment later with different filterings without the need of recomputing the search.


TODO not to forget 
- clean containers 
- remove tower token from config
- add prep db module
