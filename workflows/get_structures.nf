#!/bin/bash nextflow
include { split_if_contains } from '../modules/functions.nf'
include { MMSEQS_SEARCH;FILTER_HITS;FETCH_STRUCTURES_AF2DB; }  from '../modules/structures.nf'


workflow GET_UNIPROT_STRUCTURES {

  take: 
    seqs
    target_db
    min_id_filter
    min_cov_filter

    main:
        // 1. Find sequence hits in the target database
        MMSEQS_SEARCH(seqs,target_db.collect())

        // 2. Create the template file and obtain the best mmseqs hit
        mmseqs_hits = FILTER_HITS(MMSEQS_SEARCH.out.hits.filter{ it[2].size()>0 }, min_id_filter, min_cov_filter)

        // 3. Download the structures
        //FETCH_STRUCTURES_AF2DB(mmseqs_hits)

}