#!/bin/bash nextflow
include { split_if_contains } from '../modules/functions.nf'
include { MMSEQS_SEARCH;FILTER_HITS;FETCH_STRUCTURES_UNIPROT; }  from '../modules/structures.nf'


workflow GET_UNIPROT_STRUCTURES {

  take: 
    seqs
    target_db

    main:
  // 1. Find sequence hits in the target database
    MMSEQS_SEARCH(seqs,target_db.collect())


  // 2. Create the template file and obtain the best mmseqs hit
    mmseqs_hits = FILTER_HITS(MMSEQS_SEARCH.out.hits.filter{ it[2].size()>0 })


//   // 3. Download the structures
//   //    and cut them according to the positions the hits (extract only the real matching chunk)
//   FETCH_STRUCTURES_UNIPROT(mmseqs_hits)
//   fetched_structures = FETCH_STRUCTURES_UNIPROT.out.fetched_structures.map{ it -> [split_if_contains(it[0], "-ref", 0) , it[2], it[3], it[4]] }
//                                          .transpose()
//                                          .map{ it -> [ it[0], split_if_contains(it[3].baseName, "_ref", 0), it[1], it[2], it[3] ]}


}