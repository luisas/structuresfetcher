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
        FILTER_HITS(MMSEQS_SEARCH.out.hits.filter{ it[2].size()>0 }, min_id_filter, min_cov_filter)

        // 3. Download the structures
        FETCH_STRUCTURES_AF2DB(FILTER_HITS.out.filtered_hits)
}


workflow GET_PDB_STRUCTURES {

  take: 
    seqs
    target_db
    min_id_filter
    min_cov_filter

  main: 
    // 1. Find sequence hits in PDB
    MMSEQS_SEARCH(refs_ch,target_db.collect())

    // 2. Create the template file and obtain the best mmseqs hit
    // The mmseqs hit comes with start and end positions of the match, which I will later cut.
    mmseqs_hits = TEMPLATE_FROM_DB_HITS(MMSEQS_SEARCH.out.hits.filter{ it[2].size()>0 })
    
    // 3. Download the structures
    //    and cut them according to the positions the hits (extract only the real matching chunk)
    FETCH_STRUCTURES(mmseqs_hits)
    fetched_structures = FETCH_STRUCTURES.out.fetched_structures.map{ it -> [split_if_contains(it[0], "-ref", 0) , it[2], it[3], it[4]] }
                                            .transpose()
                                            .map{ it -> [ it[0], split_if_contains(it[3].baseName, "_ref", 0), it[1], it[2], it[3] ]}

    // 3b. Download the fastas
    FETCH_FASTA(mmseqs_hits)
    fetched_fastas = FETCH_FASTA.out.fastas.map{ it -> [split_if_contains(it[0], "-ref", 0) , it[2], it[3], it[4]] }
                                            .transpose()
                                            .map{ it -> [ it[0], split_if_contains(it[3].baseName, "_ref", 0), it[1], it[2], it[3] ]}

    structures = fetched_structures.join(fetched_fastas, by: [0,1]).map{ it -> [ it[0], it[1], it[2], it[3], it[4], it[7]]}

    // 4. Prep PDB (extract_from_pdb)
    PREP_STRUCTURES(structures)
    preprocessed_structures = PREP_STRUCTURES.out.structures.map{ it -> [ it[0], it[1], it[4]]}
    reference_structures = ADD_HEADER(preprocessed_structures)
    RENAME_PDB(reference_structures)

}