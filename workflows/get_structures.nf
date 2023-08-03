#!/bin/bash nextflow
include { split_if_contains } from '../modules/functions.nf'
include { MMSEQS_SEARCH;FILTER_HITS;FETCH_STRUCTURES_AF2DB; ADD_HEADER;ADD_HEADER as ADD_HEADER_2;FILTER_STRUCTURES;  FETCH_STRUCTURES_PDB; FETCH_FASTA; PREP_STRUCTURES; RENAME_PDB}  from '../modules/structures.nf'


workflow GET_UNIPROT_STRUCTURES {

  take: 
    seqs
    target_db
    min_id_filter
    min_cov_filter

    main:

        // 1. Find sequence hits in the target database
        MMSEQS_SEARCH(seqs,target_db.collect())

        MMSEQS_SEARCH.out.hits.view()
        // 2. Create the template file and obtain the best mmseqs hit
        FILTER_HITS(MMSEQS_SEARCH.out.hits.filter{ it[2].size()>0 }, min_id_filter, min_cov_filter)

        // 3. Download the structures (and cut them according to the positions the hits)
        FETCH_STRUCTURES_AF2DB(FILTER_HITS.out.filtered_hits)

        // Check that the sequence of the fetched structures is not shifted
        //FILTER_STRUCTURES(FETCH_STRUCTURES_AF2DB.out.fetched_structures)

        // 4. Add the header to the structures
        ADD_HEADER_2(FETCH_STRUCTURES_AF2DB.out.fetched_structures)
        
   emit: 
        structures = ADD_HEADER_2.out.pdb
       
}


workflow GET_PDB_STRUCTURES {

  take: 
    seqs
    target_db
    min_id_filter
    min_cov_filter

  main: 
    // 1. Find sequence hits in PDB
    MMSEQS_SEARCH(seqs,target_db.collect())

    // 2. Create the template file and obtain the best mmseqs hit
    // The mmseqs hit comes with start and end positions of the match, which I will later cut.
    //mmseqs_hits = TEMPLATE_FROM_DB_HITS(MMSEQS_SEARCH.out.hits.filter{ it[2].size()>0 })
    FILTER_HITS(MMSEQS_SEARCH.out.hits.filter{ it[2].size()>0 }, min_id_filter, min_cov_filter)

    // 3. Download the structures
    //    and cut them according to the positions the hits (extract only the real matching chunk)
    FETCH_STRUCTURES_PDB(FILTER_HITS.out.filtered_hits)
    fetched_structures = FETCH_STRUCTURES_PDB.out.fetched_structures.map{ it -> [split_if_contains(it[0], "-ref", 0) , it[1], it[2], it[3], it[4], it[5], it[6]] }
                                            .transpose()
                                            .map{ it -> [ it[0], split_if_contains(it[6].baseName, "_ref", 0), it[1], it[2], it[3], it[4], it[5], it[6] ]}
    // 3b. Download the fastas
    FETCH_FASTA(FILTER_HITS.out.filtered_hits)
    fetched_fastas = FETCH_FASTA.out.fastas.map{ it -> [split_if_contains(it[0], "-ref", 0) , it[1], it[2], it[3], it[4], it[5], it[6]] }
                                            .transpose()
                                            .map{ it -> [ it[0], split_if_contains(it[6].baseName, "_ref", 0), it[1], it[2], it[3], it[4], it[5], it[6] ]}


    structures = fetched_structures.join(fetched_fastas, by: [0,1,2,3,4]).map{ it -> [ it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7], it[10]]}
    // 4. Prep PDB (extract_from_pdb)
    PREP_STRUCTURES(structures)
    correct_structures =  PREP_STRUCTURES.out.structures.filter{ it[8].size()>0 }
    correct_structures.view()
    preprocessed_structures = correct_structures.map{ it -> [ it[0], it[1],it[3],it[4], it[7]]}.groupTuple(by :[0,1,2,3])
    
    reference_structures = ADD_HEADER(preprocessed_structures).transpose().unique()
    reference_structures = reference_structures.map{ it -> [it[0], it[1], split_if_contains(it[4].baseName, "_ref", 0), it[2],  it[3], it[4]]}
    RENAME_PDB(reference_structures)


            
   emit: 
        structures = RENAME_PDB.out.pdb

}
