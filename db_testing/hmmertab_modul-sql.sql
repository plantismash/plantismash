SELECT  
                    gbacc_protid_mapping.protein_id, 
                    hmmer_results.env_from,
                    hmmer_results.env_to, 
                    hmmer_results.domain_score, 
                    hmmer_results.domain_evalue
                 FROM 
                    asmash.hmmer_results, 
                    asmash.gbacc_protid_mapping
                 WHERE 
                    gbacc_protid_mapping.protein_id = hmmer_results.acc AND
                    gbacc_protid_mapping.genbank_id = 'Y16952.3' AND
                    hmmer_results.query_name = 'Condensation'
                 ORDER BY gbacc_protid_mapping.protein_id ASC, hmmer_results.domain_num ASC;