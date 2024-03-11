-- mutclust_dat.sql  
-- 09/10/2019      

-----------------------------------------------------------------------------------------------------------
-- Oracle SQL views:
--  management (mgmt) views:    mutclust/sql/mutclust_mgmt.sql
--  data views:                 mutclust/sql/mutclust_dat.sql

-- Runs (i.e. selecting correct mgmt/data for specific runs):
-- Moved to run directory e.g.  ~/Dropbox/bioinf/neofun/run/mutclust/dt12a/mc03/sql/mc03.sql

-- Note: moved from sql_41.sql on 09/10/19
-----------------------------------------------------------------------------------------------------------
-- [1] MUTCLUST_DAT (missense and silent) - FULL COUNTS i.e. hotspot residues included
        -- local:       ~/Dropbox/bioinf/neofun/run/mutclust/dt12a/mutdat
        -- pchuckle:    /home/ucbtshf/run/neofun/dt12a/mutclust/mutdat

        --  MFC format
        -- SUPERFAMILY_ID	FUNFAM_NUMBER	EF	GENE	FUN_UP	REP_ID	PDB_CODE	CHAIN_CODE	AA_POS	PDB_RES_NAME	MFC_MUT_COUNT
        -- 3.30.70.330	43708	99	gene_na	up_na	2dguA00	2dgu	A	-1	354	3
        -- 3.30.70.330	43381	99	gene_na	up_na	2disA01	2dis	A	-1	76	3
        -- 3.30.70.330	43708	99	gene_na	up_na	2dguA00	2dgu	A	-1	397	2
        -- 3.30.70.330	43381	99	gene_na	up_na	2disA01	2dis	A	-1	66	1

        -- 14/10/19 - Extended to include number of SwissProt relatives in FunFam 
        -- SUPERFAMILY_ID	FUNFAM_NUMBER	EF	GENE	FUN_UP	REP_ID	PDB_CODE	CHAIN_CODE	AA_POS	PDB_RES_NAME    NUM_SWISSPROT_IN_FF MFC_MUT_COUNT   
        
        -- (a) MutDat (all)
            DROP VIEW   funvar_tx.VW_MUTCLUST_DAT;
            CREATE VIEW funvar_tx.VW_MUTCLUST_DAT
            AS
                SELECT
                    data_source,
                    vm_synonymous,
                    variant_class,
                    sf_id SUPERFAMILY_ID,
                    ff_id FUNFAM_NUMBER,
                    '99' EF,
                    vm_gene GENE,
                    vm_uniprot_accession FUN_UP,
                    rep_id REP_ID,
                    SUBSTR( rep_id, 1, 4 ) PDB_CODE,
                    SUBSTR( rep_id, 5, 1 ) CHAIN_CODE,
                    vm_seq_no AA_POS,
                    mm_output_pdb_res PDB_RES_NAME,
                    num_ff_members NUM_SWISSPROT_IN_FF,
                    COUNT( * ) MFC_MUT_COUNT_RESIDUE

                FROM
                    funvar_tx.mvw_map_core_funfam_pdb  mcc

                LEFT JOIN
                    funvar_tx.vw_funfam_swiss_member_count smc
                    ON 
                        mcc.sf_id = smc.superfamily_id 
                    AND
                        mcc.ff_id = smc.funfam_number 
                
                WHERE 
                        mm_output_pdb_res <> '-' 

                GROUP BY 
                    data_source, 
                    vm_synonymous, 
                    variant_class, 
                    sf_id, 
                    ff_id, 
                    '99', 
                    vm_gene, 
                    vm_uniprot_accession, 
                    rep_id, 
                    SUBSTR( rep_id, 1, 4 ), 
                    SUBSTR( rep_id, 5, 1 ), 
                    vm_seq_no, 
                    mm_output_pdb_res,
                    num_ff_members
                
                ORDER BY
                    sf_id, ff_id, vm_gene, vm_seq_no;
                    
        -- (b) Missense
           SELECT
                superfamily_id,
                funfam_number,
                ef,
                gene,
                fun_up,
                rep_id,
                pdb_code,
                chain_code,
                aa_pos,
                pdb_res_name,
                num_swissprot_in_ff,
                mfc_mut_count_residue
            FROM
                funvar_tx.vw_mutclust_dat
                
            WHERE 
                data_source = 'TCGA'
                AND
                vm_synonymous = 'FALSE'
                AND
                variant_class = 'Missense_Mutation';

        -- (c) Silent
            SELECT
                superfamily_id,
                funfam_number,
                ef,
                gene,
                fun_up,
                rep_id,
                pdb_code,
                chain_code,
                aa_pos,
                pdb_res_name,
                num_swissprot_in_ff,
                mfc_mut_count_residue
            FROM
                funvar_tx.vw_mutclust_dat
                
            WHERE 
                data_source = 'TCGA'
                AND
                vm_synonymous = 'TRUE'
                AND
                variant_class = 'Silent';

