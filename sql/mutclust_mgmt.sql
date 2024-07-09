-- mutclust_mgmt.sql  
-- 09/10/2019      
-- 12 06 2024: include AlphaFold version in §2d

-----------------------------------------------------------------------------------------------------------
-- Oracle SQL views (for running MutClust)
--  management (mgmt) views:    mutclust/sql/mutclust_mgmt.sql
--  data views:                 mutclust/sql/mutclust_dat.sql

-- Runs (i.e. selecting correct mgmt/data for specific runs):
-- Moved to run directory e.g.  ~/Dropbox/bioinf/neofun/run/mutclust/dt12a/mc03/sql/mc03.sql

-- Note: moved from sql_41.sql on 09/10/19
-----------------------------------------------------------------------------------------------------------

    -----------------------------------------------------------------------------------------------------------
    -- [1]  Unique position counts view
    --      This view is for flattening all FunFam residue hotspots to 1 mutation - i.e.
    --      only intersted in where mutations can happen in structure
    --      [Note Oct'19: This was planned to be used in MutClust, but instead the R code has been 
    --      changed to clarify clusters based on residue counts or mutation counts (on those residues) 
    -----------------------------------------------------------------------------------------------------------
    DROP VIEW funvar_tx.VW_MUTCLUST_UNQ_MUT_POS;
    CREATE VIEW funvar_tx.VW_MUTCLUST_UNQ_MUT_POS
        AS
        SELECT
            DISTINCT
            sf_id,
            ff_id,
            99 ef,
            '<gene>' gene,
            '<uniprot>' fun_up,
            rep_id,
            SUBSTR(rep_id, 1,4) pdb_code,
            SUBSTR(rep_id, 5,1) chain_code,
            mm_output_pdb_res,
            vm_synonymous,
            variant_class

        FROM
            funvar_tx.mvw_map_core_funfam_pdb 

        WHERE 
            mm_output_pdb_res <> '-' 
        
        AND 
            data_source='TCGA'
            
        ORDER BY sf_id, ff_id, rep_id;

    -----------------------------------------------------------------------------------------------------------
    -- [2]  MUTCLUST_MGMT (missense, silent and the combined (for 'matched' runs))
    --      NOTE some VM_SYNONYMOUS are var_class missense!  
    --      * usse VM definition - i.e. protein-level call *
            -- SELECT vm_synonymous, variant_class, count(*)
            -- FROM funvar_tx.mvw_map_core_funfam_pdb
            -- GROUP BY vm_synonymous, variant_class;
            -- TRUE	Missense_Mutation	75
            -- TRUE	Silent	99095
            -- FALSE	Missense_Mutation	245098
            -- FALSE	Silent	30
            
            -- NOTE: AlphaFold version section (d)
    -----------------------------------------------------------------------------------------------------------
        -- (a) Create VIEW vw_mutclust_mgmt_missense
            DROP VIEW funvar_tx.vw_mutclust_mgmt_missense;
        
            CREATE VIEW funvar_tx.vw_mutclust_mgmt_missense 
            AS
            SELECT 
                ROW_NUMBER( ) OVER( ORDER BY MFC_MUT_COUNT_MISSENSE DESC ) TASKID,
                superfamily_id,
                funfam_number,
                ef,
                gene,
                fun_up,
                rep_id,
                pdb_code,
                chain_code,
                num_swissprot_in_ff,
                MFC_MUT_COUNT_MISSENSE
            FROM
                (   SELECT
                        superfamily_id,
                        funfam_number,
                        99 ef,
                        '<gene>' gene,
                        '<uniprot>' fun_up,
                        rep_id,
                        pdb_code,
                        chain_code,
                        num_swissprot_in_ff,
                        SUM( mfc_mut_count_residue ) MFC_MUT_COUNT_MISSENSE                
                    FROM
                        funvar_tx.vw_mutclust_dat
                    WHERE 
                        data_source = 'TCGA'
                        AND
                        vm_synonymous = 'FALSE'
                        AND
                        variant_class = 'Missense_Mutation'
                    GROUP BY
                        superfamily_id,
                        funfam_number,
                        99,
                        '<gene>',
                        '<uniprot>',
                        rep_id,
                        pdb_code,
                        chain_code,
                        num_swissprot_in_ff
                    HAVING COUNT( * ) > 1      
                    ORDER BY COUNT( * ) DESC
                );
            
        -- (b) Create VIEW vw_mutclust_mgmt_silent
            DROP VIEW funvar_tx.vw_mutclust_mgmt_silent;
        
            CREATE VIEW funvar_tx.vw_mutclust_mgmt_silent 
            AS
            SELECT 
                ROW_NUMBER( ) OVER( ORDER BY MFC_MUT_COUNT_SILENT DESC ) TASKID,
                superfamily_id,
                funfam_number,
                ef,
                gene,
                fun_up,
                rep_id,
                pdb_code,
                chain_code,
                num_swissprot_in_ff,
                MFC_MUT_COUNT_SILENT
            FROM
                (   SELECT
                        superfamily_id,
                        funfam_number,
                        99 ef,
                        '<gene>' gene,
                        '<uniprot>' fun_up,
                        rep_id,
                        pdb_code,
                        chain_code,
                        num_swissprot_in_ff,
                        SUM( mfc_mut_count_residue ) MFC_MUT_COUNT_SILENT                
                    FROM
                        funvar_tx.vw_mutclust_dat
                    WHERE 
                        data_source = 'TCGA'
                        AND
                        vm_synonymous = 'TRUE'
                        AND 
                        variant_class = 'Silent'
                    GROUP BY
                        superfamily_id,
                        funfam_number,
                        99,
                        '<gene>',
                        '<uniprot>',
                        rep_id,
                        pdb_code,
                        chain_code,
                        num_swissprot_in_ff
                    HAVING COUNT( * ) > 1        
                    ORDER BY COUNT( * ) DESC
                );


        -- (c) combine into one (so TASK_IDs are in *same order* for silent and missense)
        DROP VIEW funvar_tx.vw_mutclust_mgmt_matched;
        
        CREATE VIEW funvar_tx.vw_mutclust_mgmt_matched
        AS 

        SELECT
            mis.taskid,
            mis.superfamily_id,
            mis.funfam_number,
            mis.ef,
            mis.gene,
            mis.fun_up,
            mis.rep_id,
            mis.pdb_code,
            mis.chain_code,
            mis.num_swissprot_in_ff,
            mis.mfc_mut_count_missense,
            silent.taskid silent_taskid,
            silent.superfamily_id silent_superfamily_id,
            silent.funfam_number silent_funfam_number,
            silent.rep_id silent_rep_id,
            silent.pdb_code silent_pdb_code,
            silent.chain_code silent_chain_code,
            silent.num_swissprot_in_ff silent_num_swissprot_in_ff,
            silent.mfc_mut_count_silent
        FROM
            funvar_tx.vw_mutclust_mgmt_missense mis

        LEFT JOIN
            funvar_tx.vw_mutclust_mgmt_silent silent 
            ON 
            mis.superfamily_id = silent.superfamily_id
            AND 
            mis.funfam_number =  silent.funfam_number
            AND 
            mis.rep_id = silent.rep_id;

        -- (d) AlphaFold version of (a)
        -- 12 06 2024
            CREATE OR REPLACE VIEW funvar_tx.vw_mutclust_mgmt_af_pre 
            AS
              SELECT 
                ROW_NUMBER( ) OVER( ORDER BY MFC_MUT_COUNT_MISSENSE DESC ) TASKID_rownum,
				-1 TASKID,
                superfamily_id,
                funfam_number,
                ef,
                gene,
                fun_up,
                rep_id,
                cath_rep_or_uniprot uniprot_acc,
				-- SUBSTR(rep_id, 1, INSTR( rep_id, '_' ) -1 ) up,
				REPLACE( SUBSTR(rep_id, INSTR( rep_id, '_' ) +1 ), '_', '-') up_aa_range ,
                -- make rep_id the "PDB_CODE" for AF data
				-- pdb_code,
				rep_id pdb_code,
                chain_code,
                num_swissprot_in_ff,
                MFC_MUT_COUNT_MISSENSE MFC_MUT_COUNT,
				-- MutClust cluster job params
				5 		MIN_RAD,
				5 		MAX_RAD,
				10000	NO_TRIALS
            FROM
                (   SELECT
                        superfamily_id,
                        funfam_number,
                        99 ef,
                        '<gene>' gene,
                        '<uniprot>' fun_up,
                        rep_id,
                        cath_rep_or_uniprot,
                        pdb_code,
                        chain_code,
                        num_swissprot_in_ff,
                        SUM( mfc_mut_count_residue ) MFC_MUT_COUNT_MISSENSE                
                    FROM
                        funvar_tx.vw_mutclust_dat_af
                    WHERE 
                        data_source = 'TCGA'
                        AND vm_synonymous = 'FALSE'
                        AND variant_class = 'Missense_Mutation'
						AND rep_source_id = 'uniprot'
                    GROUP BY
                        superfamily_id,
                        funfam_number,
                        99,
                        '<gene>',
                        '<uniprot>',
                        rep_id,
                        cath_rep_or_uniprot,
                        pdb_code,
                        chain_code,
                        num_swissprot_in_ff
                    HAVING COUNT( * ) > 1      
                    ORDER BY COUNT( * ) DESC
                );

    -------------------------------------------------------------------------
    -- [3] MUTCLUST MGMT files (missense and silent) 
    -------------------------------------------------------------------------
        -- MISSENSE (for 10k and also 1k and 100 test versions)
        -- local:       ~/Dropbox/bioinf/neofun/run/mutclust/dt12a/mutdat/mfc_mgmt_swisscount_missense_5_10000_20191014.tsv
        -- pchuckle:    /home/ucbtshf/run/neofun/dt12a/mutclust/mutdat/mfc_mgmt_swisscount_missense_5_10000_20191014.tsv 
            SELECT
                TASKID,
                SUPERFAMILY_ID,
                FUNFAM_NUMBER,
                EF,
                GENE,
                FUN_UP,
                REP_ID,
                PDB_CODE,
                CHAIN_CODE,
                num_swissprot_in_ff,
                mfc_mut_count_missense MFC_MUT_COUNT,
                5 MIN_RAD,
                5 MAX_RAD,
                10000 NO_TRIALS
            FROM
                funvar_tx.vw_mutclust_mgmt_matched;
            
        -- SILENT (for 10k and also 1k and 100 test versions)
        -- local:       ~/Dropbox/bioinf/neofun/run/mutclust/dt12a/mutdat/mfc_mgmt_swisscount_silent_5_10000_20191014.tsv
        -- pchuckle:    /home/ucbtshf/run/neofun/dt12a/mutclust/mutdat/mfc_mgmt_swisscount_silent_5_10000_20191014.tsv
            SELECT
                TASKID,
                SUPERFAMILY_ID,
                FUNFAM_NUMBER,
                EF,
                GENE,
                FUN_UP,
                REP_ID,
                PDB_CODE,
                CHAIN_CODE,
                num_swissprot_in_ff,
                mfc_mut_count_silent MFC_MUT_COUNT,
                5 MIN_RAD,
                5 MAX_RAD,
                10000 NO_TRIALS
            FROM
                funvar_tx.vw_mutclust_mgmt_matched;

    -------------------------------------------------------------------------
    -- [4] MUTCLUST MGMT imports 
    -------------------------------------------------------------------------
    -- 27/10/19
        UPDATE FUNVAR_ADMIN.mutclust_mgmt
        SET group_id  ='mc05'
        --select * from FUNVAR_ADMIN.mutclust_mgmt
        WHERE group_id ='dummy';

        UPDATE FUNVAR_ADMIN.mutclust_mgmt
        SET run_matched_silent='mc06'
        --select * from FUNVAR_ADMIN.mutclust_mgmt
        WHERE group_id ='mc05';

        UPDATE FUNVAR_ADMIN.mutclust_mgmt
        SET group_id  ='mc06'
        --select * from FUNVAR_ADMIN.mutclust_mgmt
        WHERE group_id ='dummy';

        UPDATE FUNVAR_ADMIN.mutclust_mgmt
        SET run_matched_silent='mc06'
        --select * from FUNVAR_ADMIN.mutclust_mgmt
        WHERE group_id ='mc06';

        UPDATE FUNVAR_ADMIN.mutclust_mgmt
        SET mfc_mut_count=0
        --select * from FUNVAR_ADMIN.mutclust_mgmt
        WHERE group_id ='mc06' and mfc_mut_count is null;

    -- 31/10/19
        SELECT DISTINCT group_id FROM FUNVAR_ADMIN.mutclust_mgmt;
        -- mfc_mgmt_swisscount_missense_5_10000_20191014.tsv
        
        UPDATE FUNVAR_ADMIN.mutclust_mgmt
        SET group_id ='mc07'
        --select * FROM FUNVAR_ADMIN.mutclust_mgmt
        WHERE group_id='dummy';

        UPDATE FUNVAR_ADMIN.mutclust_mgmt
        SET run_matched_silent='mc08'
        --select * FROM FUNVAR_ADMIN.mutclust_mgmt
        WHERE group_id ='mc07';

        UPDATE FUNVAR_ADMIN.mutclust_mgmt
        SET group_id  ='mc08'
        --select * FROM FUNVAR_ADMIN.mutclust_mgmt
        WHERE group_id ='dummy';

        UPDATE FUNVAR_ADMIN.mutclust_mgmt
        SET run_matched_silent='mc08'
        --select * FROM FUNVAR_ADMIN.mutclust_mgmt
        WHERE group_id ='mc08';

        UPDATE FUNVAR_ADMIN.mutclust_mgmt
        SET mfc_mut_count=0
        --select * FROM FUNVAR_ADMIN.mutclust_mgmt
        WHERE group_id ='mc08' and mfc_mut_count is null;
    

