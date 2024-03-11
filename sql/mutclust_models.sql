-- MutClust models VIEWS
-- 31/05/2020

-- (1) mutclust_dat_models
    DROP VIEW funvar_tx.vw_mutclust_dat_models;
    CREATE VIEW funvar_tx.vw_mutclust_dat_models
    AS
      SELECT
            --data_source
            mcf.data_source,
            --vm_synonymous
            mcf.vm_synonymous,
            --variant_class
            mcf.variant_class,
            --superfamily_id
            mcf.superfamily_id,
            --funfam_number
            mcf.funfam_number,
            --ef
            '99' AS ef,
            --gene
            mcf.vm_gene AS gene,
            --fun_up
            mcf.vm_uniprot_accession AS fun_up,
            --rep_id
            CONCAT(mcf.rep_id, '.pdb') rep_id,
            --pdb_code
            --mcf.rep_source_id AS pdb_code,
            mcf.rep_id AS pdb_code,
            --chain_code
            '-' AS chain_code,
            --aa_pos
            vm_seq_no AS aa_pos,
            --pdb_res_name,
            vm_seq_no AS pdb_res_name,
            smc.num_ff_members AS num_swissprot_in_ff,
            -- Count is for each FunFam rep residue ** per gene and uniprot **
            -- struct_clust_runner.R does SUM as :
                --# Count mutations per residue (i.e. over all genes) 
	            --dt_data_mutfreq[ ,  MFC_MUT_COUNT_SUM := sum( MFC_MUT_COUNT_RESIDUE ), by = PDB_RES_NAME ];
            COUNT(*) AS mfc_mut_count_residue
        FROM
            --funvar_tx.vw_mutclust_dat
            FUNVAR_TX.mvw_map_core_funfam mcf

        LEFT JOIN
                funvar_tx.vw_funfam_swiss_member_count smc
            ON 
                mcf.superfamily_id = smc.superfamily_id 
            AND
                mcf.funfam_number = smc.funfam_number 
        WHERE 
            --  mcf.data_source = 'TCGA'
                --AND variant_class IN ( 'Missense_Mutation', 'Silent' )    -- TCGA_mutdat_sequence_reps_sql007_S1_20200529.tsv
            -- AND mcf.variant_class IN ( 'Missense_Mutation' )    -- TCGA_mutdat_sequence_reps_missense_sql007_S1_20200529.tsv
                --AND variant_class IN ( 'Silent' ) -- TCGA_mutdat_sequence_reps_silent_sql007_S1_20200529.tsv
                --AND  vm_aa_change NOT LIKE '%/*'   -- E.g. downstream missense vars E/* etc - no good for cath-mutation-mapper.pl [NB none in Tx404 data)
                mcf.rep_source_id = 'uniprot'     
        GROUP BY 
            mcf.data_source, mcf.vm_synonymous, mcf.variant_class, mcf.superfamily_id, mcf.funfam_number, 
            '99', mcf.vm_gene, mcf.vm_uniprot_accession, CONCAT(mcf.rep_id, '.pdb'), mcf.rep_id, 
            '-', vm_seq_no, smc.num_ff_members;

-- (2) mutclust_mgmt_missense_models
    DROP VIEW funvar_tx.vw_mutclust_mgmt_missense_mod;
    CREATE VIEW funvar_tx.vw_mutclust_mgmt_missense_mod
    AS
    SELECT 
        ROW_NUMBER( ) OVER( ORDER BY MFC_MUT_COUNT_MISSENSE DESC ) TASKID,
        superfamily_id,
        funfam_number,
        ef,
        gene,
        fun_up,
        rep_id,
        SUBSTR(pdb_code, 1,INSTR(pdb_code,'_')-1) rep_up,
        SUBSTR(pdb_code, INSTR(pdb_code,'_')+1 )  rep_aa_range,
        pdb_code,    -- rep_id and pdb_code are same for models
        chain_code,         -- If '-' then struct_clust_load_data func: sc.load_MutFamClust_freq_tsv will not filter using chain as models don't have chain ID
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
                funvar_tx.vw_mutclust_dat_models
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

-- (3) Generate MGMT file with radii and trials counts
    -- (a) ALL missense
        SELECT
            taskid,
            superfamily_id,
            funfam_number,
            ef,
            gene,
            fun_up,
            rep_id,
            pdb_code,
            chain_code,
            num_swissprot_in_ff,
            mfc_mut_count_missense MFC_MUT_COUNT,
            5 MIN_RAD,
            5 MAX_RAD,	
            10000 NO_TRIALS
        FROM
            funvar_tx.vw_mutclust_mgmt_missense_mod;

    -- (b) Generate MGMT file for EXACT MATCH (amino acid range of rep <-> model) for missense
    -- i.e. where amino start/end pos are same as CATH FF sequence rep.
        SELECT
            -- If doing subset on HPC, need to use an index for the subset (2nd col) 
            -- or HPC will only run say, 3 out of 16 jobs in a batch if 13 have no models...
            taskid taskid_idx,
            ROW_NUMBER( ) OVER( ORDER BY mfc_mut_count_missense DESC ) TASKID,

            superfamily_id,
            funfam_number,
            ef,
            gene,
            fun_up,
            rep_id,
            --rep_up,
            --rep_aa_range,
            pdb_code,
            --CONCAT(rep_id,'.pdb') pdb_code,
            chain_code,
            num_swissprot_in_ff,
            mfc_mut_count_missense MFC_MUT_COUNT,
            --mod.exact_match
            5 MIN_RAD,
            5 MAX_RAD,	
            10000 NO_TRIALS
        FROM
            funvar_tx.vw_mutclust_mgmt_missense_mod
        INNER JOIN 
            funvar_import.models mod
            ON rep_id = mod.model_file
        ORDER BY mfc_mut_count_missense DESC;

    -- (4) Generate MUTDAT file 
        -- (a) Missense TCGA
            -- FORMAT:
            --~/Dropbox/bioinf/neofun/run/mutclust/dt12a/mutdat$ head mfc_mutdat_missense_20191005.tsv
            -- SUPERFAMILY_ID	FUNFAM_NUMBER	EF	GENE	FUN_UP	REP_ID	PDB_CODE	CHAIN_CODE	AA_POS	PDB_RES_NAME	MFC_MUT_COUNT
            -- 3.30.70.270	50164	99	POLH	Q9Y253	3mr2A01	3mr2	A	159	159	1
            -- 3.40.50.2000	95564	99	PYGM	P11217	1em6A01	1em6	A	478	477	2
            -- 2.30.42.10	23703	99	PTPN4	P29074	2cs5A01	2cs5	A	572	73	1
            -- 3.80.10.10	105889	99	NTRK1	P04629	2ifgA01	2ifg	A	110	110	1
            SELECT
                --data_source,
                --vm_synonymous,
                --variant_class,
                superfamily_id,
                funfam_number,
                99 ef,
                gene,
                fun_up,
                rep_id,
                pdb_code,
                chain_code,
                aa_pos,
                pdb_res_name,
                --num_swissprot_in_ff,
                mfc_mut_count_residue
            FROM
                funvar_tx.vw_mutclust_dat_models
            WHERE 
                data_source = 'TCGA'
                AND
                vm_synonymous = 'FALSE'
                AND
                variant_class = 'Missense_Mutation';

-- (5) How many FunFams have UniProt rep?
    -- 6,629
    SELECT COUNT (*) FROM
    (SELECT DISTINCT 
        mm.superfamily_id,
        mm.funfam_number
        FROM
        funvar_tx.vw_mutclust_mgmt_missense_mod mm);

    -- How many EXACT matches for start/end pos of rep/model?
    -- 2,160
    SELECT COUNT (*) FROM
    (SELECT DISTINCT 
        superfamily_id,
        funfam_number
        FROM
        (
        --... <use SQL 6> ...
        ))

    -- Update EXACT_MATCH on models
    -- 2,198 rows updated.
    UPDATE     
    funvar_import.models
    SET exact_match =1 
    WHERE model_file IN
    (
    SELECT DISTINCT model_file FROM(
         (
        --... <use SQL 6> ...
        ))
--------------------------------------------------------
-- (6) Reps and Models and aa ranges...
--------------------------------------------------------
    SELECT
        --taskid,
        mm.superfamily_id,
        mm.funfam_number,
        mm.ef,
        mm.gene,
        mm.fun_up,
        
        mm.rep_id,
        mm.rep_up,
        mm.rep_aa_range,
        SUBSTR(mm.rep_aa_range,1,INSTR(mm.rep_aa_range,'_')-1) rep_aa_low,
        SUBSTR(mm.rep_aa_range,INSTR(mm.rep_aa_range,'_')+1) rep_aa_high,
        
        mod.model_uniprot,
        mod.startpos,
        mod.endpos,
        mod.model_file,
        mod.import_id,
        
        mm.pdb_code,
        mm.chain_code,
        mm.num_swissprot_in_ff,
        mm.mfc_mut_count_missense
    FROM
        funvar_tx.vw_mutclust_mgmt_missense_mod mm
    LEFT JOIN 
        FUNVAR_IMPORT.models mod
        ON mm.rep_up = mod.model_uniprot

    WHERE 
  
        -- -- EXACT MATCH
        -- mod.startpos =  SUBSTR(mm.rep_aa_range,1,INSTR(mm.rep_aa_range,'_')-1)
        -- AND
        -- mod.endpos = SUBSTR(mm.rep_aa_range,INSTR(mm.rep_aa_range,'_')+1) 

        -- REASONABLE match (i.e. model end is not before FF start or start after end...)
        -- model end after funfam start_pos
        mod.endpos > SUBSTR(mm.rep_aa_range,1,INSTR(mm.rep_aa_range,'_')-1)
        AND 
        -- model stat before funfam end_pos
        mod.startpos < SUBSTR(mm.rep_aa_range,INSTR(mm.rep_aa_range,'_')+1) 

    ORDER BY
        mm.superfamily_id,
        mm.funfam_number,
        mm.rep_up,
        SUBSTR(mm.rep_aa_range,1,INSTR(mm.rep_aa_range,'_')-1),
        mod.startpos
        ;


--------------------------------------------------------
-- (7) COPY PDB files to TEST directory
--------------------------------------------------------
    SELECT
        taskid,
        CONCAT( CONCAT( CONCAT('cp /home/ucbtshf/run/neofun/dt12b/models/', rep_id), '.pdb' ), ' /home/ucbtshf/run/neofun/dt12b/models_test' ) CMD,
        mfc_mut_count_missense MFC_MUT_COUNT
    FROM
        funvar_tx.vw_mutclust_mgmt_missense_mod
    INNER JOIN 
        funvar_import.models mod
        ON CONCAT( rep_id, '.pdb') = mod.model_file
    ORDER BY mfc_mut_count_missense DESC;
        