-- mutclust_dat_af.sql   
-- 12 06 2024
-- AlphaFold rep version of mutclust dat funvar_tx.VW_MUTCLUST_DAT_AF
-- includes FS-3D check to ensure FunFam has SC90 sites
-- NOTE: This has been changed to LEFT JOIN as FSITE tables presently only have
-- sites for PDB reps

-- ***************************************************************************************************************************
-- 20 06 2024 NOTE: needed to update FunVar/Oracle-db 
-- Sequence funfam reps were being returned for each of the MEMBER UNIPROT RANGES, not the single funfam rep_id range!
-- needed to import table with member_ids for non-human uniprots, if they were reps for a FunFam with >0 human proteins - see:
-- /Users/ash/Dropbox/bioinf/neofun/run/mutclust_af/dat/cath_v_4_2_funfam_with_member_ids_all_ff_with_a_human_protein.sql
-- FUNVAR_IMPORT.CATH_FUNFAM_SEQ_REPS (import_id "CFSR" on 20 06 2024)
-- ***************************************************************************************************************************

-----------------------------------------------------------------------------------------------------------
-- Oracle SQL views: [**SEE NOTE ABOVE 20/06/24**]
--  management (mgmt) views:    mutclust/sql/mutclust_mgmt.sql
--  data views:                 mutclust/sql/mutclust_dat.sql
--  data view AF:               mutclust/sql/mutclust_dat_af.sql

-- Runs (i.e. selecting correct mgmt/data for specific runs):
-- Moved to run directory e.g.  ~/Dropbox/bioinf/neofun/run/mutclust/dt12a/mc03/sql/mc03.sql
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
        
      -- (a) MutDat AF (all)
        -- 12 06 2024
           	--DROP VIEW   funvar_tx.VW_MUTCLUST_DAT_AF;
        	CREATE OR REPLACE VIEW funvar_tx.VW_MUTCLUST_DAT_AF AS
                SELECT
                    mcc.data_source,
                    mcc.vm_synonymous,
                    mcc.variant_class,
                    mcc.SUPERFAMILY_ID,
                    mcc.FUNFAM_NUMBER,
                    '99' EF,
                    mcc.vm_gene GENE,
                    mcc.vm_uniprot_accession FUN_UP,
					mcc.rep_source_id,
                    mcc.rep_id REP_ID,
					-- Display either PDB rep as is or bare "uniprot rep" for AF models
					CASE 
						WHEN mcc.rep_source_id ='uniprot' THEN 
							SUBSTR(mcc.rep_id, 1, INSTR( mcc.rep_id, '_' ) -1 ) 
	   					ELSE 
							'NA use VW_MUTCLUST_DAT'
					END AS cath_rep_or_uniprot,
					-- Legacy "PDB_CODE" will default to bare uniprot REP for AF
					CASE 
						WHEN mcc.rep_source_id ='uniprot' THEN 
							SUBSTR(mcc.rep_id, 1, INSTR( mcc.rep_id, '_' ) -1 ) 
	   					ELSE 
							'NA use VW_MUTCLUST_DAT'
					END AS pdb_code,
					-- Legacy "PDB_CHAIN" will default to "-" for AF
					CASE 
						WHEN mcc.rep_source_id ='uniprot' THEN 
							'-'
	   					ELSE 
							'NA use VW_MUTCLUST_DAT'
					END AS chain_code,
                    
                    mcc.vm_seq_no AA_POS,
                    mcc.vm_seq_no PDB_RES_NAME,
                    smc.num_ff_members NUM_SWISSPROT_IN_FF,
					COUNT( DISTINCT mutation_id ) MFC_MUT_COUNT_RESIDUE,
                    COUNT( * ) ROW_COUNT

                FROM
					funvar_tx.mvw_map_core_funfam mcc
                    -- funvar_tx.mvw_map_core_funfam_pdb  mcc
					-- 20 06 2024
                LEFT JOIN 
                (
                    SELECT DISTINCT
						'Y' AS SC90_fsite_data_found, 
                        sf_id,
                        ff_id,
                        rep_id
                    FROM
                        funvar_import.fsites_dist fsd
                    WHERE
                        fsd.import_id = 'IFD03'
						AND fsd.fsite_type = 'SCONS_90'
						--( SELECT import_id FROM FUNVAR_ADMIN.imports WHERE type = 'FSD_DIST' AND archive = 0 ) AND radius = 0
                ) fsd
                    ON fsd.sf_id = mcc.SUPERFAMILY_ID
                    AND fsd.ff_id = to_number( mcc.FUNFAM_NUMBER )
                    AND fsd.rep_id = mcc.rep_id
                LEFT JOIN
                    funvar_tx.vw_funfam_swiss_member_count smc
                    ON 
                         mcc.SUPERFAMILY_ID = smc.superfamily_id 
                    AND
                        mcc.FUNFAM_NUMBER = smc.funfam_number 
                
                WHERE 
					(
						(variant_class = 'Missense_Mutation' AND vm_synonymous='FALSE')
						-- OR
						-- (variant_class = 'Silent' AND vm_synonymous='TRUE')
					)	
					AND
					( 
						-- (data_source = 'Tx' AND data_import_id IN ('ITX03','ITX04'))
						-- OR
						(data_source = 'TCGA' AND data_import_id IN ('ITC02'))
					)
						AND mcc.rep_source_id = 'uniprot'
                     
                GROUP BY 
                    mcc.data_source, 
                    mcc.vm_synonymous, 
                    mcc.variant_class, 
					mcc.SUPERFAMILY_ID,
                    mcc.FUNFAM_NUMBER,
                    '99', 
                    mcc.vm_gene, 
                    mcc.vm_uniprot_accession, 
					mcc.rep_source_id,
                    mcc.rep_id, 
					CASE 
						WHEN mcc.rep_source_id ='uniprot' THEN 
							SUBSTR(mcc.rep_id, 1, INSTR( mcc.rep_id, '_' ) -1 ) 
	   					ELSE 
							mcc.rep_id 
						END,
                    CASE 
						WHEN mcc.rep_source_id ='uniprot' THEN 
							SUBSTR(mcc.rep_id, 1, INSTR( mcc.rep_id, '_' ) -1 ) 
	   					ELSE 
							'NA use VW_MUTCLUST_DAT'
					END,
					-- Legacy "PDB_CHAIN" will default to "-" for AF
					CASE 
						WHEN mcc.rep_source_id ='uniprot' THEN 
							'-'
	   					ELSE 
							'NA use VW_MUTCLUST_DAT'
					END,
                    mcc.vm_seq_no, 
                    mcc.vm_seq_no,
                    smc.num_ff_members
                
                ORDER BY
					mcc.SUPERFAMILY_ID,
                    mcc.FUNFAM_NUMBER,
                    mcc.rep_id, 
                    mcc.vm_seq_no;

-- CATH and UNIPROT reps for AF processing
-- COUNT by FunFam repid (cath and uniprot)
-- Ash / 12 03 2024
SELECT
	data_source,
	CASE WHEN data_import_id ='ITX04' THEN 'Y'
	   ELSE 'N' END AS TXP,
	data_import_id,
	vm_import_id,	
	rep_source_id,
	
	--  include for: modules/neofun/af/mutation_counts_sql_and_results/mutation_counts_rep_ids_inc_singletons_variant_class.tsv
	--	variant_class,
	
	superfamily_id,
	funfam_number,
	rep_id,
	CASE WHEN rep_source_id ='uniprot' THEN 
		SUBSTR(rep_id, 1, INSTR( rep_id, '_' ) -1 ) 
	   ELSE rep_id END AS cath_rep_or_uniprot,
	
	COUNT( DISTINCT rep_id ) AS num_ff_reps,
	COUNT( DISTINCT mutation_id ) AS num_muts,
	COUNT( * ) AS num_rows,
	
--	vm_synonymous,
	ffr_import_id
	
FROM
	funvar_tx.mvw_map_core_funfam
WHERE
	-- INCLUDE ALL SNPS
	(
		(variant_class = 'Missense_Mutation' AND vm_synonymous='FALSE')
		OR
		(variant_class = 'Silent' AND vm_synonymous='TRUE')
	)	
	AND
	( 
		(data_source = 'Tx' AND data_import_id IN ('ITX03','ITX04'))
		OR
		(data_source = 'TCGA' AND data_import_id IN ('ITC02'))
	)
	--	AND rep_source_id = 'uniprot'
GROUP BY
	data_source, 
	CASE WHEN data_import_id ='ITX04' THEN 'Y' ELSE 'N' END, 
	data_import_id, 
	vm_import_id, 
	ffr_import_id, 
--	vm_synonymous, 

--  include for: modules/neofun/af/mutation_counts_sql_and_results/mutation_counts_rep_ids_inc_singletons_variant_class.tsv
--	variant_class, 

	rep_source_id,

	superfamily_id, 
	funfam_number, 
	rep_id, 
	CASE WHEN rep_source_id ='uniprot' THEN SUBSTR(rep_id, 1, INSTR( rep_id, '_' ) -1 ) ELSE rep_id END

-- remove for: modules/neofun/af/mutation_counts_sql_and_results/mutation_counts_rep_ids_inc_singletons_variant_class.tsv
HAVING
	COUNT( DISTINCT mutation_id ) >1

ORDER BY
	data_source,
	data_import_id,
	vm_import_id,
	ffr_import_id,

	--  include for: modules/neofun/af/mutation_counts_sql_and_results/mutation_counts_rep_ids_inc_singletons_variant_class.tsv
	--	variant_class,
	
	rep_source_id,
	COUNT( DISTINCT mutation_id ) DESC

--	superfamily_id,
--	funfam_number,
--	
--	rep_id,

--	vm_synonymous,
	
	;
	
	-- mut props etc
--	mutation_id,
--	vm_uniprot_accession,
--	aa_range_low,
--	vm_seq_no,
--	aa_range_high,
--	vm_uniprot_aa,
--	vm_aa_change,
--	ffr_gene_name,
--	source_hugo_symbol,
--	vm_gene,
--	ffr_gene_id,
--	vm_gene_acc,
--	vm_refseq_gene_acc,
--	vm_refseq_transcript,
--	variant_type,
--	vm_change_type,
--	ffr_member_id,
--	ffr_num_funfam_ranges,
--	vm_note,
--hid
	
	
	