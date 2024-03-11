-- sql/mutclust_analyse_counts.sql
-- Oct'19

-----------------------------------------------------------------------------------------------------------
-- Oracle SQL views (for ANALYSING MutClust runs)

-- Includes some developmental SELECT only SQL
-----------------------------------------------------------------------------------------------------------

    -- (1) Exploratory SQL
        -- (i) Muts per funfam rep by missense or silent
        -- Note: * EXISTING VW_MUTCLUST_... VIEWS should cover this!! *
            SELECT
                data_source,
                variant_type,
                
                sf_id,
                ff_id,
                rep_id,
            
                vm_synonymous,
                --mm_output_seq_pos,
                count(*) num_mut_ff
            
            FROM
                funvar_tx.mvw_map_core_funfam_pdb
            group by data_source, sf_id, ff_id, rep_id, variant_type, vm_synonymous
            order by rep_id, vm_synonymous;

        -- (ii) What is mutdat sending 
        -- Note: ** Helpful regexp to deal with ora-01722-invalid-number: https://stackoverflow.com/questions/12549029/sql-error-ora-01722-invalid-number
            SELECT 
                *
            FROM
                funvar_tx.vw_mutclust_dat
            where variant_class = 'Missense_Mutation' and vm_synonymous = 'FALSE' AND pdb_res_name <> '-'
            and rep_id = '117eA00'
            ORDER BY rep_id, CAST( regexp_replace( pdb_res_name, '[^0-9]+', '' ) as number )

    
    -- (2) MutClusts per FunFam (4 measures) - MS or SI plus mut_count, gene count and dom length
        -- Row-wise by run_id - i.e. MS then SI (where a pair exists)
        --31/10/19: Note applies to mc05 and mc06 but pnorm() used subsequently....
        SELECT
            mca.groupid,
            mca.runid,
            mca.taskid,
            mcm.superfamily_id,
            mcm.funfam_number,
            mcm.rep_id,
            mcm.mfc_mut_count,
            mcm.num_swissprot_in_ff,
            cd.atom_length,
            sum (case when sig_num_mut_res ='Y' then 1 else 0 end ) as sig_num_mut_res,
            sum (case when sig_mut_count_sum ='Y' then 1 else 0 end ) as sig_mut_count_sum,
            sum (case when sig_weighted_mut_sum ='Y' then 1 else 0 end ) as sig_weighted_mut_sum,
            sum (case when sig_ff_weighted_mut_sum ='Y' then 1 else 0 end ) as sig_ff_weighted_mut_sum,
            mca.radius,
            mca.z_sig_level,
            cd.cath_id
        FROM
            funvar_import.mutclust_analysis mca
        INNER JOIN
            FUNVAR_ADMIN.mutclust_mgmt mcm
            ON mcm.group_id = mca.groupid AND mcm.taskid = mca.taskid
        INNER JOIN 
            funvar_import.CATH_DOMAIN cd
            ON cd.domain_id = mcm.rep_id 
        GROUP BY mca.groupid, mca.runid, mca.taskid, mcm.superfamily_id, mcm.funfam_number, 
        mcm.rep_id, mcm.mfc_mut_count, mcm.num_swissprot_in_ff, cd.atom_length, mca.radius, 
        mca.z_sig_level, cd.cath_id 
        ORDER BY runid;

    -- (3) MutClusts per FunFam summary VIEW 
        --31/10/19: Note applies to mc05 and mc06 but pnorm() used subsequently....
        --(i) VIEW updated 31/10 to use new stats based on pnorm and multiple testing correction (mc07 onwards)
        -- Full export of (available) mc05/mc06 from original view:
        -- ~/woofgit/tracerx/neofun/analysis/NFE/MutClust_MATCHED_analysis_mc05_mc06_by_ROW_20191030.xlsx
            DROP VIEW funvar_tx.VW_MUTCLUST_ANALYSIS_SUMMARY ;
            CREATE VIEW funvar_tx.VW_MUTCLUST_ANALYSIS_SUMMARY 
            AS
               SELECT
                    mca.groupid,
                    mca.runid,
                    mca.taskid,
                    mcm.superfamily_id,
                    mcm.funfam_number,
                    mcm.rep_id,
                    mcm.mfc_mut_count,
                    mcm.num_swissprot_in_ff,
                    (CASE WHEN cgc.cgc_genes IS NULL THEN 0 ELSE cgc.cgc_genes END) num_cgc_genes,
                    cd.atom_length,
                    count(*) row_count,
                   -- sum (case when sig_num_mut_res ='Y' then 1 else 0 end ) as sig_num_mut_res,
                    SUM (case when mca.num_mut_res_p_corr_sig ='Y' then 1 else 0 end ) as sig_num_mut_res_p_corr_sig,
                    
                    --sum (case when sig_mut_count_sum ='Y' then 1 else 0 end ) as sig_mut_count_sum,
                    SUM (case when mca.mut_count_sum_p_corr_sig ='Y' then 1 else 0 end ) as sig_mut_count_sum_p_corr,
                    
                    --sum (case when sig_weighted_mut_sum ='Y' then 1 else 0 end ) as sig_weighted_mut_sum,
                    SUM (case when mca.weighted_mut_sum_p_corr_sig ='Y' then 1 else 0 end ) as sig_weighted_mut_sum_p_corr,
                    
                   --sum (case when sig_ff_weighted_mut_sum ='Y' then 1 else 0 end ) as sig_ff_weighted_mut_sum,
                    SUM (case when mca.ff_weighted_mut_sum_p_corr_sig ='Y' then 1 else 0 end ) as sig_ff_weighted_mut_sum_p_corr,
                    
                    ff.name funfam_name,                    
                    cd.cath_id,
                    mca.radius,
                    mca.p_info,
                    mca.p_correct_method,
                    mca.p_corr_cutoff
                
                FROM
                    funvar_import.mutclust_analysis mca

                INNER JOIN
                    FUNVAR_ADMIN.mutclust_mgmt mcm
                    ON mcm.group_id = mca.groupid AND mcm.taskid = mca.taskid

                INNER JOIN 
                    funvar_import.CATH_DOMAIN cd
                    ON cd.domain_id = mcm.rep_id 

                LEFT JOIN 
                    funvar_import.cath_funfam ff
                    ON ff.superfamily_id = mcm.superfamily_id
                    AND ff.funfam_number = mcm.funfam_number 
                
                -- COSMIC CGC
                LEFT JOIN 
                    funvar_tx.vw_funfam_cgc_gene_count cgc
                    ON cgc.superfamily_id = mcm.superfamily_id
                    AND cgc.funfam_number = mcm.funfam_number
                
                WHERE 
                    -- Old stats calcs in these...
                    groupid NOT IN ( 'mc01','mc02', 'mc03', 'mc04', 'mc05', 'mc06' ) 
                
                GROUP BY mca.groupid, mca.runid, mca.taskid, mcm.superfamily_id, mcm.funfam_number, 
                        mcm.rep_id, mcm.mfc_mut_count, mcm.num_swissprot_in_ff, (CASE WHEN cgc.cgc_genes IS NULL THEN 0 ELSE cgc.cgc_genes END), cd.atom_length, 
                        ff.name, cd.cath_id, mca.radius, mca.p_info, mca.p_correct_method, mca.p_corr_cutoff 
                    
                ORDER BY runid, groupid;

        -- (ii) Join MS and SI using matched IDs
        -- i.e. Column-wise by TASKID (which matches FFs as MS and SI)
        -- Used on mc05 and mc06 (OLD VERSION)
                SELECT
                    ms.groupid,
                    si.groupid,
                    ms.runid,
                    si.runid,
                    ms.taskid,
                    si.taskid,
                    ms.superfamily_id,
                    --si.superfamily_id,
                    ms.funfam_number,
                    --si.funfam_number,
                    
                    ms.rep_id,
                    si.rep_id SI_rep_id,
                    ms.atom_length dom_length,
                    ms.num_swissprot_in_ff,
                    
                    -- Missense muts
                    ms.mfc_mut_count,
                    ROUND( ms.mfc_mut_count /  ms.atom_length, 1 ) avg_ms_per_res,                 
                    -- Silent muts
                    si.mfc_mut_count SI_mfc_mut_count,
                    ROUND( si.mfc_mut_count /  si.atom_length, 1 ) avg_si_per_res,
                    
                    
                    ms.sig_weighted_mut_sum,
                    si.sig_weighted_mut_sum SI_sig_weighted_mut_sum,
                    ms.sig_mut_count_sum,
                    si.sig_mut_count_sum SI_sig_mut_count_sum,
                    ms.sig_num_mut_res,
                    si.sig_num_mut_res SI_sig_num_mut_res,

                    ms.radius,
                    ms.z_sig_level,
                    
                    --ms.cath_id,
                    si.cath_id
                FROM
                    funvar_tx.vw_mutclust_analysis_summary ms 
                
                LEFT JOIN 
                    funvar_tx.vw_mutclust_analysis_summary si
                    ON ms.taskid = si.taskid 

                WHERE 
                    ms.groupid ='mc05' AND si.groupid = 'mc06' ;

        -- (iii) Join MS and SI using matched IDs
        -- i.e. Column-wise by TASKID (which matches FFs as MS and SI)
        -- mc07 and mc08 (and onwards (NEW VERSION))
        -- Updated 07/11/19 to use 'preferred' p-value coumns: MUT_COUNT_SUM_P_CORR, WEIGHTED_MUT_SUM_P_CORR
        -- Used to test HPC runs r/MutFamClust_group_HPC.R - i.e. match between (say) mc07 and tst09_od2 for tasks 100-107
            -- Alternatively use select * from funvar_tx.mutclust_analysis where groupid in ('mc07', 'tst09_od2') and runid = 'R0104MS'....;
                SELECT
                    ms.groupid,
                    si.groupid,
                    ms.runid,
                    si.runid,
                    ms.taskid,
                    si.taskid,
                    ms.superfamily_id,
                    --si.superfamily_id,
                    ms.funfam_number,
                    --si.funfam_number,
                    
                    ms.rep_id,
                    si.rep_id SI_rep_id,
                    ms.atom_length dom_length,
                    ms.num_swissprot_in_ff,
                    ms.num_cgc_genes,
                    
                    -- Missense muts
                    ms.mfc_mut_count,
                    ROUND( ms.mfc_mut_count /  ms.atom_length, 1 ) avg_ms_per_res,                 
                    -- Silent muts
                    si.mfc_mut_count SI_mfc_mut_count,
                    ROUND( si.mfc_mut_count /  si.atom_length, 1 ) avg_si_per_res,
                    
                    
                    ms.sig_weighted_mut_sum_p_corr  MS_sig_weighted_mut_sum_p_corr,
                    si.sig_weighted_mut_sum_p_corr  SI_sig_weighted_mut_sum_p_corr,
                    
                    ms.sig_mut_count_sum_p_corr     MS_sig_mut_count_sum_p_corr,
                    si.sig_mut_count_sum_p_corr     SI_sig_mut_count_sum_p_corr,
                    
                    ms.sig_num_mut_res_p_corr_sig   MS_sig_num_mut_res_p_corr_sig,
                    si.sig_num_mut_res_p_corr_sig   SI_sig_num_mut_res_p_corr_sig,

                    ms.funfam_name,
                    ms.cath_id,
                    ms.radius,
                    ms.p_info,
                    ms.p_correct_method,
                    ms.p_corr_cutoff

                FROM
                    funvar_tx.vw_mutclust_analysis_summary ms 
                
                LEFT JOIN 
                    funvar_tx.vw_mutclust_analysis_summary si
                    ON ms.taskid = si.taskid 

                WHERE 
                    --ms.groupid ='mc07' AND si.groupid = 'mc08' ;
                    ms.groupid ='mc07' AND si.groupid = 'tst09_od2' AND ms.taskid >=100 AND ms.taskid <=107;


            -- (iii) Get Medians, averages and relevant counts
            -- OLD VERSION mc05 mc06
                SELECT
                    ms.groupid,
                    si.groupid,
                    count(ms.groupid) num_ffs,
                    count(si.groupid) si_num_ffs,
            
                    SUM( ms.sig_num_mut_res ) sum_ms_num_mut_res,
                    SUM( si.sig_num_mut_res ) sum_si_num_mut_res,
                    
                    SUM( ms.sig_mut_count_sum ) sum_ms_sig_mut_count_sum,
                    SUM( si.sig_mut_count_sum ) sum_si_sig_mut_count_sum,
                    
                    SUM( ms.sig_weighted_mut_sum ) sum_ms_sig_weighted_mut_sum,
                    SUM( si.sig_weighted_mut_sum ) sum_si_sig_weighted_mut_sum,
                    
                    
                    MEDIAN( ms.mfc_mut_count )  median_ms_mut_count,
                    MEDIAN( si.mfc_mut_count ) median_si_mut_count,
                    ROUND( AVG( ms.num_swissprot_in_ff ), 0) avg_swiss_inf_ff,
                    ROUND( AVG( ms.atom_length ), 0 ) avg_dom_len
                    
        
                FROM
                    funvar_tx.vw_mutclust_analysis_summary ms 
                
                LEFT JOIN 
                    funvar_tx.vw_mutclust_analysis_summary si
                    ON ms.taskid = si.taskid 


                WHERE 
                    ms.groupid ='mc05'                
                    AND si.groupid = 'mc06' 
                                    
                    -- MORE MS THAN SI 
                    AND ms.sig_num_mut_res > si.sig_num_mut_res
                --AND MS.SIG_WEIGHTED_MUT_SUM >  si.sig_weighted_mut_sum
                    
                    -- MORE SI THAN MS
                    --AND SI.SIG_WEIGHTED_MUT_SUM >  ms.sig_weighted_mut_sum 
                    group by ms.groupid, si.groupid; 
            
    -- (4) MutClusts - analysis summary expanded out to include MutClust residues (for NFEs!)
    -- 01/11/2019
        -- Latest missense 10k MutClusts
        SELECT
            mas.groupid,
            mas.runid,
            mas.taskid,
            mas.superfamily_id,
            mas.funfam_number,
            mas.rep_id,
            
            mca.residue,
            mca.mut_count_sum_p,
            mca.mut_count_sum_p_corr,
            mca.mut_count_sum_p_corr_sig,
            mca.weighted_mut_sum_p,
            mca.weighted_mut_sum_p_corr,
            mca.weighted_mut_sum_p_corr_sig,
            
            mas.mfc_mut_count,
            mas.num_swissprot_in_ff,
            mas.num_cgc_genes,
            mas.atom_length,
            mas.row_count,
            --mas.sig_num_mut_res_p_corr_sig,
            mas.sig_mut_count_sum_p_corr,
            mas.sig_weighted_mut_sum_p_corr,
            --mas.sig_ff_weighted_mut_sum_p_corr,
            mas.funfam_name,
            mas.cath_id,
            mas.radius,
            --mas.p_info,
            mas.p_correct_method,
            mas.p_corr_cutoff
        FROM
            funvar_tx.vw_mutclust_analysis_summary mas

        INNER JOIN
            FUNVAR_IMPORT.mutclust_analysis mca
            ON mca.groupid = mas.groupid
            AND mca.runid = mas.runid

        WHERE 
            --mas.groupid ='mc07'
            ( mca.mut_count_sum_p_corr < 0.05 OR mca.weighted_mut_sum_p_corr < 0.05 )
        ORDER BY mas.groupid, mas.runid, mca.weighted_mut_sum_p_corr;
            
            