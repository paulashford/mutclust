# MutClust_results_analysis_runner.R
# 30/10/2019

    # This runs analysis of MutClust by setting up various parameters required and 
    # passing to relevant functions
    # MATCHED pairs of missene / silent are grouped here too, if available

    # Analysis script dependencies
    # MutClust_analysis.R
    # struct_clust_analysis_functions.R

    # ----------------------------------------------------------------------
    # For collation of mc07 and mc08 batches from subdirectories (Nov 2019)
    # see: sh/MutClust_collate_pchuckle_batches.sh
    # ----------------------------------------------------------------------
    
    # wd
    setwd( "~/woofgit/mutclust" ); 
    source( "r/MutClust_analysis.R" );

    # Files & paths
    dir_mutdat <- "/Users/ash/Dropbox/bioinf/neofun/run/mutclust/dt12a/mutdat";
    groups_dir <- "/Users/ash/Dropbox/bioinf/neofun/run/mutclust/dt12a";

    # MATCHED 1k randomisation group ~124 matched
        # Group ID / mgmt / batch / radius / trials / sig
        # num_trials  <- 1000L;
        # radius_used <- 5;
        # # mc05 (MS)
        #     groupID     <- "mc05";
        #     mgmt_file   <- file.path( dir_mutdat, "mfc_mgmt_swisscount_missense_5_1000_20191014.tsv" );
        #     batchCode   <- "MS";
        # # mc06 (SI)
        #     groupID     <- "mc06";
        #     mgmt_file   <- file.path( dir_mutdat, "mfc_mgmt_swisscount_silent_5_1000_20191014.tsv" );
        #     batchCode   <- "SI";

    # MATCHED 10k randomisation group ~100
        #https://serverfault.com/questions/141773/what-is-archive-mode-in-rsync
        # get data from pchuckle
        # Move complete (RData exists) from hpc directories to mc07
       
        # ** see sh/MutClust_collate_pchuckle_batches.sh **
        # ** see also _MutClust_CS_cluster_tmux_cores_summary.xlsx **
        
        # Group ID / mgmt / batch / radius / trials / sig
        num_trials  <- 10000L;
        radius_used <- 5;
        # mc07 (MS)
            groupID     <- "mc07";
            mgmt_file   <- file.path( dir_mutdat, "mfc_mgmt_swisscount_missense_5_10000_20191014.tsv" );
            batchCode   <- "MS";
           

        # # mc08 (SI)
        #     groupID     <- "mc08";
        #     mgmt_file   <- file.path( dir_mutdat, "mfc_mgmt_swisscount_silent_5_10000_20191014.tsv" );
        #     batchCode   <- "SI";
        #     #  # Only 1-500 submitted NOTE NOT COMPLETE
        #   **  taskID_list <- c(1:500)

    # TEST (e.g. for HPC parallel tst runs
    # Group ID / mgmt / batch / radius / trials / sig
        # num_trials  <- 10000L;
        # radius_used <- 5;
        # # orengodev2 / /db/neofun/mutclust/run/tst09_od2 / tst09_od2 (based on mc07 (MS) )
        #     groups_dir <- "~/Dropbox/bioinf/neofun/run/mutclust/dt12a/tst09_od2";
        #     groupID     <- "tst09_od2";
        #     mgmt_file   <- file.path( dir_mutdat, "mfc_mgmt_swisscount_missense_5_10000_20191014.tsv" );
        #     batchCode   <- "MS";
        #    ## Only 1-500 submitted NOTE NOT COMPLETE
        #    #arg_task_ID="range:100:107"
        #     taskID_list <- c(100:107)
        # # mc08 (SI)
        #     groupID     <- "mc08";
        #     mgmt_file   <- file.path( dir_mutdat, "mfc_mgmt_swisscount_silent_5_10000_20191014.tsv" );
        #     batchCode   <- "SI";
        #     #  # Only 1-500 submitted NOTE NOT COMPLETE
        #     taskID_list <- c(1:500)

    # TEST changes in commit [4111246] for data.table key speed improvements (~50% - 85%)
        
        #  **************************************************************************
        # MAIN VALIDATION [tbc]
        #   /home/ucbtshf/run/neofun/dt12a/mutclust/mc07/hpc_multicore_runs/mc07_mco_adhoc
        #   mc07_R0009-12
        # vs
        #   /db/neofun/mutclust/run/tst12_od2
        #   mc07_R0009-12

        #  **************************************************************************

        # # tst10_od2 [original]
        # # /db/neofun/mutclust/run/tst10_od2
        # # struct_clust_randomise.R [7419ecf]  (original / mc07/hpc_runs...)
        # # rsync -avhP orengodev2_remote:/db/neofun/mutclust/run/tst10_od2/ /Users/ash/Dropbox/bioinf/neofun/run/mutclust/dt12a/tst10_od2
        # num_trials  <- 100L;
        # radius_used <- 5;
        # groups_dir  <- "/Users/ash/Dropbox/bioinf/neofun/run/mutclust/dt12a";
        # groupID     <- "tst10_od2";
        # mgmt_file   <- file.path( dir_mutdat, "mfc_mgmt_swisscount_missense_5_100_20191014.tsv" );
        # #batchCode   <- "7419ecf";
        # batchCode   <- "MS";
        # taskID_list <- c( 500:507 );

        # # tst11_od2 [updated code]
        # # /db/neofun/mutclust/run/tst11_od2
        # # struct_clust_randomise.R [4111246]
        # # rsync -avhP orengodev2_remote:/db/neofun/mutclust/run/tst11_od2/ /Users/ash/Dropbox/bioinf/neofun/run/mutclust/dt12a/tst11_od2
        # num_trials  <- 100L;
        # radius_used <- 5;
        # groups_dir  <- "/Users/ash/Dropbox/bioinf/neofun/run/mutclust/dt12a";
        # groupID     <- "tst11_od2";
        # mgmt_file   <- file.path( dir_mutdat, "mfc_mgmt_swisscount_missense_5_100_20191014.tsv" );
        # #batchCode   <- "4111246";
        # batchCode   <- "MS";
        # taskID_list <- c( 500:507 );
        
        # # Import both to Oracle MUTCLUST_ANALYSIS and MUT_CLUST_MGMT
        # # comparison
        # SELECT * FROM funvar_tx.vw_mutclust_analysis_summary
        # WHERE groupid IN ('tst10_od2', 'tst11_od2')
        # ORDER BY superfamily_id, funfam_number, rep_id, groupid;
        # # Some differences as only 100 trials - the p-values are very similar via 'spot-checking':
        # SELECT
        # groupid, runid, index_num, residue, radius, num_mut_res, mfc_mut_count_sum,
        # num_mut_res_p_corr, mut_count_sum_p_corr, weighted_mut_sum_p, p_correct_method, p_corr_cutoff, p_info
        # FROM
        # funvar_import.mutclust_analysis
        # WHERE groupid IN ('tst10_od2', 'tst11_od2')
        # ORDER BY runid,  residue,groupid;


    # Test at significance....
        # Multiple testing corrections (see R/ ?p.adjust)
        # z_corrections <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
        p_correct_method <- "BH";
        p_corr_cutoff <- 0.05;

        # Tasks to run
        taskID_list <- c( "ALL" );
        #taskID_list <- c( 10:5115 );
        sca.analyse_mutclust_runs(  mgmt_file,
                                    groups_dir,
                                    groupID,
                                    batchCode,
                                    radius_used,
                                    num_trials,
                                    p_correct_method, 
                                    p_corr_cutoff, 
                                    taskID_list 
                                );

    # Import for "NFE3" on 19/11/19
        UPDATE funvar_import.mutclust_analysis
        SET groupid='mc07_nfe1_nfe2'
        WHERE groupid='mc07';
        # Imported ~/Dropbox/bioinf/neofun/run/mutclust/dt12a/MutClust_analysis_group_mc07_MS_r_5_trials_10000_pval-cutoff_0.05_pval-correction_BH.tsv
        SELECT groupid,count(*)
        from funvar_import.mutclust_analysis
        group by groupid;
        # mc05	826
        # mc06	1034
        # mc07	53619
        # mc07_nfe1_nfe2	9723
        # mc08	24165
        # tst09_od2	1111
        # tst10_od2	522
        # tst11_od2	522

          
     