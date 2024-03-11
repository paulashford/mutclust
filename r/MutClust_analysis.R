# MutClust_analysis.R
# 23/10/2019
# Main analysis of results returned from randomisations trials
# Renamed from: MutClust_analysis_by_size.R
# Placed into callable function, see: MutClust_analysis_example.R for calls

# Note: Only single significance value and radius can be passed per call, to simplify 
# nested loops in previous versions

# Analyse the results of a group of MutClust runs to identify significant clusters using the 4 metrics:
#	1) num_mut_res	: How many mutated residues nearby?
#	2) MFC_MUT_COUNT_SUM is summed mutations for residue over all genes
#	3) WEIGHTED_MUT_SUM applies a Hill curve normalisation to 'compress' (with theta=2, exp=3 means ~ [0,6])
# 4) FF_WEIGHTED_MUT_SUM divides the weighted sum by number of swissprot human paralogues in FunFam to
#		account for bias from FunFams with many members

# From randomisation trial results, significant residues are identified 
# as: (actual counts) > mean(all trials) +  z * sd(all trials) 
# for each of the 4 metrics

# Oct '19 Multiple changes reflect 4 calculations now returned by MutClust (e.g. Hill normalised)
# See /Users/ash/woofgit/mutclust/r/struct_clust_randomise.R

  sca.analyse_mutclust_runs <- function(  mgmt_file,
                                          groups_dir,
                                          groupID,
                                          batchCode = "MS",
                                          radius_used = 5,
                                          num_trials = 10000L,
                                          p_correct_method = "BH", 
                                          p_corr_cutoff  = 0.05,
                                          taskID_list = "ALL" ){

    # ----------
    # Initialise
    # ----------
      library( data.table );
      library( bio3d );
      library( fields );
      library( plyr );
      
      #library(Biostrings);
      # 27/06/18 new MutClust git directory
      setwd( '~/woofgit/mutclust' );
      source( 'r/MutClust_analysis_functions.R' );
      source( 'r/struct_clust_util.R' );
      source( 'r/struct_clust_load_data.R' );
      #source( '/Users/ash/woofgit/tracerx/neofun/script/R/neofun_util.R' );

      # Where all MutClust RData files are...
      group_dir <- file.path( groups_dir, groupID );
      dat_file_tail <- paste( "_", radius_used, "-", radius_used, "_", as.character( num_trials ), "_dtscoresALL.RData", sep = "" );

    # *********************************
    # Get the management file
    # *********************************
      dt_mgmt <- data.table( read.table( mgmt_file, header = TRUE, sep = "\t", quote = "", dec = ".", stringsAsFactors = FALSE ) );

    # *********************************
    # Loop tasks
    # *********************************
      # Either all tasks in mgmt file, or the passed in list
      if ( taskID_list[1] == "ALL" ){
        taskID_list <- sort( unique( dt_mgmt$TASKID ) );
      }
    
      for ( iTask in taskID_list  ) {
      
        # Clear any existing data
        if ( exists( "dt_scoresALL" ) ) { rm( dt_scoresALL ) };

        # Get task info by filtering the mgmt table
        dt_task <- dt_mgmt[ TASKID == iTask ];
        # Single row?
        
        if ( nrow( dt_task ) != 1 ){
          stop( writeLines( paste0( "#ERROR: Unable to select single TASKID from group management tsv file for task number ", iTask ) ) );
        }
  
        # Generate runID based on taskID & batchCode
        runID <- paste0( "R", gsub( " ", "0", sprintf( "%4s", iTask ) ), batchCode );

        # Echo some info for this run
        writeLines( paste0( ">> MutClust_analysis for taskID: ", iTask, 
                    " groupID: ", groupID, 
                    " runID: ", runID,
                    " batchCode:", batchCode,
                    " radius: ", radius_used
                    )
                  );
  
        # Load MutClust results data (observations & random trials table)
        #eg: mc06_R0151SI_5-5_1000_dtscoresALL.RData
        dat_file <- file.path( group_dir, paste0( groupID, "_",  runID, dat_file_tail ) );
        writeLines( paste0( "Using results file: ", dat_file ) );
  
        if ( file.exists( dat_file ) ) 
          {
            # MutClust file found, so load and analyse...
            load( dat_file );
    
            # Found file but is empty - something went wrong in MutClust
            if ( nrow( dt_scoresALL ) == 0 ){
              writeLines( paste0( "#ERROR: No suitable MFC data to analyse for ", dat_file ) ) ;
            }
    
            # Get means and std dev for MFC data across all 4 measures 
            # Calculates p-values and cut-offs using pnorm()
            dt_meansd <- sca.getMeanDT( dt_scoresALL, radius_used, p_correct_method, p_corr_cutoff );

            # Which residues have more mutated residues than expected: obs_muts > (mean + z_sig * SD)
            ## sca.getHighOutliers( dt_meansd, z_sig_level, type = "num_mut_res" );
            # These are now flagged in new table format e.g. 
            # dt_meansd[ weighted_mut_sum_p_corr_sig == "Y" ]"weighted_mut_sum_p_corr_sig == "Y"
            #
            # For TEST using Chimera
            #paste( dt_meansd[ sig_weighted_mut_sum == 'Y', residue ], collapse = ", " );
            
            # For analysis, just residues with significance by any of the 4 scores
            #dt_mutclust_analysis_task <- dt_meansd[ sig_num_mut_res == 'Y' | sig_mut_count_sum == 'Y' | sig_weighted_mut_sum == 'Y' | sig_ff_weighted_mut_sum == 'Y' ];

            if ( nrow ( dt_meansd ) > 0 ) {
              # Add essential ID columns and multiple testing correction method/cut-off info
              dt_meansd[ , groupid := groupID ];
              dt_meansd[ , runid := runID ];
              dt_meansd[ , taskid := iTask ];
              dt_meansd[ , p_correct_method := p_correct_method ];
              dt_meansd[ , p_corr_cutoff := p_corr_cutoff ];
            
              # Create or update analysis table
              if ( !is.defined( dt_mutclust_analysis ) ) {
                dt_mutclust_analysis <- dt_meansd;
              } else {
                dt_mutclust_analysis <- rbind ( dt_mutclust_analysis, dt_meansd );
              }

            } # nrow (dt_mutclust_analysis_task) > 0

          } else {
            print ( paste0( "MutClust file not found: ", dat_file ) );
          }

      } # iTask


    # Save full dt as RData
    file_mutclust_analysis <- file.path( groups_dir, 
                                        paste0( "MutClust_analysis_group_", groupID, 
                                        "_", batchCode, 
                                        "_r_", radius_used, 
                                        "_trials_", num_trials, 
                                        "_pval-cutoff_", p_corr_cutoff,
                                        "_pval-correction_", p_correct_method
                                        )                                        
                                      );
    save( dt_mutclust_analysis, file = paste0( file_mutclust_analysis, ".RData" ) );

    # Save tsv form for Oracle import
    # Note issues with small values exporting scientific notation "9.371231421123e-05" etc
    write.table( dt_mutclust_analysis,
					file = paste0( file_mutclust_analysis, ".tsv" ), 
					quote = FALSE, row.names = FALSE, sep = "\t", fileEncoding = "UTF-8" );

    write.table(  dt_mutclust_analysis[ ,
                                        .(  trialNo,
                                            index,  
                                            residue,
                                            radius, 
                                            num_mut_res,                   
                                            MFC_MUT_COUNT_SUM,             
                                            WEIGHTED_MUT_SUM = format( WEIGHTED_MUT_SUM, scientific = FALSE ),
                                            FF_WEIGHTED_MUT_SUM = format( FF_WEIGHTED_MUT_SUM, scientific = FALSE ),
                                            num_mut_res_mean = format( num_mut_res_mean, scientific = FALSE ),
                                            num_mut_res_sd = format( num_mut_res_sd, scientific = FALSE ),
                                            mut_count_sum_mean = format( mut_count_sum_mean, scientific = FALSE ),
                                            mut_count_sum_sd = format( mut_count_sum_sd, scientific = FALSE ) ,  
                                            weighted_mut_sum_mean = format( weighted_mut_sum_mean, scientific = FALSE ),         
                                            weighted_mut_sum_sd = format( weighted_mut_sum_sd, scientific = FALSE ),
                                            ff_weighted_mut_sum_mean = format( ff_weighted_mut_sum_mean, scientific = FALSE ),      
                                            ff_weighted_mut_sum_sd = format( ff_weighted_mut_sum_sd, scientific = FALSE ),        
                                            num_mut_res_p = format( num_mut_res_p, scientific = FALSE ),
                                            num_mut_res_p_corr= format( num_mut_res_p_corr, scientific = FALSE ),            
                                            num_mut_res_p_corr_sig,
                                            mut_count_sum_p = format( mut_count_sum_p, scientific = FALSE ),         
                                            mut_count_sum_p_corr = format( mut_count_sum_p_corr, scientific = FALSE ),
                                            mut_count_sum_p_corr_sig,
                                            weighted_mut_sum_p = format( weighted_mut_sum_p, scientific = FALSE ),
                                            weighted_mut_sum_p_corr = format( weighted_mut_sum_p_corr, scientific = FALSE ), 
                                            weighted_mut_sum_p_corr_sig, 
                                            ff_weighted_mut_sum_p = format( ff_weighted_mut_sum_p, scientific = FALSE ), 
                                            ff_weighted_mut_sum_p_corr = format( ff_weighted_mut_sum_p_corr, scientific = FALSE ), 
                                            ff_weighted_mut_sum_p_corr_sig,
                                            p_info,
                                            groupid,
                                            runid,
                                            taskid, 
                                            p_correct_method,
                                            p_corr_cutoff 
                                        )
                                      ],
          file = paste0( file_mutclust_analysis, ".tsv" ), 
					quote = FALSE, row.names = FALSE, sep = "\t", fileEncoding = "UTF-8" );         

}