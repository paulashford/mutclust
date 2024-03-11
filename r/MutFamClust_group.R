# MutFamClust_group
# Ash 04/02/16
# Run MutFamClust for a group of mutfams
# Uses qsub tasks to create multiple jobs as per mgmt tsv file

# Oct'19 - updates for explicit input args & clarity

# Takes  arguments:
# 1: MGMT .tsv file for GroupID
# 2: taskID  (indexes 1st column of this mgmt tsv)
# 3: PDB DIR
# 4: Explicit mutdat file
# 5: GroupID e.g. "mc03" for multiple FunFams with same mut data
# 6: Type/batch code - e.g. "MS" missense, "SY" synonymous/silent
# 7: Output dir - where to write all files

library( data.table );
# Data table concurrent threads (for cluster)
# https://jangorecki.gitlab.io/data.cube/library/data.table/html/openmp-utils.html
setDTthreads( 1 );

source( 'struct_clust_main.R' );

# Arguments input
  args <- commandArgs( trailingOnly = TRUE );

  if ( length( args ) != 7 ) {
    stop( "Args required - 1: MGMT file, 2: taskID, 3: PDB_DIR, 4: MUTDAT file, 5: GroupID, 6: Type/batch code (MS: Missense; SI: Silent), 7: Output dir" );
  }

  # Mgmt file
  arg_mgmt_file   <- args[ 1 ];

  # Task ID 
  arg_task_ID     <- args[ 2 ]; 

  # PDB dir
  arg_pdb_dir     <- args[ 3 ];

  # MutDat
  arg_mutdat_file  <- args[ 4 ];

  # GroupID
  arg_groupID     <- args[ 5 ];

  # Type/batch code (e.g.)
  arg_batchCode   <- args[ 6 ];

  # Output dir
  arg_out_dir     <- args[ 7 ];

# Basic validation
    # No mgmt file... 
  if ( file.exists( arg_mgmt_file ) == 0 ){
    stop( paste0( "File ", arg_mgmt_file, " not found." ) );
  }

  # No mutdat file
  if ( file.exists( arg_mutdat_file ) == 0 ){
    stop( paste0( "File ", arg_mutdat_file, " not found." ) );
  }

# Read mgmt file
  dt_mgmt <- data.table( read.table(  arg_mgmt_file, 
                                      header = TRUE, 
                                      sep = "\t", 
                                      quote = "\"",
                                      dec = ".", 
                                      stringsAsFactors = FALSE ) );
                                      
  # Task ID Get just row for specific TASKID
  dt_mgmt       <- dt_mgmt[ TASKID == arg_task_ID ];
  
  # If TASK_ID not in mgmt file:
  if ( nrow( dt_mgmt ) == 0 ){
    stop( paste0( "No task-specific row found in mgmt file ", arg_mgmt_file, " for TASKID ", arg_task_ID, "." ) );
  }

# Extract relevant vars from mgmt file
# Header (Oct '19)
# TASKID	SUPERFAMILY_ID	FUNFAM_NUMBER	EF	GENE	FUN_UP	REP_ID	PDB_CODE	CHAIN_CODE	NUM_SWISSPROT_IN_FF	MFC_MUT_COUNT	MIN_RAD	MAX_RAD	NO_TRIALS
  taskID                <- dt_mgmt[ , TASKID ];
  sf_ID                 <- dt_mgmt[ , SUPERFAMILY_ID ];
  ff_num                <- dt_mgmt[ , FUNFAM_NUMBER ];
  ef                    <- dt_mgmt[ , EF ];
  gene                  <- dt_mgmt[ , GENE ];     # FunFam GENE - i.e. of rep (not currently used)
  fun_up                <- dt_mgmt[ , FUN_UP ];   # FunFam UniProt - i.e. rep (not currently used)
  rep_ID                <- dt_mgmt[ , REP_ID ];
  pdb_code              <- dt_mgmt[ , PDB_CODE ];
  pdb_chain             <- dt_mgmt[ , CHAIN_CODE ];
  num_swissprot_in_ff   <- dt_mgmt[ , NUM_SWISSPROT_IN_FF ];  # (Human) SwissProt members of FunFam
  mfc_mut_count         <- dt_mgmt[ , MFC_MUT_COUNT ];
  MIN_RAD               <- dt_mgmt[ , MIN_RAD ];
  MAX_RAD               <- dt_mgmt[ , MAX_RAD ];
  NO_TRIALS             <- dt_mgmt[ , NO_TRIALS ];

# Generate runID based on taskID & batchCode
  runID <- paste( "R", gsub( " ", "0", sprintf( "%4s", taskID ) ), arg_batchCode, sep = "" );

# Output textual information on this MutFamClust run
  info_out_file <- file.path( arg_out_dir, paste( arg_groupID, "_", runID, "_info.txt", sep = "" ) );
  info_out      <- paste( "MutClust info - ",
                          "groupID:", arg_groupID,
                          "; batch:", arg_batchCode,
                          "; task:", taskID,
                          "; sfam:", sf_ID , 
                          "; ff_num:", ff_num ,
                          "; PDB_code:", pdb_code, 
                          "; pdb_chain:", pdb_chain, 
                          "; rep_ID:", rep_ID, 
                          "; num_Swiss:", num_swissprot_in_ff,
                          "; mfc_mut_count:", mfc_mut_count,
                          "; MGMT_FILE:", arg_mgmt_file,
                          "; MUTDAT_FILE:", arg_mutdat_file,
                          "; PDB_DIR:", arg_pdb_dir, 
                          "; OUT_DIR:", arg_out_dir,
                          sep = ""
                          );
  write( info_out, info_out_file );

# Call MutClust run for this MutFam...
  dt_scoresALL <- mfc.singleMutFamClust(    arg_out_dir, 
                                            arg_groupID, 
                                            sf_ID, 
                                            ff_num,
                                            pdb_code, 
                                            pdb_chain, 
                                            rep_ID,
                                            num_swissprot_in_ff,
                                            mfc_mut_count, 
                                            MIN_RAD, 
                                            MAX_RAD, 
                                            NO_TRIALS, 
                                            arg_pdb_dir, 
                                            arg_mutdat_file,
                                            info_out
                                          );

  rdat_file <- file.path( arg_out_dir, paste0( arg_groupID, "_", runID, "_", MIN_RAD, "-", MAX_RAD, "_", NO_TRIALS, "_dtscoresALL.RData" ) );
  save( dt_scoresALL, file = rdat_file );
