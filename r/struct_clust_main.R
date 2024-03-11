# Struct cluster main
# Ash 10/03/2015

# 27/05/2015
# Update to allow command line argument
# [1] DATA_DIR : The root dir for the PDB file
# Some lab notes...
# https://www.evernote.com/shard/s8/nl/959078/330c33c5-02e7-4b60-b383-2c2af5e75d83/

# Ash: 07/12/2015
# Updates re: MutFamClust - e.g. calling new views, splitting parameters & runs
# of struct_clust and MutFamClust.

# Ash: 28/01/2016
# Adding method for scripted calling (i.e. large runs in pipeline)

# 04/02/16: See MFC_exampe_MGMT.tsv for run management format.

# Oct'19 - updates for explicit input args & clarity

# -----------------------------------------
# MutFamClust main call (single MutFam)
# Expected calling fn: MutFamClust_group.R
# -----------------------------------------
mfc.singleMutFamClust <- function(  out_dir, 
                                    groupID,
                                    sf_ID, 
                                    ff_num, 
                                    PDB_code, 
                                    PDB_chain, 
                                    rep_ID,
                                    num_swissprot_in_ff,
                                    mfc_mut_count,
                                    MIN_RAD, 
                                    MAX_RAD, 
                                    MAX_TRIALS, 
                                    pdb_dir, 
                                    mutdat_file,
                                    info_string ){
  # parameters
  # out_dir  : Directoy for this group's MutFamClust runs
  # groupID    : groupID for MutFamClust this group (of MutFamClust runs) 
  # MIN/MAX_RAD : Min and max radii tested
  # MAX_TRIALS : Number of randomisation runs

    # Load libraries
    require("data.table");
    require(bio3d);
    # Using rdist to find distances between sets of atoms
    require(fields);
    #require(ggplot2);
    require(stringr);

    # GLOBAL VARS
    # For Hill-function 'compression' of mutation hotspots
    HILL_THETA  <<- 2;
    HILL_EXP    <<- 3;

    VERBOSE_OUTPUT  <<- FALSE;

    # Source required files
    source( "struct_clust_util.R" );
    source( "struct_clust_load_data.R" );
    source( "struct_clust_obs_muts.R" );
    source( "struct_clust_randomise.R" );
    source( "struct_clust_runner.R" );

    dt_data_mutfreq <- sc.load_MutFamClust_freq_tsv( mutdat_file, sf_ID, ff_num, PDB_code, PDB_chain );

    if ( length( unique( dt_data_mutfreq$PDB_CODE ) ) !=1 ){
        #print("More than 1 PDB id returned")
        return( -11 );
      } else {
        pdb_xyz <- sc.load_pdb( file.path( pdb_dir, rep_ID ) );
    }

    #info_string <- paste0( groupID, "|", sf_ID, "|", ff_num, "|", rep_ID, "|", pdb_dir, "|", mutdat_file, "|", MIN_RAD, "|", MAX_RAD, "|", MAX_TRIALS );
    dt_scoresALL <- sc.struct_clust_runner( dt_data_mutfreq, pdb_xyz, MAX_TRIALS, num_swissprot_in_ff, info_string, RADIUS = MIN_RAD );
    
    return( dt_scoresALL )
    
}
