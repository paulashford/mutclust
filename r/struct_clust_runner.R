# structural_clustering/struct_clust_runner.R
# Ash 10/03/2015

# Test
# H1: For a given cut-off radius, low frequency (/"other") mutations are clustered around a particular 'hot-spot' mutation more than would be expected by chance.
# H0: The distribution of low freq  (/"other") mutations around a given 'hot-spot' mutation is no different to that expected by chance.

# 29/04/2015
# Make this into a function and allow passing of aminos_specified_list so can run for any bespoke
# list of residues, including whole domains.

# Jan/Feb 16: MutFamClust runs

# Oct'19: 	As have working git repo: mods to remove old FGFR and run_type="MutFamClust" (this was function call param 1)
# 			Also removed 'amino_specified...' type variables - this now just works on whatever mutations are in the mutdat file.
# 			Clusters: within r Angstoms calculates i) num_mut_res; ii) num_muts; iii) weighted_num_muts; iv) ff_weighted_num_muts
#			Significance by comparison with ranomisation trials using same cluster calculations (see struct_clust_randomise.R)

# Define sets
sc.struct_clust_runner <- function( dt_data_mutfreq, pdb_xyz, MAX_TRIALS, num_swissprot_in_ff, info_string, RADIUS ){
	# Oct '19 calls & vars simplified to remove redundant (confusing) naming from previous versions
 
	writeLines( paste( "START: sc.struct_clust_runner with MAX_TRIALS: ", MAX_TRIALS, "\n", info_string, "\n", sep = "" ) );

	# What is the set of all residues in the PDB?
	pdb_residues <- as.integer( unique( pdb_xyz$resno ) );

	# What are all mutated residues in the mutdat file for this dom rep pdb?
	mutated_residues <-  sort( unique( as.integer( dt_data_mutfreq$PDB_RES_NAME ) ) );

	# Which mutants overlap with PDB residues?
	# Given that mut_mapper etc already tests this, the mutdat file should only contain dom pdb residues,if using a standard FunFam rep.
	mutated_pdb_residues <- intersect( pdb_residues, mutated_residues );
	# This just ensures that mutfreq data only includes mutations that can be tested on this PDB
	dt_data_mutfreq <- dt_data_mutfreq[ PDB_RES_NAME %in% mutated_pdb_residues ];

	# For cluster calculations
		# Clusters: within r Angstoms calculates:
		# i) num_mut_res; 
		# ii) num_muts; 
		# iii) weighted_num_muts; 
		# iv) ff_weighted_num_muts
	
	# num_mut_res: how many mutated residues?
	num_mut_res <- length( mutated_pdb_residues );

	# num_muts: how many mutations in total summing across all obs mutated residues?
	num_muts 	<- sum( dt_data_mutfreq$MFC_MUT_COUNT );

	# Count mutations per residue (i.e. over all genes) 
	# SQL equivalent as test: /Users/ash/woofgit/mutclust/sql/mutclust_mut_count_by_residue.sql
	dt_data_mutfreq[ ,  MFC_MUT_COUNT_SUM := sum( MFC_MUT_COUNT_RESIDUE ), by = PDB_RES_NAME ];
	
	# Mutations by (unique) residue
	# Unique residue list in this dt then used for random sampling in the permutation trials
	dt_pdbres_mut	<- unique( pdb_xyz[ , .( resno ) ] );
	dt_pdbres_mut[ , index := .I ];
					
	dt_pdbres_mut <- merge(	dt_pdbres_mut,
							unique( dt_data_mutfreq[ , .( PDB_RES_NAME = as.numeric( PDB_RES_NAME ), MFC_MUT_COUNT_SUM ) ] ),
							by.x = c( "resno" ),
							by.y = c( "PDB_RES_NAME" ),
							all.x = TRUE
						);
	# Remove NAs
	dt_pdbres_mut[ is.na( MFC_MUT_COUNT_SUM ), MFC_MUT_COUNT_SUM := 0 ];

	# Rename merge column name
	setnames( dt_pdbres_mut, "resno", "PDB_RES_NAME" );

	# Weighted mut_count using Hill function to normalise
	dt_pdbres_mut[ , WEIGHTED_MUT_SUM 	:= mut_weight( MFC_MUT_COUNT_SUM, theta = HILL_THETA, exponent = HILL_EXP ) ];
	# The weighted mut_count 'normalised' by num. human FunFam members
	dt_pdbres_mut[ , FF_WEIGHTED_MUT_SUM := WEIGHTED_MUT_SUM / num_swissprot_in_ff ];

	# Overall probability of a mutant residue, based on obervations of actual mutants
	# Oct '19 - using explicit re-sampling of residue positions now
	##p_mut <- num_mut_res / length( pdb_residues );

	# Initialise data.table to hold all exp obs muts and randomisation trials for different radii from given residues.
	dt_scores <- data.table( 	"trialNo" = 0, "index" = 0, "residue" = 0, "radius" = 0, "num_mut_res" = 0, 
								"MFC_MUT_COUNT_SUM" = 0, "WEIGHTED_MUT_SUM" = 0, "FF_WEIGHTED_MUT_SUM" = 0 );
	dt_scores <- dt_scores[ 0 ];

	# Loop radii for testing (single radii runs now!)
	#for (r in MIN_RAD : MAX_RAD){
	r <- RADIUS; 
		# OBSERVED MUTATED RESIDUES AND LOCAL MUTATION COUNT IN SPHERICAL CLUSTERS
	writeLines( " " );
	writeLines( "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" );
	writeLines( paste0( "@CALL: sc.clust_obs_muts using radius: ", r, 
						"-> num PDB residues: ", length( pdb_residues ), "\n",
						"-> num mutated residues (mutdat): ", length( mutated_residues ), "\n", 
						"-> num mutated PDB residues: ", num_mut_res, "\n",
						"-> num mutations (total): ", num_muts, "\n",
						"All PDB residues: ", paste( pdb_residues, collapse = "," ), "\n",
						"Mutated PDB residues: ", paste( mutated_pdb_residues, collapse = "," ), "\n"
					)
				);
	writeLines( "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" );
	dt_expScores <- sc.clust_obs_muts( pdb_xyz, r, mutated_pdb_residues, dt_pdbres_mut, num_swissprot_in_ff ) ;

	# TEST CLUSTERING OF RANDOMISED MUTATIONS AROUND EACH MUTATED RESIDUE
	writeLines( " " );
	writeLines( "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" );
	writeLines(	paste0("@CALL: sc.run_random_mut_trials using radius: ",r, "\n",
						"-> MAX_TRIALS: ", MAX_TRIALS, "\n",
						"-> num_swissprot_in_ff: ", num_swissprot_in_ff, "\n"
					) 
				);
	writeLines( "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" );
	dt_randScores  <- sc.run_random_mut_trials( pdb_xyz, r, mutated_pdb_residues, dt_pdbres_mut, MAX_TRIALS, num_swissprot_in_ff )
	
	dt_scores <- rbind( dt_scores, rbind ( dt_expScores, dt_randScores ) );

	writeLines( paste( "FINISH: sc.struct_clust_runner with MAX_TRIALS: ", MAX_TRIALS, "\n", info_string, "\n", sep = "" ) );

	# Tidy vars up and return full data.table of exp./obs. and randomised trials
	rm( dt_randScores, dt_expScores, r );
	setkey( dt_scores, trialNo, radius, residue, index );
	return( dt_scores  )

}
