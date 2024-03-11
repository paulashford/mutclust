# structural_clustering/struct_clust_randomise.R
# Ash March 2015

# Oct '19 - Modified random trials to use mutdat file as a parameter, rather than just an overall 
# probability of mutation (p_mut)
#	1) num_mut_res	: How many mutated residues nearby?
#	2) MFC_MUT_COUNT_SUM is summed mutations for residue over all genes
#	3) WEIGHTED_MUT_SUM applies a Hill curve normalisation to 'compress' (with theta=2, exp=3 means ~ [0,6])
# 	4) FF_WEIGHTED_MUT_SUM divides the weighted sum by number of swissprot human paralogues in FunFam to
#		account for bias from FunFams with many members
#	A random trial samples from *all* PDB residues, with the new ordering of residues 'acquiring' mutations (or lack of) 
# 	at the original residue index. Essentially, the mutation count, keyed by residue number, is now keyed by the new randomly sampled
# 	residue number.  
#   Clusters: within r Angstoms calculates i) num_mut_res; ii) num_muts; iii) weighted_num_muts; iv) ff_weighted_num_muts

# 14 Nov '19 - further speed improvements by properly using data.table keys
# See: /Users/ash/woofgit/mutclust/r/MutClust_dev_speed_test_rand_in_vs_key.R

# Previous random model (see git for code!)
		# 26/06/2018: Allow passing of mut_sum (how many muts on all mut residues in PDB?) and calculation of cluster size
		# Also supress log file lines

		# Performs randomised mutation of residues (simple - 0: Not, 1: Mutated)
		# Creates a data.table recording results of randomisation trials.
		# For each residue in list residues_of_interest
		# 	A trial,t, is:
		#		1) Randomly mutate (sampling from Uniform distribution) each amino position in pdb.
		# 		2) Map these mutants back to the structure pdb.
		# 		3) Count how many random mutants are within rA of residue
		#   	4) Update randScore[[h]] array with mutant count for trial t
		# https://www.evernote.com/shard/s8/nl/959078/e5e98dcf-6797-4ab2-9ce9-7b7cf02ae048/

sc.run_random_mut_trials <- function( pdb_xyz, r, mutated_pdb_residues, dt_pdbres_mut, no_trials, num_swissprot_in_ff ) {

	# Create a data.table to store the randomisation scores
	dt_scoresRAND <- data.table(	"trialNo" = 0, "index" = 0, "residue" = 0, "radius" = 0, "num_mut_res" = 0, 
								"MFC_MUT_COUNT_SUM" = 0, "WEIGHTED_MUT_SUM" = 0, "FF_WEIGHTED_MUT_SUM" = 0 );
	dt_scoresRAND <- dt_scoresRAND[ 0 ];

			# The test or map (I think...)
			#dt_resno_idx[resno %in% mutated_pdb_residues & sampled %in% mutated_pdb_residues]

	writeLines( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" );
	writeLines( paste0( "(~) ***** Permutation trials - radius: ", r, ", num_trials: ", no_trials, "num_swissprot_in_ff: ", num_swissprot_in_ff, " ******" ) );

	start_time <- Sys.time();
	writeLines( paste0( "start_time: ", start_time ) );

	# Loop random trials
	for ( rand_trial in 1:no_trials ){
	
		writeLines( paste0( "(~) -> Trial: ", rand_trial ) );
		
		# Randomly sample the PDB_RES and place into new sample_ column
		dt_pdbres_mut[ , sample_PDB_RES_NAME := sample( PDB_RES_NAME ) ];
		setkey( dt_pdbres_mut, sample_PDB_RES_NAME );

		# What are mutated residues *in this randomised sample*?
		rand_mutated_residues <- dt_pdbres_mut[ MFC_MUT_COUNT_SUM > 0, sample_PDB_RES_NAME ];

		# Loop mutated residues *as observed in mutdat*
		for ( this_res in mutated_pdb_residues ){
		
			# Get cluster around this residue for given radius
			cluster_res_rad <- sc.calc_struct_cluster_residues( h = this_res, r = r, pdb = pdb_xyz ) ;
			if ( VERBOSE_OUTPUT ){
				writeLines( paste( "(~) -> All residues within ",r,"A of residue ", this_res ) );
				print.table( cluster_res_rad );
			}

			# Find intersection of cluster with randomised mutantss (leave test residue in as this is how mutations are counted)
			rand_cluster_intersect <- intersect( rand_mutated_residues, cluster_res_rad );
			
			# Analogous calculations to observed, but using the sampled residue numbering...
			rand_num_mut_res 			<- length( rand_cluster_intersect );

			# Sum of mutations in cluster
			#rand_mut_count_sum 		<- dt_pdbres_mut[ sample_PDB_RES_NAME %in% rand_cluster_intersect, sum( MFC_MUT_COUNT_SUM ) ];
			rand_mut_count_sum 			<- dt_pdbres_mut[ .( rand_cluster_intersect ), sum( MFC_MUT_COUNT_SUM ) ];

			# Sum of hill-weighted mutation counts
			#rand_weighted_mut_sum 		<- dt_pdbres_mut[ sample_PDB_RES_NAME %in% rand_cluster_intersect, sum( WEIGHTED_MUT_SUM ) ] ;
			rand_weighted_mut_sum 		<- dt_pdbres_mut[ .( rand_cluster_intersect ), sum( WEIGHTED_MUT_SUM ) ] ;

			# Sum of FF weighted...
			#rand_ff_weighted_mut_sum	<- dt_pdbres_mut[ sample_PDB_RES_NAME %in% rand_cluster_intersect, sum( FF_WEIGHTED_MUT_SUM ) ];
			rand_ff_weighted_mut_sum	<- dt_pdbres_mut[ .( rand_cluster_intersect ), sum( FF_WEIGHTED_MUT_SUM ) ];

			if ( VERBOSE_OUTPUT ){
				writeLines( paste0( "(~) Randomised mutant residues in cluster C(h=",this_res,", r=",r,"): ", paste( rand_cluster_intersect, collapse = ", " ), "\n" ) );
				writeLines( paste0( "(~) Summing muts from all ", rand_num_mut_res, " mutated randomised residues in cluster C(h=",this_res,", r=",r,"):", "\n",
									"(~) -> rand mut_count sum = ", rand_mut_count_sum, "\n",
									"(~) -> rand Hill weighted (theta = ", HILL_THETA, ", exp = ", HILL_EXP,") mut_count sum = ", rand_weighted_mut_sum, "\n",
									"(~) -> rand FF weighted Hill sum ( human paralogues = ", num_swissprot_in_ff, ") = ", rand_ff_weighted_mut_sum, "\n"	
							)
						);
			}
			# Update the randomised scores table
			# PDB resno index (internal)
			#rand_res_index <- dt_pdbres_mut[ sample_PDB_RES_NAME == this_res, index ];
			
			# Note 26/11/19 - this is causing crash on occasion if PDB has negative numbered residues
			#rand_res_index <- dt_pdbres_mut[ this_res, index ];
			rand_res_index <- dt_pdbres_mut[ .(this_res), index ];

			dt_scoresRAND <- rbind( dt_scoresRAND, as.list( c( rand_trial, rand_res_index, this_res, r, rand_num_mut_res, rand_mut_count_sum, rand_weighted_mut_sum, rand_ff_weighted_mut_sum ) ) );
			
		}
		# [Remove excess I/O]
	 	##writeLines( paste0( "(~) ... end of trial: ", rand_trial ) );
	}
	
	setkey(dt_scoresRAND, trialNo, radius, residue );
	
	end_time <- Sys.time();
	writeLines( paste0( "end_time: ", end_time ) );
	writeLines( paste0( "end_time - start_time: ", end_time - start_time )  );

	return( dt_scoresRAND );

}
