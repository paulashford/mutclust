# structural_clustering/struct_clust_obs_muts.R
# Ash March 2015
# Counts mutations occuring on pdb structure within radius r for set of aminos_to_test

# 26/06/2018: Allow calculation of cluster size and passing in of dt_pdbres_mut to do so

# Oct'19: As have working git repo have clarified this and removed aminos_to_test param; function call uses 
# 			all obs mutants in mutdat file
#   Clusters: within r Angstoms calculates i) num_mut_res; ii) num_muts; iii) weighted_num_muts; iv) ff_weighted_num_muts
#	Significance by comparison with ranomisation trials using same cluster calculations (see struct_clust_randomise.R)

sc.clust_obs_muts <- function( pdb_xyz, r, mutated_pdb_residues, dt_pdbres_mut, num_swissprot_in_ff ){

	# Create a data.table to store the observed scores
	dt_scoresOBS <- data.table(	"trialNo" = 0, "index" = 0, "residue" = 0, "radius" = 0, "num_mut_res" = 0, 
								"MFC_MUT_COUNT_SUM" = 0, "WEIGHTED_MUT_SUM" = 0, "FF_WEIGHTED_MUT_SUM" = 0 );
	dt_scoresOBS <- dt_scoresOBS[ 0 ];

	writeLines( "==========================================================================================" );

	#print.default('Start: Cluster - observed mut data');
	writeLines( paste( "***** Experimentally observed muts - radius: ", r, " ******" ) );
	writeLines("Here are the mutated PDB residues: ");
	print.table( mutated_pdb_residues );

	for ( h in mutated_pdb_residues ){

		writeLines(paste("RESIDUE: ", h));

		cluster_res_rad <- sc.calc_struct_cluster_residues ( h = h, r = r, pdb = pdb_xyz ) ;
		writeLines( paste( "All residues within ",r,"A of residue ", h ) );
		print.table( cluster_res_rad );

		# Find intersection of cluster with known mutants (leave test residue in as this is how mutations are counted)
		cluster_intersect <- intersect( mutated_pdb_residues, cluster_res_rad );
		num_mut_res <- length( cluster_intersect );
		
		# Sum of mutations in cluster
		mut_count_sum 		<- sum ( dt_pdbres_mut[ PDB_RES_NAME %in% cluster_intersect, MFC_MUT_COUNT_SUM ] );

		# Sum of hill-weighted mutation counts
		weighted_mut_sum 	<- sum ( dt_pdbres_mut[ PDB_RES_NAME %in% cluster_intersect, WEIGHTED_MUT_SUM ] );

		# Sum of FF weighted...
		ff_weighted_mut_sum	<- sum ( dt_pdbres_mut[ PDB_RES_NAME %in% cluster_intersect, FF_WEIGHTED_MUT_SUM ] );
		
		writeLines( paste0( "Known mutants residues in cluster C(h=",h,",r=",r,"): ", paste( cluster_intersect, collapse = ", " ), "\n" ) );
		writeLines( paste0( "Summing muts from all ", num_mut_res, " mutated residues in cluster C(h=",h,",r=",r,"):", "\n",
							"-> mut_count sum = ", mut_count_sum, "\n",
							"-> Hill weighted (theta = ", HILL_THETA, ", exp = ", HILL_EXP,") mut_count sum = ", weighted_mut_sum, "\n",
							"-> FF weighted Hill sum ( human paralogues = ", num_swissprot_in_ff, ") = ", ff_weighted_mut_sum, "\n"	
						)
					);

		# PDB resno index (internal)
		res_index <- dt_pdbres_mut[ PDB_RES_NAME == h, index ];
		dt_scoresOBS <- rbind( dt_scoresOBS, as.list( c( -1, res_index, h , r, num_mut_res, mut_count_sum, weighted_mut_sum, ff_weighted_mut_sum  ) ) );
	
	}
	writeLines( "==========================================================================================" );

	# Return nice data.table with our scoring of experimentally observed mutants near the hotspots
	setkey( dt_scoresOBS, trialNo, radius, residue );
	#print.table(t(dt_scoresOBS));
	return( dt_scoresOBS );

}
