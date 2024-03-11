# MutClust_analysis_functions.R
# 23/10/2019 
# Renamed from: structural_clustering/analysis_utils.R
# Removed old sections (e.g. iGraph related)

# Changes reflect 4 calculations now returned by MutClust (e.g. Hill normalised)
# See /Users/ash/woofgit/mutclust/r/struct_clust_randomise.R
	
	# Get means and standard deviations for 4 metrics based on randomisation trials for input radius of clusters
	# Includes p-vals and correction over all p-vals
	sca.getMeanDT <- function( dt_scoresALL, this_radius, p_correct_method = "BH", p_corr_cutoff = 0.05 ){

		setkey( dt_scoresALL, residue );

		# Means and standard deviations of the random trials (trialNo ==-1 are actual observations)
		dt_meansd 	<-  dt_scoresALL[ ( trialNo > 0 & radius == this_radius ),
											.( 	num_mut_res_mean = round( mean( num_mut_res, na.rm = TRUE ), 3 ), 
												num_mut_res_sd   = round( sd( num_mut_res, na.rm = TRUE ), 3 ),						
												mut_count_sum_mean = round( mean( MFC_MUT_COUNT_SUM, na.rm = TRUE ), 3 ), 
												mut_count_sum_sd   = round( sd( MFC_MUT_COUNT_SUM, na.rm = TRUE ), 3 ),
												weighted_mut_sum_mean = round( mean( WEIGHTED_MUT_SUM, na.rm = TRUE ), 6 ), 
												weighted_mut_sum_sd = round( sd( WEIGHTED_MUT_SUM, na.rm = TRUE ), 6 ),
												ff_weighted_mut_sum_mean = round( mean( FF_WEIGHTED_MUT_SUM, na.rm = TRUE ), 8 ), 
												ff_weighted_mut_sum_sd = round( sd( FF_WEIGHTED_MUT_SUM, na.rm = TRUE ), 8 ) 
											), 
											by = c( "residue" ) 
									];

		# Join mean/sd to scores table (so we have obs mut counts too!)
		setkey( dt_meansd, residue );
		dt_meansd <-  dt_scoresALL[ trialNo == -1 ][ dt_meansd ];

		# Use pnorm, assuming normal distributiton of randomisations on MutClusts
		# num_mut_res
		dt_meansd[ , num_mut_res_p 		:= pnorm( num_mut_res, num_mut_res_mean, num_mut_res_sd, lower.tail = FALSE ) ];
		# Adjustment over p-values for testing multiple atoms in same PDB in many counts
		dt_meansd[ , num_mut_res_p_corr := p.adjust( dt_meansd[ , num_mut_res_p ], method = p_correct_method )  ];
		# Simple significance flag 
		dt_meansd[ , num_mut_res_p_corr_sig := "N" ];
		dt_meansd[ num_mut_res_p_corr < p_corr_cutoff, num_mut_res_p_corr_sig := "Y" ];

		# mut_count_sum
		dt_meansd[ , mut_count_sum_p := pnorm( MFC_MUT_COUNT_SUM, mut_count_sum_mean, mut_count_sum_sd, lower.tail = FALSE ) ];
		# Adjustment over p-values for testing multiple atoms in same PDB in many counts
		dt_meansd[ , mut_count_sum_p_corr := p.adjust( dt_meansd[ , mut_count_sum_p ], method = p_correct_method )  ];
		# Simple significance flag 
		dt_meansd[ , mut_count_sum_p_corr_sig := "N" ];
		dt_meansd[ mut_count_sum_p_corr < p_corr_cutoff, mut_count_sum_p_corr_sig := "Y" ];

		# weighted_mut_sum
		dt_meansd[ , weighted_mut_sum_p := pnorm( WEIGHTED_MUT_SUM, weighted_mut_sum_mean, weighted_mut_sum_sd, lower.tail = FALSE ) ];
		# Adjustment over p-values for testing multiple atoms in same PDB in many counts
		dt_meansd[ , weighted_mut_sum_p_corr := p.adjust( dt_meansd[ , weighted_mut_sum_p ], method = p_correct_method )  ];
		# Simple significance flag 
		dt_meansd[ , weighted_mut_sum_p_corr_sig := "N" ];
		dt_meansd[ weighted_mut_sum_p_corr < p_corr_cutoff, weighted_mut_sum_p_corr_sig := "Y" ];

		# ff_weighted_mut_sum
		dt_meansd[ , ff_weighted_mut_sum_p := pnorm( FF_WEIGHTED_MUT_SUM, ff_weighted_mut_sum_mean, ff_weighted_mut_sum_sd, lower.tail = FALSE ) ];
		# Adjustment over p-values for testing multiple atoms in same PDB in many counts
		dt_meansd[ , ff_weighted_mut_sum_p_corr := p.adjust( dt_meansd[ , ff_weighted_mut_sum_p ], method = p_correct_method )  ];
		# Simple significance flag 
		dt_meansd[ , ff_weighted_mut_sum_p_corr_sig := "N" ];
		dt_meansd[ ff_weighted_mut_sum_p_corr < p_corr_cutoff, ff_weighted_mut_sum_p_corr_sig := "Y" ];

		
		# Include correction method
		dt_meansd[ , p_info := paste0( "pnorm_", p_correct_method,"_cutoff_sig_", p_corr_cutoff ) ];

		return( dt_meansd );
	}

	# Get high outliers by type for given significance level
	sca.getHighOutliers <- function( dt_meansd, sigmaz, type )
	{
		# SUPERCEDED by explicit pnorm() calculations in sca.getMeanDT 

		# if ( type == "num_mut_res" ){
		# 	dt_meansd[ trialNo == -1, sig_num_mut_res := 'N' ] ;
		# 	dt_meansd[ trialNo == -1 & num_mut_res > ( num_mut_res_mean + ( sigmaz * num_mut_res_sd ) ), sig_num_mut_res := 'Y' ];
		# 	#return( meansd_randDT[ trialNo == -1 ][ meansd_randDT ][ num_mut_res > ( num_mut_res_mean + ( sigmaz * num_mut_res_sd ) ) ] );
		# }
		# if ( type == "mut_count_sum" ){
		# 	dt_meansd[ trialNo == -1, sig_mut_count_sum := 'N' ] ;
		# 	dt_meansd[ trialNo == -1 & MFC_MUT_COUNT_SUM > ( mut_count_sum_mean + ( sigmaz * mut_count_sum_sd ) ), sig_mut_count_sum := 'Y' ];
		# 	#return( meansd_randDT[ trialNo == -1 ][ meansd_randDT ][ MFC_MUT_COUNT_SUM > ( mut_count_sum_mean + sigmaz * mut_count_sum_sd ) ] );
		# }
		# if ( type == "weighted_mut_sum" ){
		# 	dt_meansd[ trialNo == -1, sig_weighted_mut_sum := 'N' ] ;
		# 	dt_meansd[ trialNo == -1 & WEIGHTED_MUT_SUM > ( weighted_mut_sum_mean + ( sigmaz * weighted_mut_sum_sd ) ), sig_weighted_mut_sum := 'Y' ];
		# 	#return( meansd_randDT[ trialNo == -1 ][ meansd_randDT ][ WEIGHTED_MUT_SUM > ( weighted_mut_sum_mean + sigmaz * weighted_mut_sum_sd ) ] );
		# }

		# if ( type == "ff_weighted_mut_sum" ){
		# 	dt_meansd[ trialNo == -1, sig_ff_weighted_mut_sum := 'N' ] ;
		# 	dt_meansd[ trialNo == -1 & FF_WEIGHTED_MUT_SUM > ( ff_weighted_mut_sum_mean + ( sigmaz * ff_weighted_mut_sum_sd ) ), sig_ff_weighted_mut_sum := 'Y' ];
		# 	#return( meansd_randDT[ trialNo == -1 ][ meansd_randDT ][ FF_WEIGHTED_MUT_SUM > ( ff_weighted_mut_sum_mean + sigmaz * ff_weighted_mut_sum_sd ) ] );
		# }
		
		#return( paste0( "Can't find type ", type,". Maybe you should have a beer instead." ) );
		#return( dt_meansd )
		return( "SUPERCEDED by explicit pnorm() calculations in sca.getMeanDT(...)" )  
	}

	# Handy object existence checker - works without quotes to see if variable defined without
	# 1st just throwing Error: object '....' not found!
	# From https://stackoverflow.com/questions/9368900/how-to-check-if-object-variable-is-defined-in-r
	is.defined <- function( sym ) 
	{
		sym <- deparse( substitute( sym ) )
		env <- parent.frame( )
		exists( sym, env )
	}



