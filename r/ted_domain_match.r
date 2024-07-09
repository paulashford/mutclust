# mutclust AF TED100 domain matching
# 14 06 2024

library(tidyverse)
source('/Users/ash/git/mutclust/r/struct_clust_load_data.R')

# TED domain matching 
# TED dom file
ted100_file		<- '/Users/ash/data/funvar_pipeline/neofun_ted_doms/ash_ted100_hits'
df_ted			<- read_delim(ted100_file)
# MutClust pre MGMT file
# mcpre_file		<- '/Users/ash/Dropbox/bioinf/neofun/run/mutclust_af/tst_mcaf02/mfc_mgmt_swisscount_missense_af_pre_20240614.tsv'
mcpre_file		<- '/Users/ash/Dropbox/bioinf/neofun/run/mutclust_af/dat/mfc_mgmt_swisscount_missense_af_pre_20240614.tsv'
df_mcpre		<- read_delim(mcpre_file)

# MutClust FunFams Unique domains
df_mc_counts1 	<- df_mcpre %>%
					select(UNIPROT_ACC, SUPERFAMILY_ID, FUNFAM_NUMBER,  MFC_MUT_COUNT) %>%
					group_by(UNIPROT_ACC, SUPERFAMILY_ID, FUNFAM_NUMBER) %>%
					mutate(max_mfc = max(MFC_MUT_COUNT) ) %>%
					select(-MFC_MUT_COUNT) %>%
					# arrange(SUPERFAMILY_ID, FUNFAM_NUMBER, UNIPROT_ACC) %>%
					unique() %>%
					arrange(desc(max_mfc))

# MERGE MC AND TED
# separate ted_id
df_ted	<- df_ted %>%
	separate(ted_id, sep = "-", into = c("AF", "uniprot", "F", "model"), remove = FALSE )

# A0A061IFB9_17_310
# join on uniprot initially
df_merge	<- left_join(df_mcpre,	
							df_ted,
							by = c("UNIPROT_ACC" = "uniprot"),
							keep = TRUE,
							na_matches = "never",
							multiple = "all",
							relationship = "many-to-many"
						)

# Check if in TED100 or not prior to filters...
df_missing_or_not <- df_merge %>% 
	# filter(is.na(ted_id) ) %>%
	filter(!is.na(ted_id) ) %>%
	select(UNIPROT_ACC, SUPERFAMILY_ID, FUNFAM_NUMBER,  MFC_MUT_COUNT) %>%
	group_by(UNIPROT_ACC, SUPERFAMILY_ID, FUNFAM_NUMBER) %>%
	mutate(max_mfc = max(MFC_MUT_COUNT)) %>%
	select(-MFC_MUT_COUNT) %>%
	unique() %>%
	arrange(desc(max_mfc))
View(df_missing_or_not)

# check ranges
df_merge1	<- df_merge %>%
					filter(num_segments == 1) %>%
					separate(UP_AA_RANGE, sep = "-", 
							into = c("UP_AA_RANGE_LOW", "UP_AA_RANGE_HIGH"), 
							remove = FALSE ) %>%
					separate(chopping, sep = "-", 
							into = c("chopping_low", "chopping_high"), 
							remove = FALSE )
					# mutate(match_res_count = )
						# group_by(UNIPROT_ACC)
					
# calc amino acid overlap in ranges and what proportion of our FunFam rep that represents
df_merge1 	<- df_merge1 %>%
				mutate( chopping_low = as.numeric(chopping_low) ) %>%
				mutate( chopping_high = as.numeric(chopping_high) ) %>%
				mutate( UP_AA_RANGE_LOW = as.numeric(UP_AA_RANGE_LOW) ) %>%
				mutate( UP_AA_RANGE_HIGH = as.numeric(UP_AA_RANGE_HIGH) ) %>%
				rowwise() %>%
				mutate( aa_overlap = sum(between( UP_AA_RANGE_LOW:UP_AA_RANGE_HIGH, chopping_low, chopping_high ) ) ) %>%
				mutate( aa_prop_overlap = aa_overlap / (UP_AA_RANGE_HIGH - UP_AA_RANGE_LOW + 1) )	
View(df_merge1)

df_merge_filt	<- df_merge1 %>%
						group_by(UNIPROT_ACC) %>%
						slice_max(order_by = aa_overlap, n = 1, with_ties = FALSE) %>%
						filter(aa_overlap >0) %>%
						filter(aa_prop_overlap >0.5 )
View(df_merge_filt)
# pLDDT >=90 only (1,028 FunFam rep domains -> TED domains)
df_merge_filt_90	<- filter(df_merge_filt, plddt >=90)
df_merge_filt_80	<- filter(df_merge_filt, plddt >=80)
	
View(df_merge_filt_90)

# df_filt_counts <- df_merge_filt %>% 
df_filt_counts_90 <- df_merge_filt_90 %>% 
	select(UNIPROT_ACC, SUPERFAMILY_ID, FUNFAM_NUMBER,  MFC_MUT_COUNT) %>%
	group_by(UNIPROT_ACC, SUPERFAMILY_ID, FUNFAM_NUMBER) %>%
	mutate(max_mfc = max(MFC_MUT_COUNT)) %>%
	select(-MFC_MUT_COUNT) %>%
	unique() %>%
	arrange(desc(max_mfc))
View(df_filt_counts_90)

# write doms to get doms!  
write_delim( select( ungroup(df_merge_filt), ted_id), 
			file = "/Users/ash/Dropbox/bioinf/neofun/run/mutclust_af/tst_mcaf02/ted100_match_2270"
			)
# ** These TEDs extracted by Nico: mutclust_af/ash_models_2770 **

# Create MutClust MGMT for
# TED100, match FunFamv4.2 rep with amino seq overlap >0.5; pLDDT >=90
# colnames(df_mcpre)
#  [1] "TASKID_ROWNUM"       "TASKID"              "SUPERFAMILY_ID"     
#  [4] "FUNFAM_NUMBER"       "EF"                  "GENE"               
#  [7] "FUN_UP"              "REP_ID"              "UNIPROT_ACC"        
# [10] "UP_AA_RANGE"         "PDB_CODE"            "CHAIN_CODE"         
# [13] "NUM_SWISSPROT_IN_FF" "MFC_MUT_COUNT"       "MIN_RAD"            
# [16] "MAX_RAD"             "NO_TRIALS"  
# Initial sanity-check
df_mgmt_out <- df_merge_filt_90 %>%
	select( TASKID_ROWNUM, TASKID, 
			SUPERFAMILY_ID, FUNFAM_NUMBER, 
			EF, GENE, FUN_UP, REP_ID, UNIPROT_ACC, 
			UP_AA_RANGE, PDB_CODE, CHAIN_CODE, NUM_SWISSPROT_IN_FF, MFC_MUT_COUNT, 
			MIN_RAD, MAX_RAD, NO_TRIALS, 
			ted_id, chopping)
# now remove grouping (UNIPROT_ACC) and any extraneous cols...
df_mgmt_out1 <- df_mgmt_out %>%
	ungroup() %>%
	select( TASKID, SUPERFAMILY_ID, FUNFAM_NUMBER, EF, GENE, FUN_UP, 
			REP_ID = ted_id, 
			PDB_CODE, CHAIN_CODE, 
			NUM_SWISSPROT_IN_FF, MFC_MUT_COUNT, 
			MIN_RAD, MAX_RAD, NO_TRIALS
			) %>%
	mutate( REP_ID = paste0(REP_ID, ".pdb") ) 
# Add sequential TASKID for the job scheduler
n_tasks 		<- nrow(df_mgmt_out1)
df_mgmt_out1	<- mutate(df_mgmt_out1, TASKID = 1:n_tasks )

# write mgmt!
mgmt_file_name 	<- '/Users/ash/Dropbox/bioinf/neofun/run/mutclust_af/tst_mcaf02/mfc_mgmt_swisscount_missense_af_20240614.tsv'
write_delim( df_mgmt_out1, 
			file = mgmt_file_name,
			delim = "\t",
			quote = "none",
			na = "NA"
			)
