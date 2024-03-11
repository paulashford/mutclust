# structural_clustering/struct_clust_load_data.R
# Ash 10/03/2015

# 28/05/2015
# Made a more generic load function & now uses Oracle FGFR db connection for frequency reps.

# 07/12/2015
# Amendments for MutFamClust - sc.load_MutFamClust_freq

# 25/01/2016
# MutFamClust for cs_cluster (MutFamClust^2?!)
# Allow tsv-driven mutation data due to issues with tunnel from cluster to
# gene 3d01.
# New method: sc.load_MutFamClust_freq_tsv

# *************
# Load data
# *************

sc.bioconductor <-function(){

  source("http://bioconductor.org/biocLite.R")
  biocLite()
  library(Biostrings)
}

# ------------------------------
# Get PDB structure...
# ------------------------------
sc.load_pdb <-function( pdb_file ){
 # rm(pdb_xyz);
  # if (missing(data.dir)==T){
  #   pdb_file <- pdb_file;
  # } else{
  #   pdb_file <- paste(data.dir, pdb_file, sep="/");
  # }
  # if (data.dir==""){
  #   pdb_file <- pdb_file;
  # }else{
  #   #pdb_file <- paste(data.dir, pdb_file, sep="/");
  #   pdb_file <- file.path(data.dir, pdb_file);
  # }

  pdb <- read.pdb(pdb_file);

  # Don't pick up HETATOMs!
  pdb_xyz <-data.table(pdb$atom[pdb$atom$type=="ATOM",c("resno", "resid", "elety","x","y","z")]);
  setkey(pdb_xyz,resno);
  rm(pdb);
  return(pdb_xyz);
}

# ------------------------------
# FGFR3 specific amino mutant freqs with alignment to F3active/apo models
# ------------------------------
sc.load_fgfr3_spec_amino_freq <- function (dbConn){
    data_combo <-data.table(dbGetQuery(dbConn,"SELECT * FROM FGFR.VW_MUT_FREQ_6_F3_SPEC_AM_NOTOT WHERE FGFR3_AA>396"));
    setkey(data_combo,"GLOBAL_AL");
    data_combo <- sc.add_manual_combo_mut_changes(data_combo);
    return(data_combo);
}

# ------------------------------
# Family wide COMBO: Global Align/Cancer muts/Inher muts
# ------------------------------
sc.load_fgfr_combo_muts <- function (dbConn)
# 23/07/15 VW_CANCER_AND_INH_MUTS_GAL updated to use alignment v3
# includes FGFR3 models of active/inactive
  {

    # Note: This only returns germline mutations where there is also a cancer mutation in FGFRs1-4
    data_combo <-data.table(dbGetQuery(dbConn,"SELECT * FROM FGFR.VW_CANCER_AND_INH_MUTS_GAL"));
    setkey(data_combo,"GLOBAL_AL");
    data_combo <- sc.add_manual_combo_mut_changes(data_combo);
    return(data_combo);
  }

sc.add_manual_combo_mut_changes <- function(data_combo)
# At present, it's easier to add odd bespoke mutation info here than change the Oracle db
{
  data_combo[data_combo[FGFR3_AA==647],"TOT_CANCER_MUTS"]=1;
  data_combo[data_combo[FGFR3_AA==647],"FGFR3_CANCER"]=1;
  return(data_combo);
}


# ------------------------------
# Family wide mutation frequencies
# ------------------------------
sc.load_fgfr_fam_mut_freq <- function (dbConn)
  {
   # data_mutfreq_all <-data.table(dbGetQuery(dbConn,"SELECT * FROM FGFR.VW_MUT_FREQ_REP_GAL_1_OV_TOTS"));
   # 01/07/2015 - v6 of freq rep includes FGFR4 COSMIC, changes to BioMuta, and
    data_mutfreq_all <-data.table(dbGetQuery(dbConn,"SELECT * FROM FGFR.VW_MUT_FREQ_REP_6"));
    setkey(data_mutfreq_all,"GLOBAL_AL");
    data_mutfreq_all <- sc.add_manual_mut_changes(data_mutfreq_all);
    return(data_mutfreq_all);
  }

  # -------------------------------------
  # Load MutFamClust mutation frequencies (Oracle)
  # -------------------------------------
  sc.load_MutFamClust_freq <- function (dbConn, SUPERFAMILY_ID, FUNFAM_NUMBER, PDB_CODE, PDB_CHAIN)
    {
      # Ash 07/12/2015 load data for MutFamClust
      sql <- paste("SELECT * FROM FGFR.VW_MF_MUTFAMCLUST_PANCAN_PDB
                    WHERE SUPERFAMILY_ID='", SUPERFAMILY_ID,"' AND FUNFAM_NUMBER='", FUNFAM_NUMBER ,"' AND PDB_CODE='", PDB_CODE, "' AND CHAIN_CODE='", CHAIN_CODE, "'",
                    sep="");

      load_MutFamClust_freq <- data.table(dbGetQuery(dbConn,sql));
      #setkey(data_mutfreq_all,"GLOBAL_AL");
      return(load_MutFamClust_freq);
    }

# -------------------------------------
# Load MutFamClust mutation frequencies (tsv)
# -------------------------------------
# 25/01/2016 Load MutFam freqs from flat file (rather than Oracle db if not av.)
sc.load_MutFamClust_freq_tsv <- function( tsv_file, superfamily_id, funfam_number, PDB_code, PDB_chain )
    {
      require(data.table);

      mutFamClust_freq <- data.table( 
                            read.table( tsv_file, 
                                        header = TRUE, 
                                        sep = "\t", 
                                        quote = "\"",
                                        dec = ".", 
                                        stringsAsFactors = FALSE,
                                        as.is = TRUE,
                                        colClasses=c( "PDB_CODE" = "character" ) 
                                      ) 
                                    );
      # 11/06/2020 If no PDB_chain (e.g. when using model structure) then don't include in filter
      if ( PDB_chain == '-' ){
        #mutFamClust_freq <- data.table(read.csv(file=tsv_file, sep="\t", stringsAsFactors=FALSE, as.is=TRUE, header=TRUE, colClasses=c("PDB_CODE" = "character")));
        mutFamClust_filter <- mutFamClust_freq[ SUPERFAMILY_ID == superfamily_id & 
                                                FUNFAM_NUMBER == funfam_number & 
                                                PDB_CODE == PDB_code 
                                              ];

      } else{
        #mutFamClust_freq <- data.table(read.csv(file=tsv_file, sep="\t", stringsAsFactors=FALSE, as.is=TRUE, header=TRUE, colClasses=c("PDB_CODE" = "character")));
        mutFamClust_filter <- mutFamClust_freq[ SUPERFAMILY_ID == superfamily_id & 
                                              FUNFAM_NUMBER == funfam_number & 
                                              PDB_CODE == PDB_code & 
                                              PDB_chain == PDB_chain ];
      }
      
      return(mutFamClust_filter);
    }

sc.add_manual_mut_changes <- function(data_mutfreq_all)
# At present, it's easier to add odd bespoke mutation info here than change the Oracle db
{
  data_mutfreq_all[data_mutfreq_all[FGFR3_AA==647],"FAMILY_TOT"]=1;
  data_mutfreq_all[data_mutfreq_all[FGFR3_AA==647],"FGFR3_TOT"]=1;
  return(data_mutfreq_all);
}

# This reduced version is used for multi reporting
sc.load_fgfr_fam_mut_freq_small <- function (dbConn)
  {

    data_mutfreq_small <-data.table(dbGetQuery(dbConn,"SELECT * FROM FGFR.VW_MUT_FREQ_REP_6_SMALL WHERE FGFR3_AA>396"));
    setkey(data_mutfreq_small,"GLOBAL_AL");
    data_mutfreq_small <- sc.add_manual_mut_changes(data_mutfreq_small);
    return(data_mutfreq_small);
  }

# ------------------------------
# Family wide mutations- non-cancer (e.g. humsavar)
# ------------------------------
sc.load_fgfr_fam_muts_non_cancer <- function (dbConn)
  {

    # 08/07/2015 - put this in a view!  Except this always uses FGFR3 AA so alignment not correct for FGFR1,2 & 4!
    # data_mutfreq_non_cancer <-data.table(dbGetQuery(dbConn,"SELECT GAL.GLOBAL_AL, GAL.FGFR3_AA, MHV.GENE_SYMBOL, MHV.SWISS_PROT, REGEXP_REPLACE(MHV.AA_CHANGE,'(p.)*[A-za-z]','') AA_POS, MHV.AA_CHANGE, MHV.VARIANT_TYPE, MHV.FTID, MHV.DBSNP_OR_DISEASE_ID, MHV.DISEASE_NAME, MHV.IS_CANCER FROM FGFR.VW_FGFR_FAMILY_ALIGNMENT_2 GAL LEFT JOIN STAGING.STA_MUT_HUMVAR MHV ON CAST(REGEXP_REPLACE(MHV.AA_CHANGE,'(p.)*[A-za-z]','') AS NUMBER) = GAL.FGFR3_AA WHERE MHV.IMPORT_ID='S00051'AND MHV.GENE_SYMBOL LIKE 'FGFR%' AND MHV.IS_CANCER=-1 ORDER BY MHV.GENE_SYMBOL,GLOBAL_AL"));
    setkey(data_mutfreq_non_cancer,"GLOBAL_AL");
    return(data_mutfreq_non_cancer);
  }

# ------------------------------
# Family wide mutation frequencies - non-cancer (e.g. humsavar)
# ------------------------------

# 09/07/2015 - ** KNOWN BUG - THIS DOES NOT DO CORRECT ALIGNMENT AND MAPPING ONTO FGFR3 FOR PLOTS ***
# USE: VW_MUT_HUMSAVAR_FAM_FULL OR VW_MUT_HUMSAVAR_FAM_SIMPLE IN FUTURE
#


# 08/07/2015 - ** Not sure if this is useful - reporting on germline is just done by marking diseases at a given position - there's no way to get a frequency in the
# way of cancer genomics, hence commented and not wrapped in VIEW **
# On second thoughts, the group by does stop multiple lines being reported
# This was used to generate the 1st germline mutation plots for FGFR1-4 to match the cancer mut_freqs plots
# Uses FGFR3 AA so alignment not correct for FGFR1,2 & 4!
sc.load_fgfr_fam_mut_freq_non_cancer <- function (dbConn)
  {
#     data_mutfreq_non_cancer <-data.table(dbGetQuery(dbConn,"  SELECT CAST(GAL.GLOBAL_AL AS NUMBER) AS GAL,
#       GAL.FGFR3_AA, MHV.GENE_SYMBOL,
#       REGEXP_REPLACE(MHV.AA_CHANGE,'(p.)*[A-za-z]','') AA_POS, count(*) freq
# FROM FGFR.VW_FGFR_FAMILY_ALIGNMENT_2 GAL LEFT JOIN STAGING.STA_MUT_HUMVAR MHV ON CAST(REGEXP_REPLACE(MHV.AA_CHANGE,'(p.)*[A-za-z]','') AS NUMBER) = GAL.FGFR3_AA
# WHERE MHV.IMPORT_ID='S00051'AND MHV.GENE_SYMBOL LIKE 'FGFR%' AND MHV.IS_CANCER=-1
# group by cast(GAL.GLOBAL_AL as number), GAL.FGFR3_AA, MHV.GENE_SYMBOL, REGEXP_REPLACE(MHV.AA_CHANGE,'(p.)*[A-za-z]','')
# ORDER BY MHV.GENE_SYMBOL,cast(GLOBAL_AL as number)"));
#     setkey(data_mutfreq_non_cancer,"GAL");
#     return(data_mutfreq_non_cancer);
  return("defunct/legacy");
  }


# ------------------------------
# Harshnira's mutations of interest
# ------------------------------
sc.load_muts_of_interest_h <- function (data.dir)
  {
  data_Hmuts <- data.table(read.csv(paste(data.dir,"table1_nospace_apr15_new.csv",sep="/")));
  # Specific muts for experiment only...
  data_HmutsSpecific <- data.table(read.csv(paste(data.dir,"table1_specific_muts_nospace.csv",sep="/")));
  # Tidy data
  data_Hmuts[,"category" := "All_Table1"];
  data_HmutsSpecific[,"category" := "Specific_exp_muts"];
  data_mutsOfInterest <- rbind(data_Hmuts,data_HmutsSpecific);

  rm(data_Hmuts, data_HmutsSpecific);
  return(data_mutsOfInterest);
  }


    # ------------------------------
    # ClustalW alignment of FGFR1-4
    # ------------------------------
sc.load_clustalw_alignment <-function(data.dir)
{
    data_clustalw <- data.table(read.csv(paste(data.dir,"fgfr1-4_global_clustalw.csv",sep="/")));
    setkey(data_clustalw,GLOBAL_AL);
    return(data_clustalw);
}
