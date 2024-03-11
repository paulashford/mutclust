# Utility functions for the structural clustering and null randomisation models
# Ash 10th March 2015

sc.calc_struct_cluster_residues <- function(h, r, pdb_xyz) {
# Nuts and bolts function for struct clustering
# Given the pdb, this fn returns a set of resnos that have atoms within
# a radius r of hotspot (=residue to search aroung) h.
# h is a residue in the pdb.

	if (h %in% unique(pdb_xyz$resno) == FALSE){
		return ("Error: Search residue (h) is not in pdb.");
	}

	# What are atom (x,y,z) of search residue?
	h_atoms <- pdb_xyz[resno ==h];

	# Which atoms in pdb are within rA of any hotspot(=search) atom?
	# * Note there is also a $calpha logical in pdb, which may be an alternative to all-atom...
	atom_dist = fields.rdist.near(pdb_xyz[,.(x,y,z)], h_atoms[,.(x,y,z)], delta=r);

	# From the returned index (of seq no of pdb within r) what atoms records are they?
	# Note: If only one atom, seem to get a vector for $ind, not an array, hence
	# this test to avoid error...
	if (is.null(dim(atom_dist$ind))){
			cluster_atoms = pdb_xyz[unique(atom_dist$ind[[1]]),];
		}else{
			cluster_atoms = pdb_xyz[unique(atom_dist$ind[,1]),];
		};

	# Return the set of residues from the atoms:
	return (unique(cluster_atoms$resno));
};

sc.rX <- function (n,p)
# Generate n observations of random binomial variable X, with P(X=1)=p
#Adapted from http://dept.stat.lsa.umich.edu/~jasoneg/Stat406/lab5.pdf
{
  U <- runif(n)

  X <- rep(0,n)

  w1 <- which (U <= p)
  X[w1] <- 1

  w2 <- which (U > p)
  X[w2] <- 0

return (X)

}

sc.scatter_muts <- function ( seq_length_n, p, tot_muts ) {
	# 26/06/2018
	# Improve MutClust background model by scattering all mutations
	# mapped to PDB to structure, rather than just marking amino as mutated
	# or not.  May need to be careful with huge hotspots (which increase overall background
	# count everywhere when scattered) but not dealing with that now.

	# Calls sc.rX multiple times, until tot_muts "used up".
	# Sums the resulting random seqs. 
	# Randomly removes element from the last iteration, so as to not exceed
	# tot_muts

	# Initialise sequence
	seq_scattered 	<- sc.rX ( seq_length_n, 0 );

	muts_scattered <- 0;
	while ( muts_scattered < tot_muts ) {

		# Run a scatter trial seq (of 0s and 1s)
		scatter_trial <- sc.rX ( seq_length_n, p );

		# Add trial to seq 
		seq_scattered <- seq_scattered + scatter_trial ;
		
		# How many muts so far?
		muts_scattered <- sum( seq_scattered );

	}

	# Now subtract off any excess mutations
	num_muts_remove <- sum ( seq_scattered ) - tot_muts ;
	if ( num_muts_remove > 0 ) {
		
		muts_remove <- sample( 1:seq_length_n, num_muts_remove );
		seq_scattered[ muts_remove ] <- seq_scattered[ muts_remove ] - 1;	
	}
	
	return ( seq_scattered );

}

# --------------------------------------------------------------------------
# Hill function mut normailization
# "Compress" range of mutations on residue hotspots
# with theta = 2 and exponent = 3 will give range [0,1) over first 6 mutations

# Method used in 
# https://www.pnas.org/content/pnas/112/40/E5486.full.pdf
# https://www.pnas.org/content/pnas/suppl/2015/09/17/1516373112.DCSupplemental/pnas.1516373112.sapp.pdf
# Suppl fig 13
# theta =2, m=3
# --------------------------------------------------------------------------
mut_weight <- function( num_muts, theta, exponent ){
	return( 
			( num_muts^exponent ) / ( ( theta^exponent ) + ( num_muts^exponent ) )             
		)
	}



# A simple function to calculate the Jaccard index for two sets a and b.
jaccard <- function(a,b)
{
	length(intersect(a,b)) / length(union(a,b))
}

# Creating cluster centroids lists... (ad-hoc)
sc.get_cluster_centroids <- function(reslist,pdb,outfile){
	pdb_clusters <-pdb[resno %in% reslist];
	# What are the mean atomic positions of each residue?
	pdb_cluster_means<-pdb_clusters[,.(mean(x),mean(y),mean(z)), by=resno]
	write.csv(pdb_cluster_means,file=outfile, row.names=F);
}


# Amino acid - get the letter given 3 char name...
aminoAcidLetter <- function(aminoThreeLetterCode) {
	require(Biostrings);
	return(names(AMINO_ACID_CODE[tolower(AMINO_ACID_CODE)==tolower(aminoThreeLetterCode)]));
}

# Basic quantile colouring
# Takes bottom two and makes greens, top two reds
# Ignores rest...
quantileColourFOLDX <- function(foldQuants, val){
	if (is.na(val)) {return("BLACK")}
	col="WHITE";
	noQuants <- length(foldQuants);

	if (noQuants>3){
			if (val>foldQuants[[noQuants]]) {return("RED")}
			if (val>foldQuants[[noQuants-1]]) {return("LIGHTRED")}
			if (val<foldQuants[[1]]) {return("GREEN")}
			if (val<foldQuants[[2]]) {return("LIGHTGREEN")}
			return("WHITE");
		}else{
			return("WHITE");
		}

	}
