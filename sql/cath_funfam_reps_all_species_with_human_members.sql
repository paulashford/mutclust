SELECT
	FF.superfamily_id,
	FF.funfam_number,
	seed_alignment,
	seed_dops_score,
	name,
	rep_id,
	inclusion_bitscore,
	inclusion_e_value,
	num_rep_annotations,
	rep_source_id,
	num_members_in_seed_aln,
	num_members_in_funfam,
	fug.sequence_md5,
	fug.member_id,
	fug.taxon_id
	
FROM
	funfam FF
INNER JOIN
	(SELECT  DISTINCT
		SUPERFAMILY_ID,
		FUNFAM_NUMBER,
		SEQUENCE_MD5,
		MEMBER_ID,
		TAXON_ID
	FROM FUNFAM_TO_UNIPROT_GENE 
	WHERE TAXON_ID=9606) FUG
	ON FF.superfamily_id = FUG.superfamily_id
	AND FF.funfam_number = FUG.funfam_number
	