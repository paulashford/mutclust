-- mutclust_funfam_swissprot_members.sql
-- 11/10/2019      

    -----------------------------------------------------------------------------------------------------------
    -- How many human paralogues make up FunFams?
    -----------------------------------------------------------------------------------------------------------
    DROP VIEW funvar_tx.vw_funfam_swiss_member_count;
    CREATE VIEW funvar_tx.vw_funfam_swiss_member_count
    AS 
        SELECT
            fsg.superfamily_id,
            fsg.funfam_number,
            ff.name ff_name,
            COUNT( * ) num_ff_members
    
        FROM funvar_tx.vw_funfam_swissprot_gene  fsg 
        INNER JOIN funvar_import.cath_funfam ff 
            ON ff.superfamily_id = fsg.superfamily_id 
            AND ff.funfam_number = fsg.funfam_number
        
        WHERE fsg.taxon_id = 9606 AND fsg.status = 'reviewed'

        GROUP BY
            fsg.superfamily_id,
            fsg.funfam_number,
            ff.name
            
        ORDER BY
            COUNT( * ) DESC,
            fsg.superfamily_id, 
            fsg.funfam_number;
    