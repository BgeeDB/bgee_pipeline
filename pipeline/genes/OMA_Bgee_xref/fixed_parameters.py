oma_entrez_url = 'https://omabrowser.org/All/oma-entrez.txt.gz'
oma_ensembl_url = 'https://omabrowser.org/All/oma-ensembl.txt.gz'
gene_ensembl_url = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz'

#SQL to query cross-reference files
sql_query = '''
      SELECT geneId, `OMA ID`, speciesId FROM bgee_genes as bgee JOIN oma_ensembl as ol ON ol.`Ensembl ID` = bgee.geneId
      UNION
      SELECT geneId, `OMA ID`, speciesId FROM bgee_genes as bgee JOIN oma_entrez as oz ON oz.`EntrezGene ID` = bgee.geneId
      UNION
      select geneId, `OMA ID`, speciesId FROM 
       (select bgee.geneId, n.GeneID as ncbi_id, speciesId  from (
       SELECT geneId, speciesId  FROM bgee_genes  where speciesId in (105023, 8364)) as bgee JOIN
        ncbi_ensembl as n ON n.Ensembl_gene_identifier = bgee.geneId) as bgee_ncbi 
       JOIN oma_entrez as oz ON oz.`EntrezGene ID` = bgee_ncbi.ncbi_id
       '''
drop_database = True
