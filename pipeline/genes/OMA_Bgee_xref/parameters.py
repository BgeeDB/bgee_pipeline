oma_entrez_file = '/Users/tarcisio/Downloads/oma-entrez.txt'
oma_ensembl_file = '/Users/tarcisio/Downloads/oma-ensembl.txt'
gene_ensembl_file = '/Users/tarcisio/Downloads/gene2ensembl.tsv'
bgee_genes_file = '/Users/tarcisio/Downloads/bgeeGenes.csv'

output_file_name = 'oma-bgee-xref'
drop_database = True

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