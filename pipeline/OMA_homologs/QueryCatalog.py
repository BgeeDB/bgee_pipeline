class QueryCatalogOMA:
    def __init__(self, prefixes: str = None) -> object:
        """The QueryCatalogOMA class constructor.

        :param prefixes: The SPARQL query prefixes. 
            Default:       
            'PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
            PREFIX dct: <http://purl.org/dc/terms/>
            PREFIX obo: <http://purl.obolibrary.org/obo/>
            PREFIX ensembl: <http://rdf.ebi.ac.uk/resource/ensembl/>
            PREFIX oma: <http://omabrowser.org/ontology/oma#>
            PREFIX orth: <http://purl.org/net/orth#>
            PREFIX sio: <http://semanticscience.org/resource/>
            PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
            PREFIX lscr: <http://purl.org/lscr#>'
        """
        if prefixes is None:
            self.prefixes = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX dct: <http://purl.org/dc/terms/>
PREFIX obo: <http://purl.obolibrary.org/obo/>
PREFIX ensembl: <http://rdf.ebi.ac.uk/resource/ensembl/>
PREFIX oma: <http://omabrowser.org/ontology/oma#>
PREFIX orth: <http://purl.org/net/orth#>
PREFIX sio: <http://semanticscience.org/resource/>
PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
PREFIX lscr: <http://purl.org/lscr#>

"""
        else:
            self.prefixes = prefixes

    def get_all_prefixes(self) -> str:
        """Get query prefixes using SPARQL language.

        :return: str : SPARQL query prefixes
        """
        return self.prefixes

    def query_orthologs_per_species(self, species1: str, species2: str,
                                    gene_xref_lscr_term_s1: str, gene_xref_lscr_term_s2: str) -> str:
        """Get the OMA SPARQL query to retrieve orthology information for a species pair by exclusively considering
            OMA genes with given cross-reference properties.

        :param species1: A first either NCBI or Ensembl species of reference for orthology relationship.
            The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/9606" (i.e. homo sapiens).
        :param species2: A second either NCBI or Ensembl species of reference for orthology relationship.
            The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/9606" (i.e. homo sapiens).
        :param gene_xref_lscr_term_s1: it is the cross-reference property for species1 genes used
            in OMA RDF store. It must be an IRI value. It can also take advantage of self.prefixes,
             e.g. "lscr:xrefNCBIGene".
        :param gene_xref_lscr_term_s2: it is the cross-reference property for species2 genes used
            in OMA RDF store. It must be an IRI value. It can also take advantage of self.prefixes,
             e.g. "lscr:xrefNCBIGene".
        :return: str : The OMA SPARQL query.
        """
        query_OMA = self.get_all_prefixes() + """
select distinct ?gene1 ?gene2  ?tax_level ?tax_id
where {
values ?species1 {<""" + species1 + """>}
values ?species2 {<""" + species2 + """>}
?cluster a orth:OrthologsCluster.
?cluster orth:hasTaxonomicRange ?taxRange.
?taxRange orth:taxRange ?tax_level;
          orth:taxRangeId ?tax_id.
?cluster orth:hasHomologousMember ?node1.
?cluster orth:hasHomologousMember ?node2. 
?node2 orth:hasHomologousMember* ?protein2. 
?node1 orth:hasHomologousMember* ?protein1.
?protein1 sio:SIO_010079 ?gene1_uri.
?gene1_uri """ + gene_xref_lscr_term_s1 + """/dct:identifier  ?gene1.
?protein1 orth:organism/obo:RO_0002162 ?species1 .
?protein2 sio:SIO_010079 ?gene2_uri.
?gene2_uri """ + gene_xref_lscr_term_s2 + """/dct:identifier ?gene2. 
?protein2 orth:organism/obo:RO_0002162 ?species2 .
filter(?node1 != ?node2  ) 
}
"""
        return query_OMA

    # species1 or species2 must be a drosophila-like species
    def query_flybase_orthologs_per_species(self, species1: str, species2: str,
                                            gene_xref_lscr_term: str) -> str:
        """Get the OMA SPARQL query to retrieve orthology information for a species pair where one of the species is
        from FlyBase.  It solely considers a cross-reference property for genes of a non-Flybase species.

        :param species1: A first (NCBI, Ensembl, or FlyBase) species of reference for orthology relationship.
            The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/9606" (i.e. homo sapiens).
            Note that either species1 or species2 must be a FlyBase species.
        :param species2: A second (NCBI, Ensembl, or FlyBase) species of reference for orthology relationship.
            The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/9606" (i.e. homo sapiens). Note
            that either species1 or species2 must be a FlyBase species.
        :param gene_xref_lscr_term: it is the cross-reference property for the non-FlyBase genes used
            in OMA RDF store. It must be an IRI value. It can also take advantage of self.prefixes,
             e.g. "lscr:xrefNCBIGene".
        :return: str : The OMA SPARQL query.
        """
        query_OMA = self.get_all_prefixes() + """
select distinct ?gene1 ?gene2 ?tax_level ?tax_id
where {
values ?species1 {<""" + species1 + """>}
values ?species2 {<""" + species2 + """>}
?cluster a orth:OrthologsCluster.
?cluster orth:hasTaxonomicRange ?taxRange.
?taxRange orth:taxRange ?tax_level;
          orth:taxRangeId ?tax_id.
?cluster orth:hasHomologousMember ?node1.
?cluster orth:hasHomologousMember ?node2. 
?node2 orth:hasHomologousMember* ?protein2. 
?node1 orth:hasHomologousMember* ?protein1.
?protein1 sio:SIO_010079 ?gene1_uri.
{?gene1_uri """ + gene_xref_lscr_term + """/dct:identifier  ?gene1.} UNION
{?gene1_uri dct:identifier  ?gene1.}
?protein1 orth:organism/obo:RO_0002162 ?species1 .
?protein2 sio:SIO_010079 ?gene2_uri.
{?gene2_uri """ + gene_xref_lscr_term + """/dct:identifier ?gene2. } UNION
{?gene2_uri dct:identifier  ?gene2.} 
?protein2 orth:organism/obo:RO_0002162 ?species2 .
filter(?node1 != ?node2  ) 
}
"""
        return query_OMA

    # species1 must be a wormbase-like species
    def query_orthologs_per_species_filter_gene_prefix(self, species1: str, species2: str, gene_xref_lscr_term_s2: str,
                                                       gene1_prefix: str, swap_gene_projection: bool = False) -> str:
        """Get the OMA SPARQL query to retrieve orthology information for a species pair. Where species1 gene
        identifiers in OMA starts with a given gene_prefix and species2 genes states a given cross-reference property.

        :param species1: A first species of reference for orthology relationship.
            The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/9606" (i.e. homo sapiens).
        :param species2: A second either NCBI or Ensembl species of reference for orthology relationship.
            The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/9606" (i.e. homo sapiens).
        :param gene1_prefix: species1 gene/protein identifier prefix.
        :param gene_xref_lscr_term_s2: it is the cross-reference property for species2 genes used
            in OMA RDF store. It must be an IRI value. It can also take advantage of self.prefixes,
             e.g. "lscr:xrefNCBIGene".
        :param swap_gene_projection: change order of orthologous gene projection, default is False.
        :return: str : The OMA SPARQL query.
        """
        if swap_gene_projection:
            gene_projection = "(?gene_2 as ?gene1) (?gene_1 as ?gene2)"
        else:
            gene_projection = "(?gene_1 as ?gene1) (?gene_2 as ?gene2)"
        query_OMA = self.get_all_prefixes() + """
select distinct """ + gene_projection + """ ?tax_level ?tax_id
where {
values ?species1 {<""" + species1 + """>}
values ?species2 {<""" + species2 + """>}
?cluster a orth:OrthologsCluster.
?cluster orth:hasTaxonomicRange ?taxRange.
?taxRange orth:taxRange ?tax_level;
          orth:taxRangeId ?tax_id.
?cluster orth:hasHomologousMember ?node1.
?cluster orth:hasHomologousMember ?node2. 
?node2 orth:hasHomologousMember* ?protein2. 
?node1 orth:hasHomologousMember* ?protein1.
?protein1 dct:identifier  ?gene_1.
?protein1 orth:organism/obo:RO_0002162 ?species1 .
?protein2 sio:SIO_010079 ?gene2_uri.
?gene2_uri """ + gene_xref_lscr_term_s2 + """/dct:identifier ?gene_2. 
?protein2 orth:organism/obo:RO_0002162 ?species2 .
filter(?node1 != ?node2)
filter regex(?gene_1, '""" + gene1_prefix + """[0-9]*$') 
}
"""
        return query_OMA

        # species1 and 2 must be a special case species
    def query_orthologs_per_species_filter_gene_prefixes(self, species1: str, species2: str,
                                                         gene1_prefix: str, gene2_prefix: str,
                                                         swap_gene_projection: bool = False) -> str:
        """Get the OMA SPARQL query to retrieve orthology information for a species pair. Where species1 gene/protein
        identifiers starts with a given gene1_prefix and species2 genes/proteins starts with a given gene2_prefix.
        The returned query is to process orthology information of species defined in the species_to_gene_id_prefix
        config.properties file parameter.

        :param species1: A first species of reference for orthology relationship.
            The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/9606" (i.e. homo sapiens).
        :param species2: A second species of reference for orthology relationship.
            The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/9606" (i.e. homo sapiens).
        :param gene1_prefix: species1 gene/protein identifier prefix.
            WARNING: If the prefix is 'ENS' referring to Ensembl species genes, use other query method.
        :param gene2_prefix: species2 gene/protein identifier prefix.
            WARNING: If the prefix is 'ENS' referring to Ensembl species genes, use other query method.
        :param swap_gene_projection: change order of orthologous gene projection, default is False.
        :return: str : The OMA SPARQL query.
        """
        if swap_gene_projection:
            gene_projection = "(?gene_2 as ?gene1) (?gene_1 as ?gene2)"
        else:
            gene_projection = "(?gene_1 as ?gene1) (?gene_2 as ?gene2)"
        query_OMA = self.get_all_prefixes() + """
    select distinct """ + gene_projection + """ ?tax_level ?tax_id
    where {
    values ?species1 {<""" + species1 + """>}
    values ?species2 {<""" + species2 + """>}
    ?cluster a orth:OrthologsCluster.
    ?cluster orth:hasTaxonomicRange ?taxRange.
    ?taxRange orth:taxRange ?tax_level;
              orth:taxRangeId ?tax_id.
    ?cluster orth:hasHomologousMember ?node1.
    ?cluster orth:hasHomologousMember ?node2. 
    ?node2 orth:hasHomologousMember* ?protein2. 
    ?node1 orth:hasHomologousMember* ?protein1.
    ?protein1 dct:identifier  ?gene_1.
    ?protein1 orth:organism/obo:RO_0002162 ?species1 .
    ?protein2 dct:identifier ?gene_2. 
    ?protein2 orth:organism/obo:RO_0002162 ?species2 .
    filter(?node1 != ?node2)
    filter regex(?gene_1, '""" + gene1_prefix + """[0-9]*$') 
    filter regex(?gene_2, '""" + gene2_prefix + """[0-9]*$') 
    }
    """
        return query_OMA

    def query_flybase_orthologs_per_species_filter_gene_prefix(self, species1: str, species2: str,
                                                               gene1_prefix: str, swap_gene_projection: bool = False
                                                               ) -> str:
        """Get the OMA SPARQL query to retrieve orthology information for a species pair. Where species1 gene/protein
        identifiers starts with a given gene1_prefix different from 'ENS' and species2 is a FlyBase species.
        The returned query is to process orthology information where one of the species is a special case and it is
        defined in the species_to_gene_id_prefix config.properties file parameter.

        :param species1: A first species of reference for orthology relationship.
            The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/9606" (i.e. homo sapiens).
        :param species2: A second Flybase species of reference for orthology relationship.
            The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/9606" (i.e. homo sapiens).
        :param gene1_prefix: species1 gene/protein identifier prefix.
            WARNING: If the prefix is 'ENS' referring to Ensembl species genes, use other query method.
        :param swap_gene_projection: change order of orthologous gene projection, default is False.
        :return: str : The OMA SPARQL query.
        """
        if swap_gene_projection:
            gene_projection = "(?gene_2 as ?gene1) (?gene_1 as ?gene2)"
        else:
            gene_projection = "(?gene_1 as ?gene1) (?gene_2 as ?gene2)"
        query_OMA = self.get_all_prefixes() + """
select distinct """ + gene_projection + """ ?tax_level ?tax_id
where {
values ?species1 {<""" + species1 + """>}
values ?species2 {<""" + species2 + """>}
?cluster a orth:OrthologsCluster.
?cluster orth:hasTaxonomicRange ?taxRange.
?taxRange orth:taxRange ?tax_level;
          orth:taxRangeId ?tax_id.
?cluster orth:hasHomologousMember ?node1.
?cluster orth:hasHomologousMember ?node2. 
?node2 orth:hasHomologousMember* ?protein2. 
?node1 orth:hasHomologousMember* ?protein1.
?protein1 dct:identifier  ?gene_1.
?protein1 orth:organism/obo:RO_0002162 ?species1 .
?protein2 sio:SIO_010079 ?gene2_uri.
?gene2_uri dct:identifier  ?gene_2.
?protein2 orth:organism/obo:RO_0002162 ?species2 .
filter(?node1 != ?node2)
filter regex(?gene_1, '""" + gene1_prefix + """[0-9]*$') 
}
"""
        return query_OMA

    def query_flybase_species_paralogs(self, species1: str, species2: str) -> str:
        """Get the OMA SPARQL query to retrieve paralogy information for a species pair. Where species1 and species2
         are both of them FlyBase species.
        :param species1: A first FlyBase species of reference for paralogy relationship.
            The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/7217".
        :param species2: A second Flybase species of reference for paralogy relationship.
            The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/7217".
        :return: str : The OMA SPARQL query.
        """
        query_OMA = self.get_all_prefixes() + """
select distinct ?gene1 ?gene2  ?tax_level ?tax_id
where {
values ?species1 {<""" + species1 + """>}
values ?species2 {<""" + species2 + """>}
?cluster_speciation a orth:OrthologsCluster.
?cluster orth:hasHomologousMember ?cluster_speciation.
?cluster_speciation orth:hasTaxonomicRange ?taxRange.
?taxRange orth:taxRange ?tax_level;
          orth:taxRangeId ?tax_id.
?cluster a orth:ParalogsCluster.
?cluster orth:hasHomologousMember ?node1.
?cluster orth:hasHomologousMember ?node2. 
?node2 orth:hasHomologousMember* ?protein2. 
?node1 orth:hasHomologousMember* ?protein1.
?protein1 sio:SIO_010079 ?gene1_uri.
?gene1_uri dct:identifier  ?gene1.
?protein1 orth:organism/obo:RO_0002162 ?species1 .
?protein2 sio:SIO_010079 ?gene2_uri.
?gene2_uri dct:identifier  ?gene2. 
?protein2 orth:organism/obo:RO_0002162 ?species2 .
filter(?node1 != ?node2  ) 
}
"""
        return query_OMA

    def query_flybase_paralogs_per_species(self, species1: str, species2: str, gene_xref_lscr_term_s2: str) -> str:
        """Get the OMA SPARQL query to retrieve paralogy information for a species pair.
        Where species1 must be a FlyBase species and species2 a non-FlyBase species.

        :param species1: A first species of reference for paralogy relationship.
                    The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/7217".
        :param species2: A second either NCBI or Ensembl species of reference for paralogy relationship.
                    The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/9606" (i.e. homo sapiens).
        :param gene_xref_lscr_term_s2: it is the cross-reference property for species2 genes used
                    in OMA RDF store. It must be an IRI value. It can also take advantage of self.prefixes,
                     e.g. "lscr:xrefNCBIGene".
        :return: str : The OMA SPARQL query.
        """
        query_OMA = self.get_all_prefixes() + """
select distinct ?gene1 ?gene2  ?tax_level ?tax_id
where {
values ?species1 {<""" + species1 + """>}
values ?species2 {<""" + species2 + """>}
?cluster_speciation a orth:OrthologsCluster.
?cluster orth:hasHomologousMember ?cluster_speciation.
?cluster_speciation orth:hasTaxonomicRange ?taxRange.
?taxRange orth:taxRange ?tax_level;
          orth:taxRangeId ?tax_id.
?cluster a orth:ParalogsCluster.
?cluster orth:hasHomologousMember ?node1.
?cluster orth:hasHomologousMember ?node2. 
?node2 orth:hasHomologousMember* ?protein2. 
?node1 orth:hasHomologousMember* ?protein1.
?protein1 sio:SIO_010079 ?gene1_uri.
?gene1_uri dct:identifier  ?gene1.
?protein1 orth:organism/obo:RO_0002162 ?species1 .
?protein2 sio:SIO_010079 ?gene2_uri.
?gene2_uri """ + gene_xref_lscr_term_s2 + """/dct:identifier ?gene2. 
?protein2 orth:organism/obo:RO_0002162 ?species2 .
filter(?node1 != ?node2  ) 
}
"""
        return query_OMA

    # This query/function is not applicable to flybase and wormbase species
    def query_paralogs_per_species(self, species1: str, species2: str,
                                   gene_xref_lscr_term_s1: str, gene_xref_lscr_term_s2: str) -> str:
        """Get the OMA SPARQL query to retrieve paralogy information for a species pair by exclusively considering
                    OMA genes with given cross-reference properties.

        :param species1: A first either NCBI or Ensembl species of reference for paralogy relationship.
                    The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/9606" (i.e. homo sapiens).
        :param species2: A second either NCBI or Ensembl species of reference for paralogy relationship.
                    The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/9606" (i.e. homo sapiens).
        :param gene_xref_lscr_term_s1: it is the cross-reference property for species1 genes used
                    in OMA RDF store. It must be an IRI value. It can also take advantage of self.prefixes,
                     e.g. "lscr:xrefNCBIGene".
        :param gene_xref_lscr_term_s2: it is the cross-reference property for species2 genes used
                    in OMA RDF store. It must be an IRI value. It can also take advantage of self.prefixes,
                     e.g. "lscr:xrefNCBIGene".
        :return: str : The OMA SPARQL query.
                """
        query_OMA = self.get_all_prefixes() + """
select distinct ?gene1 ?gene2  ?tax_level ?tax_id
where {
values ?species1 {<""" + species1 + """>}
values ?species2 {<""" + species2 + """>}
?cluster_speciation a orth:OrthologsCluster.
?cluster orth:hasHomologousMember ?cluster_speciation.
?cluster_speciation orth:hasTaxonomicRange ?taxRange.
?taxRange orth:taxRange ?tax_level;
          orth:taxRangeId ?tax_id.
?cluster a orth:ParalogsCluster.
?cluster orth:hasHomologousMember ?node1.
?cluster orth:hasHomologousMember ?node2. 
?node2 orth:hasHomologousMember* ?protein2. 
?node1 orth:hasHomologousMember* ?protein1.
?protein1 sio:SIO_010079 ?gene1_uri.
?gene1_uri """ + gene_xref_lscr_term_s1 + """/dct:identifier  ?gene1.
?protein1 orth:organism/obo:RO_0002162 ?species1 .
?protein2 sio:SIO_010079 ?gene2_uri.
?gene2_uri """ + gene_xref_lscr_term_s2 + """/dct:identifier ?gene2. 
?protein2 orth:organism/obo:RO_0002162 ?species2 .
filter(?node1 != ?node2  ) 
}
"""
        return query_OMA

    def query_flybase_paralogs_per_species_filter_gene_prefix(self, species1: str, species2: str,
                                                              gene_prefix: str = None,
                                                              swap_gene_projection: bool = False) -> str:
        """Get the OMA SPARQL query to retrieve paralogy information for a species pair. Where species1 gene/protein
        identifiers starts with a given gene1_prefix different from 'ENS' and species2 is a FlyBase species.
        The returned query is to process paralogy information where one of the species is a special case and it is
        defined in the species_to_gene_id_prefix config.properties file parameter.

        :param species1: A first species of reference for paralogy relationship.
            The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/9606" (i.e. homo sapiens).
        :param species2: A second Flybase species of reference for paralogy relationship.
            The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/9606" (i.e. homo sapiens).
        :param gene1_prefix: species1 gene/protein identifier prefix.
            WARNING: If the prefix is 'ENS' referring to Ensembl species genes, use other query method.
        :param swap_gene_projection: change order of paralogous gene projection, default is False.
        :return: str : The OMA SPARQL query.
        """
        if swap_gene_projection:
            gene_projection = "(?gene_2 as ?gene1) (?protein1 as ?gene2)"
        else:
            gene_projection = "(?protein1 as ?gene1) (?gene_2 as ?gene2)"
        query_OMA = self.get_all_prefixes() + """
select distinct """ + gene_projection + """ ?tax_level ?tax_id
where {
values ?species1 {<""" + species1 + """>}
values ?species2 {<""" + species2 + """>}
?cluster_speciation a orth:OrthologsCluster.
?cluster orth:hasHomologousMember ?cluster_speciation.
?cluster_speciation orth:hasTaxonomicRange ?taxRange.
?taxRange orth:taxRange ?tax_level;
          orth:taxRangeId ?tax_id.
?cluster a orth:ParalogsCluster.
?cluster orth:hasHomologousMember ?node1.
?cluster orth:hasHomologousMember ?node2. 
?node2 orth:hasHomologousMember* ?protein2. 
?node1 orth:hasHomologousMember* ?protein1.
#?protein1 dct:identifier  ?gene1.
?protein1 orth:organism/obo:RO_0002162 ?species1 .
?protein2 sio:SIO_010079 ?gene2_uri.
?gene2_uri dct:identifier  ?gene_2.
?protein2 orth:organism/obo:RO_0002162 ?species2 .
filter(?node1 != ?node2)
#filter regex(?gene1, '""" + gene_prefix + """[0-9]*$')
}
"""
        return query_OMA

    def query_paralogs_per_species_filter_gene_prefix(self, species1: str, species2: str, gene_xref_lscr_term_s2: str,
                                                      gene_prefix: str = None, swap_gene_projection: bool = False
                                                      ) -> str:
        """Get the OMA SPARQL query to retrieve paralogy information for a species pair. Where species1 gene
        identifiers in OMA starts with a given gene_prefix and species2 genes states a given cross-reference property.

        :param species1: A first species of reference for paralogy relationship.
            The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/9606" (i.e. homo sapiens).
        :param species2: A second either NCBI or Ensembl species of reference for paralogy relationship.
            The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/9606" (i.e. homo sapiens).
        :param gene1_prefix: species1 gene/protein identifier prefix.
        :param gene_xref_lscr_term_s2: it is the cross-reference property for species2 genes used
            in OMA RDF store. It must be an IRI value. It can also take advantage of self.prefixes,
             e.g. "lscr:xrefNCBIGene".
        :param swap_gene_projection: change order of paralogous gene projection, default is False.
        :return: str : The OMA SPARQL query.
        """
        if swap_gene_projection:
            gene_projection = "(?gene_2 as ?gene1) (?protein1 as ?gene2)"
        else:
            gene_projection = "(?protein1 as ?gene1) (?gene_2 as ?gene2)"
        query_OMA = self.get_all_prefixes() + """
select distinct """ + gene_projection + """ ?tax_level ?tax_id
where {
values ?species1 {<""" + species1 + """>}
values ?species2 {<""" + species2 + """>}
?cluster_speciation a orth:OrthologsCluster.
?cluster orth:hasHomologousMember ?cluster_speciation.
?cluster_speciation orth:hasTaxonomicRange ?taxRange.
?taxRange orth:taxRange ?tax_level;
          orth:taxRangeId ?tax_id.
?cluster a orth:ParalogsCluster.
?cluster orth:hasHomologousMember ?node1.
?cluster orth:hasHomologousMember ?node2. 
?node2 orth:hasHomologousMember* ?protein2. 
?node1 orth:hasHomologousMember* ?protein1.
#?protein1 dct:identifier  ?gene1.
?protein1 orth:organism/obo:RO_0002162 ?species1 .
?protein2 sio:SIO_010079 ?gene2_uri.
?gene2_uri """ + gene_xref_lscr_term_s2 + """/dct:identifier ?gene_2. 
?protein2 orth:organism/obo:RO_0002162 ?species2 .
filter(?node1 != ?node2)
#filter regex(?gene1, '""" + gene_prefix + """[0-9]*$')
}
"""
        return query_OMA

    def query_paralogs_per_species_filter_gene_prefixes(self, species1: str, species2: str,
                                                        gene1_prefix: str = None, gene2_prefix: str = None,
                                                        swap_gene_projection: bool = False) -> str:
        """Get the OMA SPARQL query to retrieve paralogy information for a species pair. Where species1 gene/protein
        identifiers starts with a given gene1_prefix and species2 genes/proteins starts with a given gene2_prefix.
        The returned query is to process paralogy information of species defined in the species_to_gene_id_prefix
        config.properties file parameter.

        :param species1: A first species of reference for paralogy relationship.
            The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/9606" (i.e. homo sapiens).
        :param species2: A second species of reference for paralogy relationship.
            The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/9606" (i.e. homo sapiens).
        :param gene1_prefix: species1 gene/protein identifier prefix.
            WARNING: If the prefix is 'ENS' referring to Ensembl species genes, use other query method.
        :param gene2_prefix: species2 gene/protein identifier prefix.
            WARNING: If the prefix is 'ENS' referring to Ensembl species genes, use other query method.
        :param swap_gene_projection: change order of paralogous gene projection, default is False.
        :return: str : The OMA SPARQL query.
        """
        if swap_gene_projection:
            gene_projection = "(?protein2 as ?gene1) (?protein1 as ?gene2)"
        else:
            gene_projection = "(?protein1 as ?gene1) (?protein2 as ?gene2)"
        query_OMA = self.get_all_prefixes() + """
select distinct """ + gene_projection + """ ?tax_level ?tax_id
where {
values ?species1 {<""" + species1 + """>}
values ?species2 {<""" + species2 + """>}
?cluster_speciation a orth:OrthologsCluster.
?cluster orth:hasHomologousMember ?cluster_speciation.
?cluster_speciation orth:hasTaxonomicRange ?taxRange.
?taxRange orth:taxRange ?tax_level;
          orth:taxRangeId ?tax_id.
?cluster a orth:ParalogsCluster.
?cluster orth:hasHomologousMember ?node1.
?cluster orth:hasHomologousMember ?node2. 
?node2 orth:hasHomologousMember* ?protein2. 
?node1 orth:hasHomologousMember* ?protein1.
#?protein1 dct:identifier ?gene1.
?protein1 orth:organism/obo:RO_0002162 ?species1 .
#?protein2 dct:identifier ?gene2. 
?protein2 orth:organism/obo:RO_0002162 ?species2 .
filter(?node1 != ?node2)
#filter regex(?gene1, '""" + gene1_prefix + """[0-9]*$')
#filter regex(?gene2, '""" + gene2_prefix + """[0-9]*$')
}
"""
        return query_OMA

    def get_mapper_OMA_IRI_to_prefixed_id(self, species, id_prefix) -> str:
        query = self.get_all_prefixes() + """
            select * {
            ?protein a orth:Protein.
            ?protein dct:identifier ?gene. 
            ?protein orth:organism/obo:RO_0002162 <""" + species + """> .
            filter regex(?gene, '""" + id_prefix + """[0-9]*$')
            }"""
        return query

    def get_mapper_NCBI_ID_to_ENSEMBL_ID(self, species: str) -> str:
        """Query UniProt RDF store via OMA sparql endpoint to get a mapping from gene NCBI identifiers to
         Ensembl gene ids for a given species.
              :param species:
                  The species value must be an IRI, e.g. "http://purl.uniprot.org/taxonomy/105023"
                  (i.e. Turquoise killifish).
              :return: str : The OMA SPARQL query.
              """
        query = self.get_all_prefixes() + """
select * {service<https://sparql.uniprot.org/sparql/>{
SELECT (STRAFTER(str(?xref), "http://purl.uniprot.org/geneid/") as ?ncbi_gene) 
  (STRAFTER(str(?ensemblGeneId), "http://rdf.ebi.ac.uk/resource/ensembl/") as ?ens_gene)
WHERE
{
  ?protein a up:Protein .
  ?protein up:organism <""" + species + """>.
  ?protein rdfs:seeAlso ?xref.
  ?xref up:database <http://purl.uniprot.org/database/GeneID>.
  ?protein rdfs:seeAlso ?gene.
  ?gene up:transcribedFrom ?ensemblGeneId
  
}}} """
        return query

