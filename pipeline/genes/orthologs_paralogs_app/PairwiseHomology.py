from SPARQLWrapper import SPARQLWrapper, CSV
from QueryCatalog import QueryCatalogOMA
import pandas


class PairwiseHomologyOMA:
    """A selector of the corresponding OMA SPARQL query for generating the pairwise homology relations between two
    given species.
    """
    def __init__(self, ncbi_species: list, flybase_species: list, taxon_namespace: str, ncbi_xref_rdf_property: str,
                 ensembl_xref_rdf_property: str, species2prefix: dict, query_catalog: QueryCatalogOMA = None) -> object:
        """PairwiseHomologyOMA class constructor

        :param ncbi_species: The NCBI species list.
        :param flybase_species: The FlyBase species list.
        :param taxon_namespace: The species namespace prefix.
        :param ncbi_xref_rdf_property: The NCBI cross-reference property used in the RDF store to be queried.
        :param ensembl_xref_rdf_property: The Ensembl cross-reference property used in the RDF store to be queried.
        :param species2prefix: A mapping between the species and a prefix used in their gene identifiers, for example,
         "BL" is prefix for gene identifiers of the Branchiostoma lanceolatum species in Ensembl metazoa database.
        :param query_catalog: the Query catalog object that contains the queries to be executed in the sparql endpoint.
        """
        self.ncbi_species = ncbi_species
        self.species_drosophila = flybase_species
        self.ncbi_xref_rdf_property = ncbi_xref_rdf_property
        self.ensembl_xref_rdf_property = ensembl_xref_rdf_property
        self.taxon_namespace = taxon_namespace
        self.species2prefix = {}
        if species2prefix is not None:
            for species_id in species2prefix:
                self.species2prefix.__setitem__(taxon_namespace + str(species_id),
                                            species2prefix.get(species_id))
        if query_catalog is None:
            self.query_catalog = QueryCatalogOMA()
        else:
            self.query_catalog = query_catalog
        if self.species2prefix.__len__() != 0:
            self.special_case_species = self.species2prefix.keys()
        else:
            self.special_case_species = None


    def get_query_OMA_orthologs(self, species1: str, species2: str) -> str:
        """Get the corresponding orthology related OMA SPARQL query from the self.query_catalog
         for a given pair of species.

        :param species1: A first species of reference for orthology relationship.
        :param species2: A second species of reference for orthology relationship.

        :return: str : the SPARQL query using the OMA data schema.
        """
        # Drosophila-like species gene ids starting with FbXXX are not considered as an Ensembl cross-reference in OMA.
        if species1 in self.species_drosophila or species2 in self.species_drosophila:
            if species1 in self.ncbi_species:
                query_OMA = self.query_catalog.query_flybase_orthologs_per_species(
                    species1, species2, self.ncbi_xref_rdf_property)
            else:
                if species2 in self.ncbi_species:
                    query_OMA = self.query_catalog.query_flybase_orthologs_per_species(
                    species1, species2, self.ncbi_xref_rdf_property)
                else:
                    query_OMA = self.query_catalog.query_flybase_orthologs_per_species(
                    species1, species2, self.ensembl_xref_rdf_property)
        else:
            if species1 in self.ncbi_species and species2 in self.ncbi_species:
                query_OMA = self.query_catalog.query_orthologs_per_species(
                species1, species2, self.ncbi_xref_rdf_property, self.ncbi_xref_rdf_property)
            else:
                if species1 in self.ncbi_species:
                    query_OMA = self.query_catalog.query_orthologs_per_species(
                    species1, species2, self.ncbi_xref_rdf_property, self.ensembl_xref_rdf_property)
                else:
                    if species2 in self.ncbi_species:
                        query_OMA = self.query_catalog.query_orthologs_per_species(
                        species1, species2, self.ensembl_xref_rdf_property, self.ncbi_xref_rdf_property)
                    else:
                        query_OMA = self.query_catalog.query_orthologs_per_species(
                        species1, species2, self.ensembl_xref_rdf_property, self.ensembl_xref_rdf_property)

        # Treat special case where the gene name does not start with ENS. In OMA 2020, for example,
        # WBGeneXXXXX genes are not considered as an Ensembl cross-reference, hence, it requires a different query.
        if self.special_case_species is not None:
            if species1 in self.special_case_species: #or species2 in self.special_case_species:
                #if species2 in self.special_case_species:
                #    temp = species2
                #    species2 = species1
                #    species1 = temp
                species1_prefix = self.species2prefix.get(species1)
                if species2 in self.species_drosophila:
                    query_OMA = self.query_catalog.query_flybase_orthologs_per_species_filter_gene_prefix(
                        species1, species2, species1_prefix)
                else:
                    if species2 in self.ncbi_species:
                        query_OMA = self.query_catalog.query_orthologs_per_species_filter_gene_prefix(
                            species1, species2, self.ncbi_xref_rdf_property, species1_prefix)
                    else:
                        if species2 in self.special_case_species:
                            query_OMA = self.query_catalog.query_orthologs_per_species_filter_gene_prefixes(
                                species1, species2, species1_prefix, self.species2prefix.get(species2))
                        else:
                            query_OMA = self.query_catalog.query_orthologs_per_species_filter_gene_prefix(
                                species1, species2, self.ensembl_xref_rdf_property, species1_prefix)
            elif species2 in self.special_case_species:
                    species2_prefix = self.species2prefix.get(species2)
                    if species1 in self.species_drosophila:
                        query_OMA = self.query_catalog.query_flybase_orthologs_per_species_filter_gene_prefix(
                            species2, species1, species2_prefix, True)
                    else:
                        if species1 in self.ncbi_species:
                            query_OMA = self.query_catalog.query_orthologs_per_species_filter_gene_prefix(
                                species2, species1, self.ncbi_xref_rdf_property, species2_prefix, True)
                        else:
                            if species1 in self.special_case_species:
                                query_OMA = self.query_catalog.query_orthologs_per_species_filter_gene_prefixes(
                                    species2, species1, species2_prefix, self.species2prefix.get(species1), True)
                            else:
                                query_OMA = self.query_catalog.query_orthologs_per_species_filter_gene_prefix(
                                    species2, species1, self.ensembl_xref_rdf_property, species2_prefix, True)
        return query_OMA

    # Generate the paralogous relations per pair of species and their direct speciation event before the duplication
    # as a NCBI taxon name and identifier.
    def get_query_OMA_paralogs(self, species1: str, species2: str) -> str:
        """Get the corresponding paralogy related OMA SPARQL query from the self.query_catalog
         for a given pair of species.

        :param species1: A first species of reference for paralogy relationship.
        :param species2: A second species of reference for paralogy relationship.

        :return: str : the SPARQL query using the OMA data schema.
        """
        if species1 in self.species_drosophila or species2 in self.species_drosophila:
            if species1 in self.species_drosophila and species2 in self.species_drosophila:
                query_OMA = self.query_catalog.query_flybase_species_paralogs(species1, species2)
            else:
                #if species2 in self.species_drosophila:
                #    temp = species2
                #    species2 = species1
                #    species1 = temp
                if species1 in self.species_drosophila:
                    if species2 in self.ncbi_species:
                        query_OMA = self.query_catalog.query_flybase_paralogs_per_species(species1, species2,
                                                                                      self.ncbi_xref_rdf_property)
                    else:
                        query_OMA = self.query_catalog.query_flybase_paralogs_per_species(species1, species2,
                                                                                      self.ensembl_xref_rdf_property)
                elif species2 in self.species_drosophila:
                    if species1 in self.ncbi_species:
                        query_OMA = self.query_catalog.query_flybase_paralogs_per_species(species2, species1,
                                                                                      self.ncbi_xref_rdf_property)
                    else:
                        query_OMA = self.query_catalog.query_flybase_paralogs_per_species(species2, species1,
                                                                                      self.ensembl_xref_rdf_property)
        else:
            if species1 in self.ncbi_species and species2 in self.ncbi_species:
                query_OMA = self.query_catalog.query_paralogs_per_species(
                    species1, species2, self.ncbi_xref_rdf_property, self.ncbi_xref_rdf_property)
            else:
                if species1 in self.ncbi_species:
                    query_OMA = self.query_catalog.query_paralogs_per_species(
                        species1, species2, self.ncbi_xref_rdf_property, self.ensembl_xref_rdf_property)
                else:
                    if species2 in self.ncbi_species:
                        query_OMA = self.query_catalog.query_paralogs_per_species(
                            species1, species2, self.ensembl_xref_rdf_property, self.ncbi_xref_rdf_property)
                    else:
                        query_OMA = self.query_catalog.query_paralogs_per_species(
                            species1, species2, self.ensembl_xref_rdf_property, self.ensembl_xref_rdf_property)
        # Treat special case where the gene name does not start with ENS. In OMA 2020, for example,
        # WBGeneXXXXX genes are not considered as an Ensembl cross-reference, hence, it requires a different query.
        if self.special_case_species is not None:
            if species1 in self.special_case_species: #or species2 in self.special_case_species:
                #if species2 in self.special_case_species:
                #    temp = species2
                #    species2 = species1
                #   species1 = temp
                species1_prefix = self.species2prefix.get(species1)
                if species2 in self.species_drosophila:
                    query_OMA = self.query_catalog.query_flybase_paralogs_per_species_filter_gene_prefix(
                        species1, species2, species1_prefix)
                else:
                    if species2 in self.ncbi_species:
                        query_OMA = self.query_catalog.query_paralogs_per_species_filter_gene_prefix(
                            species1, species2, self.ncbi_xref_rdf_property, species1_prefix)
                    else:
                        if species2 in self.special_case_species:
                            query_OMA = self.query_catalog.query_paralogs_per_species_filter_gene_prefixes(
                                species1, species2, species1_prefix, self.species2prefix.get(species2))
                        else:
                            query_OMA = self.query_catalog.query_paralogs_per_species_filter_gene_prefix(
                                species1, species2, self.ensembl_xref_rdf_property, species1_prefix)
            elif species2 in self.special_case_species:
                species2_prefix = self.species2prefix.get(species2)
                if species1 in self.species_drosophila:
                    query_OMA = self.query_catalog.query_flybase_paralogs_per_species_filter_gene_prefix(
                        species2, species1, species2_prefix, True)
                else:
                    if species1 in self.ncbi_species:
                        query_OMA = self.query_catalog.query_paralogs_per_species_filter_gene_prefix(
                            species2, species1, self.ncbi_xref_rdf_property, species2_prefix, True)
                    else:
                        if species1 in self.special_case_species:
                            query_OMA = self.query_catalog.query_paralogs_per_species_filter_gene_prefixes(
                                species2, species1, species2_prefix, self.species2prefix.get(species1), True)
                        else:
                            query_OMA = self.query_catalog.query_paralogs_per_species_filter_gene_prefix(
                                species2, species1, self.ensembl_xref_rdf_property, species2_prefix, True)
        return query_OMA


