import configparser
from SPARQLWrapper import SPARQLWrapper, JSON, CSV
import logging
import pandas
import io
from QueryCatalog import QueryCatalogOMA
from util import Util
from configparser import ConfigParser


class SpeciesConfigDataReader:
    """Class to read and load configuration parameters and species lists, for example, from files.
    This class only contains static methods. No internal attributes.
    """

    @staticmethod
    def get_species_from_sparql_endpoint(sparql_wrapper: SPARQLWrapper,
                                         query: str = None, query_projection: str = "species") -> list:
        """Get species from a SPARQL endpoint.

        :param sparql_wrapper: (SPARQLWrapper)
            Wrapper around an online access to a SPARQL Web entry point.
        :param query : (str, optional)
            A SPARQL query for getting species from a SPARQL endpoint. The default query retrieves species by using
            the UniProt core RDF vocabulary. Default:
            PREFIX up: <http://purl.uniprot.org/core/>
                SELECT ?species {
                ?species a up:Taxon .
                ?species up:rank up:Species .
                }
        :param query_projection : (str)
            The query projection variable name used in the query. (default is 'species')
        :return: (list)
            a species list
        """
        if query is None:
            query = """PREFIX up: <http://purl.uniprot.org/core/>
                SELECT ?species {
                ?species a up:Taxon .
                ?species up:rank up:Species .
                }
                """
        sparql_bgee = sparql_wrapper
        sparql_bgee.setQuery(query)
        sparql_bgee.setReturnFormat(JSON)
        results_bgee = sparql_bgee.query().convert()
        species_bgee_list = list(Util.convert_result_to_key_value(query_projection, results_bgee))
        return species_bgee_list

    @staticmethod
    def get_drosophila_species_from_sparql_endpoint(sparql_wrapper: SPARQLWrapper) -> object:
        """Get all drosophila like species where the NCBI id starts with 72 from a SPARQL endpoint.

         The SPARQL endpoint must use the UniProt core ontology for describing species.
        :param sparql_wrapper: (SPARQLWrapper)
          Wrapper around an online access to a SPARQL Web entry point.
        :return:
          list: a list of drosophila species
             """
        query = """
    PREFIX up: <http://purl.uniprot.org/core/>
    SELECT ?species {
    ?species a up:Taxon .
    ?species up:rank up:Species .
    filter strstarts(str(?species),'http://purl.uniprot.org/taxonomy/72')
    }
    """
        return SpeciesConfigDataReader.get_species_from_sparql_endpoint(sparql_wrapper,query)

    @staticmethod
    def get_species_from_file(file_path: str, taxon_namespace: str) -> list:
        """Get species list from a text file and add a prefix to it.

        Each species is define in a new line and the first line is not considered (i.e. the column name).

        :param file_path: the file path where the species list is stored. For each new line a species is defined, the
            1st line is not considered (i.e. the column name).
        :param taxon_namespace: it is the prefix to be added for each species from file_path
        :return:
            list : a list of species
        """
        file = open(file_path, 'r')
        lines = file.readlines()
        species = []
        count = 0
        for line in lines:
            line_strip = line.strip()
            if count != 0 and line_strip != '':
                species.append(taxon_namespace + line_strip)
            count = count + 1
        if not species.__len__():
            logging.warning("No species defined in " + file_path)
        file.close()
        return species

    @staticmethod
    def load_config_file(file_path: str = 'config.properties',
                         minimum_parameters_set: list = None) -> ConfigParser:
        """Load the configuration file.

        :param file_path: the file path of the configuration file (default is 'config.properties' in the current
         directory).
        :param minimum_parameters_set: the minimum parameters that the configuration file must contains. If None,
        the defaut is the following list:
            ['file_directory_output',
                                      'oma_sparql_endpoint',
                                      'ncbi_gene_crossref_property',
                                      'ensembl_gene_crossref_property',
                                      'species_file',
                                      'taxon_namespace']

        :return:
            ConfigParser : object built based on file at file_path.
        """
        if minimum_parameters_set is None:
            minimum_parameters_set = ['file_directory_output',
                                      'oma_sparql_endpoint',
                                      'ncbi_gene_crossref_property',
                                      'ensembl_gene_crossref_property',
                                      'species_file',
                                      'taxon_namespace']
        config = configparser.ConfigParser()
        try:
            f = open(file_path)
            f.close()
        except FileNotFoundError:
            raise ConfigFileNotFoundError(file_path + ' file was not found in the current directory.',
                                    "Run this application in the same directory as the file: " + file_path)
        config.read(file_path)
        all_parameters = []
        throw_error = False
        msg_error = 'One or more mandatory parameters are missing in the ' + file_path + ' file, see below:\n'
        for section in config.items():
            for parameter in section.__getitem__(1).items().__iter__():
                all_parameters.append(parameter[0])
        for parameter in minimum_parameters_set:
            if not all_parameters.__contains__(parameter):
                msg_error = msg_error + " " + parameter + '\n'
                throw_error = True
        if throw_error:
            raise ConfigFileError(file_path, msg_error)
        return config

class GeneIDMapperDataReader:
    """Class that provides a gene identifier mapper as a dictionary.
    This class only contains static methods. No internal attributes.
    """
    def __init__(self, sparql_wrapper: SPARQLWrapper, taxon_namespace: str, species2prefix: list, ncbi_csv_URL: str,
                 species_ncbi_to_ensembl_ids: list = None, query_catalog: QueryCatalogOMA = None) -> object:
        """PairwiseHomologyOMA class constructor


        :param species2prefix: A mapping between the species and a prefix used in their gene identifiers, for example,
         "BL" is prefix for gene identifiers of the Branchiostoma lanceolatum species in Ensembl metazoa database.
        :param query_catalog: the Query catalog object that contains the queries to be executed in the sparql endpoint.
        :param species_ncbi_to_ensembl_ids: The species to be considered in order to map the NCBI gene ids
         to Ensembl gene ids.
        """
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
        self.species_url_ncbi_to_ensembl = []
        if species_ncbi_to_ensembl_ids is not None:
            for species_id in species_ncbi_to_ensembl_ids:
                self.species_url_ncbi_to_ensembl.append(taxon_namespace + str(species_id))
        else:
            species_ncbi_to_ensembl_ids = []
        self.__set_gene_id_mapper(sparql_wrapper, ncbi_csv_URL, species_ncbi_to_ensembl_ids)

    def get_gene_id_mapper(self):
        return self.mapper

    def get_mapper_species(self):
        result = []
        if self.special_case_species is not None:
            result.extend(self.special_case_species)
        if self.species_url_ncbi_to_ensembl is not None:
            result.extend(self.species_url_ncbi_to_ensembl)
        return result

    def __set_gene_id_mapper(self, sparql_wrapper_endpoint: SPARQLWrapper, ncbi_csv_url: str, ncbi_species_ids: list):
        self.mapper = {}
        queries = []
        for species in self.species2prefix:
            queries.append(self.query_catalog.get_mapper_OMA_IRI_to_prefixed_id(species, self.species2prefix[species]))
        for query_oma in queries:
            sparql_wrapper_endpoint.setQuery(query_oma)
            sparql_wrapper_endpoint.setReturnFormat(CSV)
            results_oma = sparql_wrapper_endpoint.query().response
            csv_result_pd = pandas.read_csv(results_oma, sep=',', index_col=0, header=None, squeeze=True)
            csv_result_dict = csv_result_pd.to_dict()
            self.mapper.update(csv_result_dict)
        for species in ncbi_species_ids:
            csv_file = pandas.read_csv(ncbi_csv_url, sep='\t')
            csv_file.columns = csv_file.columns.str.replace('#', '')
            df_taxid = csv_file.query('tax_id == {}'.format(species))
            df_taxid = df_taxid.filter(items=["GeneID", "Ensembl_gene_identifier"])
            df_taxid = df_taxid.drop_duplicates()
            df_taxid = df_taxid.groupby("GeneID")['Ensembl_gene_identifier'].apply(list)
            # Recreate a dataframe by removing header and dataframe internal index (indexing the 1st column)
            # in this way, to_dict() creates a dictionary where the key is the 1st column and the value the second one.
            #df_taxid = pandas.read_csv(io.StringIO(u"" + df_taxid.to_csv(index=False, header=None)),
            #                              sep=',', header=None, index_col=0, squeeze=True)
            self.mapper.update(df_taxid.to_dict())
        return self.mapper


class ConfigFileError(Exception):
    """Exception raised for errors in the config.properties input.

    :param
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """
    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


class ConfigFileNotFoundError(FileNotFoundError):
    """Exception raised for config.properties file not found.

    :param
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message
