from datetime import datetime
import errno
import logging, os
import fcntl
import time

from typing.io import IO

from util import Util
from SPARQLWrapper import JSON, SPARQLWrapper


class HomologsFileWriter:
    """Class to write the output files that contain the orthologs and paralogs.
    """
    def __init__(self, species1: str, species2: str, taxon_namespace: str,
                 file_directory: str, file_name_prefix: str, file_index: int,
                 output_header: str) -> object:
        """The HomologsFileWriter class constructor.

        :param species1: A first species of reference for orthology or paralogy relationship.
        :param species2: A second species of reference for orthology or paralogy relationship.
        :param taxon_namespace: The taxon namespace prefix used for species1 and species2.
        :param file_directory:  The file directory path where the file will be saved.
        :param file_name_prefix: The file name prefix.
        :param file_index: The file creation index.
        :param output_header: The header of the created file.
        """
        if species1.startswith(taxon_namespace):
            self.species1 = species1
        else:
            self.species1 = taxon_namespace + species1
        if species2.startswith(taxon_namespace):
            self.species2 = species2
        else:
            self.species2 = taxon_namespace + species2
        self.taxon_namespace = taxon_namespace
        self.file_directory = file_directory
        self.file_name_prefix = file_name_prefix
        self.file_index = file_index
        self.output = output_header
        if self.output is None:
            raise ValueError("File header is not well defined.")

    def _create_homologs_output_file(self, sparql_wrapper_endpoint: SPARQLWrapper, sparql_query: str,
                                     tmp_dir: str = '', drop_duplicates: bool = False, mapper: dict = None,
                                     species_mapper: list = None, query_projection_to_map: list = None):
        """Create text files with either orthology or paralogy relationships.

        :param sparql_wrapper_endpoint: The SPARQLWrapper object that contains info about the SPARQL endpoint.
        :param sparql_query: The query to be executed by the SPARQL endpoint via the sparql_wrapper_endpoint.
        :param tmp_dir: The directory path where temporary files containing the processed species are saved. (optional)
        """
        # set the query to be executed against the OMA endpoint and set the return format to JSON
        sparql_wrapper_endpoint.setQuery(sparql_query)
        sparql_wrapper_endpoint.setReturnFormat(JSON)
        logging.debug("\n{} - Executing sparql query ...\n".format(datetime.now()) + sparql_query)
        results_OMA = sparql_wrapper_endpoint.query().convert()
        logging.debug("{} - Finished sparql query ...".format(datetime.now()))
        if results_OMA["results"]["bindings"] != []:
            #header = results_OMA["results"]["bindings"][0].keys()
            taxid1 = self.species1.split(self.taxon_namespace,1)[1]
            taxid2 = self.species2.split(self.taxon_namespace,1)[1]
            file = os.path.join(self.file_directory,  self.file_name_prefix + "_" + taxid1 + '-' + taxid2 + ".csv")
            logging.debug("{} - ".format(datetime.now()) + str(self.file_index) +' Writing file ...\n' + file)
            if species_mapper is not None and (self.species1 in species_mapper or self.species2 in species_mapper):
                Util.rewrite_results_csv(self.output, file, results_OMA, query_projection_to_map,
                                         mapper, drop_duplicates=drop_duplicates)
            else:
                Util.rewrite_results_csv(self.output, file, results_OMA, drop_duplicates=drop_duplicates)
            self.__write_pairs_tmp_file(os.path.join(tmp_dir, self.file_name_prefix + "_pairs_with_results.tmp"),
                                   [self.species1,self.species2].__str__().replace("'", "\"") + ",")
            logging.info("{} - Saved file".format(datetime.now()) + file)
        else:
            self.__write_pairs_tmp_file(os.path.join(tmp_dir,  self.file_name_prefix + "_pairs_no_results.tmp"),
                                        [self.species1, self.species2].__str__().replace("'", "\"") + ",")
            logging.warning('No results returned for query:\n ' + sparql_query)

    def __write_pairs_tmp_file(self, file_path: str, data):
        with open(file_path, "a+") as pairs_file:
            while True:
                try:
                    fcntl.flock(pairs_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    pairs_file.write(data)
                    break
                except IOError as e:
                    # raise on unrelated IOErrors
                    if e.errno != errno.EAGAIN:
                        raise
                    else:
                        time.sleep(0.1)