import multiprocessing
import os
import time
from os import path
from SPARQLWrapper import SPARQLWrapper
import logging
from DataReader import SpeciesConfigDataReader, GeneIDMapperDataReader
from PairwiseHomology import PairwiseHomologyOMA
from FileWriter import HomologsFileWriter
from util import Util
import sys, getopt
from pathlib import Path
from DataReader import ConfigFileError, ConfigFileNotFoundError
from distutils.util import strtobool
import ast
from multiprocessing import Process, Value


ORTHOLOGY_PREFIX = "orthologs"
PARALOGY_PREFIX = "paralogs"


class ProgressBar(Process):

    def __init__(self, total_max_files: int, counter,  stop_error_msg, *args, **kwargs) -> object:
        super(ProgressBar, self).__init__(*args, **kwargs)
        self.total_max_files = total_max_files
        self.counter = counter
        self.stop_error = stop_error_msg

    def run(self):
        count = self.counter.value
        while count < self.total_max_files and not self.stop_error.is_set():
            Util.printProgressBar(count, self.total_max_files, prefix='Progress:', suffix='Complete', length=50)
            time.sleep(60)
            count = self.counter.value
        if not self.stop_error.is_set():
            Util.printProgressBar(self.total_max_files, self.total_max_files,
                              prefix='Progress:', suffix='Complete', length=50)
            print("\nEnded with success.")
        else:
            print("Some process raised errors. For further information, " +
                  "check individual process logs.")
            print("You can try to rerun this application with the same template.properties file and parameters." +
                  " The application will try to continue from where it stopped. Activate the debug mode in the"
                  + " template.properties for more information about the error (log = debug).")




class GenerateHomologsParallel(Process):
    # Constants declarations

    def __init__(self, path_config_file, counter, tmp_dir, stop_main, *args, **kwargs):
        super(GenerateHomologsParallel, self).__init__(*args, **kwargs)
        self.config_file = path_config_file
        self.ORTHOLOGY_PREFIX = ORTHOLOGY_PREFIX
        self.PARALOGY_PREFIX = PARALOGY_PREFIX
        self.counter = counter
        self.tmp_dir = tmp_dir
        self.stop_main = stop_main

    def run(self):
        """
        The application main function. The runtime exceptions are treated here.
        Raised errors are saved in the log file (defined in the LOG_FILE constant).
        The log file is at the same directory of the execution of this app.
        """
        try:
            conf_dir = self.config_file
            tmp_dir = self.tmp_dir
            config_file_path = os.path.join(conf_dir, "config.properties")
            LOG_FILE = os.path.join(conf_dir, 'generate_homologs_app.log')
            # Load config.properties file and set the constants with the config file parameters values.
            config = SpeciesConfigDataReader.load_config_file(config_file_path)
            FILE_DIRECTORY = config['OUTPUT']['FILE_DIRECTORY_OUTPUT']
            OMA_SPARQL_ENDPOINT = str(config['DEFAULT']['OMA_SPARQL_ENDPOINT'])
            NCBI_GENE2ENSEMBL_URL = str(config['DEFAULT']['ncbi_gene2ensembl_URL'])
            SPECIES_FILE = config['SPECIES']['SPECIES_FILE']
            TAXON_NAMESPACE = str(config['SPECIES']['TAXON_NAMESPACE'])
            NCBI_GENE_CROSSREF_PROPERTY = str(config['RDF_VOCABULARY']['NCBI_GENE_CROSSREF_PROPERTY'])
            ENSEMBL_GENE_CROSSREF_PROPERTY = str(config['RDF_VOCABULARY']['ENSEMBL_GENE_CROSSREF_PROPERTY'])
            DEFAULT_PARAMETERS = config['DEFAULT']
            try:
                IS_REMOVE_PAIRS = strtobool(config['DEFAULT']['remove_species_pairs'])
            except:
                IS_REMOVE_PAIRS = False
            try:
                SPECIES_MAPPING = config['SPECIES_MAPPING']
                map_species = SPECIES_MAPPING.items()
            except:
                map_species = []
            try:
                NCBI_SPECIES_FILE = config['SPECIES']['NCBI_SPECIES_FILE']
            except:
                NCBI_SPECIES_FILE = None
            try:
                FLYBASE_SPECIES_FILE = config['SPECIES']['FLYBASE_SPECIES_FILE']
            except:
                FLYBASE_SPECIES_FILE = None
            try:
                DYNAMIC_VARIABLES = config['DYNAMIC_VARIABLES']
            except:
                config['DYNAMIC_VARIABLES'] = {'START_INDEX_ORTHOLOG': '0',
                                               'START_INDEX_PARALOG': '0'}
            try:
                ORTHOLOG_INDEX = int(config['DYNAMIC_VARIABLES']['START_INDEX_ORTHOLOG'])
                if ORTHOLOG_INDEX < 0:
                    ORTHOLOG_INDEX = 0
                    config['DYNAMIC_VARIABLES']['START_INDEX_ORTHOLOG'] = '0'
            except:
                config['DYNAMIC_VARIABLES']['START_INDEX_ORTHOLOG'] = '0'
                ORTHOLOG_INDEX = 0
            # the maximum number of species pairs to consider, if stated 0 or none in the config files the default set to 0.
            try:
                MAX_INDEX = int(config['DYNAMIC_VARIABLES']['MAX_INDEX'])
                if MAX_INDEX < 0:
                    MAX_INDEX = 0
            except:
                MAX_INDEX = 0
            try:
                PARALOG_INDEX = int(config['DYNAMIC_VARIABLES']['START_INDEX_PARALOG'])
                if PARALOG_INDEX < 0:
                    PARALOG_INDEX = 0
                    config['DYNAMIC_VARIABLES']['START_INDEX_PARALOG'] = '0'
            except:
                config['DYNAMIC_VARIABLES']['START_INDEX_PARALOG'] = '0'
                PARALOG_INDEX = 0
            # State the Debug level, if it is not defined in the config file, the info level is set.
            try:
                DEBUG = str(config['DEFAULT']['log'])
            except:
                DEBUG = 'info'
                config['DEFAULT']['log'] = 'info'
            if DEBUG.upper().__eq__("DEBUG"):
                logging.basicConfig(filename=LOG_FILE, level=logging.DEBUG)
            else:
                logging.basicConfig(filename=LOG_FILE, level=logging.INFO)
            logging.info('Started logging mode: ' + DEBUG)
            try:
                SPECIES_TO_GENE_ID_PREFIX = ast.literal_eval(config['SPECIES']['SPECIES_TO_GENE_ID_PREFIX'])
            except:
                config['SPECIES']['SPECIES_TO_GENE_ID_PREFIX'] = "{6239: 'WBGene', 7740: 'BL'}"
                SPECIES_TO_GENE_ID_PREFIX = {6239: 'WBGene', 7740: 'BL'}
                logging.warning(
                    "SPECIES_TO_GENE_ID_PREFIX parameter was not defined or wrong syntax.\n The parameter value"
                    + " must be a python dictionary.")
                logging.warning(
                    "A default value was defined:\nSPECIES_TO_GENE_ID_PREFIX = {6239: 'WBGene', 7740: 'BL'}")
            try:
                SPECIES_NCBI_TO_ENSEMBL_IDS = ast.literal_eval(config['SPECIES']['species_ncbi_to_ensembl'])
            except:
                config['SPECIES']['species_ncbi_to_ensembl'] = "[105023]"
                SPECIES_NCBI_TO_ENSEMBL_IDS = [105023]
                logging.warning(
                    "'species_ncbi_to_ensembl' parameter was not defined or wrong syntax.\n The parameter value"
                    + " must be a list e.g. [1,23,5].")
                logging.warning(
                    "A default value was defined:\nspecies_ncbi_to_ensembl = [105023]")
            logging.info('config.properties file loaded')
            # Create directory if does not exist
            Path(str(FILE_DIRECTORY)).mkdir(parents=True, exist_ok=True)
            # the endpoints must be defined as wrappers for executing SPARQL queries
            sparql_OMA = SPARQLWrapper(OMA_SPARQL_ENDPOINT)
            logging.info('Set SPARQL endpoints')
            # Load species file that contains all species
            species_bgee_list = SpeciesConfigDataReader.get_species_from_file(SPECIES_FILE, TAXON_NAMESPACE)
            # Load species file that contains species from NCBI
            ncbi_species = []
            if not NCBI_SPECIES_FILE is None:
                ncbi_species = SpeciesConfigDataReader.get_species_from_file(NCBI_SPECIES_FILE, TAXON_NAMESPACE)
            # Load species file that contains species from FlyBase
            flybase_species = []
            if not FLYBASE_SPECIES_FILE is None:
                flybase_species = SpeciesConfigDataReader.get_species_from_file(FLYBASE_SPECIES_FILE, TAXON_NAMESPACE)
            logging.info('Loaded species related files')
            # If necessary, map species tax ids in species lists based on mappings defined
            # in the config.properties file
            for item in map_species:
                if not item in DEFAULT_PARAMETERS.items():
                    Util.replaces_item_string_list(species_bgee_list, TAXON_NAMESPACE + item[0],
                                                   TAXON_NAMESPACE + item[1])
                    Util.replaces_item_string_list(ncbi_species, TAXON_NAMESPACE + item[0], TAXON_NAMESPACE + item[1])
                    Util.replaces_item_string_list(flybase_species, TAXON_NAMESPACE + item[0],
                                                   TAXON_NAMESPACE + item[1])
            # Create the species pair list excluding symetric pairs and pairs with the same species
            values_species = Util.get_asymmetric_pairs_from_list(species_bgee_list, False)
            if IS_REMOVE_PAIRS and ORTHOLOG_INDEX < len(values_species):
                Util.remove_pairs(os.path.join(tmp_dir, self.ORTHOLOGY_PREFIX + "_pairs_with_results.tmp"),
                                  values_species)
               # Start index is set to 0 since the size of values_species changed with the removal of some species pairs
               # ORTHOLOG_INDEX = 0

            # generate pairwise orthology information files from OMA database
            pairwise_homologs = PairwiseHomologyOMA(ncbi_species, flybase_species, TAXON_NAMESPACE,
                                                    NCBI_GENE_CROSSREF_PROPERTY, ENSEMBL_GENE_CROSSREF_PROPERTY,
                                                    SPECIES_TO_GENE_ID_PREFIX)
            geneID_mapper = GeneIDMapperDataReader(sparql_OMA, TAXON_NAMESPACE, SPECIES_TO_GENE_ID_PREFIX,
                                                   NCBI_GENE2ENSEMBL_URL, SPECIES_NCBI_TO_ENSEMBL_IDS)
            geneID_mapper_dict = geneID_mapper.get_gene_id_mapper()
            geneID_mapper_species = geneID_mapper.get_mapper_species()
            species_pairs_size = len(values_species)
            if ORTHOLOG_INDEX > species_pairs_size:
                ORTHOLOG_INDEX = species_pairs_size
            if MAX_INDEX == 0:
                MAX_INDEX = species_pairs_size
            count = ORTHOLOG_INDEX + 1
            for (species1, species2) in values_species[ORTHOLOG_INDEX:MAX_INDEX]:
                query_OMA = pairwise_homologs.get_query_OMA_orthologs(species1, species2)
                homolog_file_writer = HomologsFileWriter(species1, species2, TAXON_NAMESPACE,
                                                         FILE_DIRECTORY, self.ORTHOLOGY_PREFIX, count,
                                                         "gene1,gene2,tax_level,tax_level_id")
                if species1 in geneID_mapper.species_url_ncbi_to_ensembl or\
                    species2 in geneID_mapper.species_url_ncbi_to_ensembl:
                    homolog_file_writer._create_homologs_output_file(sparql_OMA, query_OMA, tmp_dir, False,
                                                                     geneID_mapper_dict, geneID_mapper_species,
                                                                     ["gene1", "gene2"])
                else:
                    homolog_file_writer._create_homologs_output_file(sparql_OMA, query_OMA, tmp_dir)
                config['DYNAMIC_VARIABLES']['START_INDEX_ORTHOLOG'] = str(count) 
                with open(config_file_path, 'w') as configfile:
                    config.write(configfile)
                count = count + 1
                with self.counter.get_lock():
                    self.counter.value += 1
            # Create the species pair list excluding symetric pairs, e.g. (a,b) = (b,a), thus only one of them is kept
            values_species = Util.get_asymmetric_pairs_from_list(species_bgee_list, True)
            if IS_REMOVE_PAIRS and PARALOG_INDEX < len(values_species):
                Util.remove_pairs(os.path.join(tmp_dir, self.PARALOGY_PREFIX + "_pairs_with_results.tmp"),
                                  values_species)
                # Start index is set to 0 since the size of values_species changed with the removal of some species pairs
                #PARALOG_INDEX = 0
            # generate pairwise paralogy information files from OMA database
            species_pairs_size = len(values_species)
            if PARALOG_INDEX > species_pairs_size:
                PARALOG_INDEX = species_pairs_size
            if MAX_INDEX == 0:
                MAX_INDEX = species_pairs_size
            count = PARALOG_INDEX + 1
            for (species1, species2) in values_species[PARALOG_INDEX:MAX_INDEX]:
                query_OMA = pairwise_homologs.get_query_OMA_paralogs(species1, species2)
                homolog_file_writer = HomologsFileWriter(species1, species2, TAXON_NAMESPACE,
                                                         FILE_DIRECTORY, self.PARALOGY_PREFIX, count,
                                                         "gene1,gene2,tax_level,tax_level_id")
                homolog_file_writer._create_homologs_output_file(sparql_OMA, query_OMA, tmp_dir, False,
                                                                 geneID_mapper_dict, geneID_mapper_species,
                                                                 ["gene1", "gene2"])
                config['DYNAMIC_VARIABLES']['START_INDEX_PARALOG'] = str(count)
                with open(config_file_path, 'w') as configfile:
                    config.write(configfile)
                count = count + 1
                with self.counter.get_lock():
                    self.counter.value += 1
            # reinitialize to zero the array index of species pairs, because all pairs were processed with success
           # config['DYNAMIC_VARIABLES']['START_INDEX_PARALOG'] = '0'
           # config['DYNAMIC_VARIABLES']['START_INDEX_ORTHOLOG'] = '0'
            with open(config_file_path, 'w') as configfile:
                config.write(configfile)
            logging.info('Logging finished without errors\n')
        except ConfigFileError as e:
            print()
            print(e.message)
            logging.debug("Configuration file error: " + e.message)
            logging.info('Logging finished with errors\n')
            self.stop_main.set()
        except ConfigFileNotFoundError as e:
            print()
            print(e.message)
            logging.debug("Configuration file error: " + e.message)
            logging.info('Logging finished with errors\n')
            self.stop_main.set()
        except ValueError as e:
            logging.debug("Value error: " + str(e.args[0]))
            logging.info('Logging finished with errors\n')
            print("\nUnexpected error happened, see .log file for further details at\n{}".format(conf_dir))
            self.stop_main.set()
        except FileNotFoundError as e:
            logging.debug("File not found error: {}".format(e.args[0]))
            logging.info('Logging finished with errors\n')
            print("\nUnexpected error happened, see .log file for further details at\n{}".format(conf_dir))
            self.stop_main.set()
        except:
            logging.debug("Unexpected error: " + str(sys.exc_info()[0]) + "\n" + str(sys.exc_info()[1]))
            logging.info('Logging finished with errors\n')
            print("\nUnexpected error happened, see .log file for further details at\n{}".format(conf_dir))
            self.stop_main.set()


def getArgs(argv):
    inputdir = ''
    no_ortholog = False
    no_paralog = False
    process_number = 3
    try:
        opts, args = getopt.getopt(argv, "hopc:n:", ["configdir=", "process="])
    except getopt.GetoptError:
        print('GenerateHomologsParallelApp.py -c <configuration (config.properties) and temp file directory>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('GenerateHomologsParallelApp.py -c <configuration (config.properties) and temp file directory>')
            print('-p : only output paralogs <optional>')
            print('-o : only output orthologs <optional>')
            print('-n or --process : number of processes <the default value is -n 3>')
            sys.exit()
        elif opt in ("-c", "--configdir"):
            inputdir = arg
        elif opt in ("-o"):
            no_paralog = True
        elif opt in ("-p"):
            no_ortholog = True
        elif opt in ("-n", "--process"):
            try:
                process_number = int(arg)
                if process_number > 4 or process_number < 1:
                    process_number = 3
            except:
                print("Number of processes (-n) must be an integer.")

    print('Configuration and temporary file directory is ', inputdir)
    print('In case no output directory for the created paralog and orhtolog files is provided in ',
          'template.properties file, <file_directory_output> parameter, ', inputdir, ' is used as default')
    return {'dir': inputdir, 'no_paralog': no_paralog,'no_ortholog': no_ortholog, 'process_number':process_number }

def main(argv):
    try:
        multiprocessing.set_start_method('spawn')
        stop_main = multiprocessing.Event()
        stop_error_msg = multiprocessing.Event()
        counter = Value('i', 0)
        args_dict = getArgs(argv)
        dir_path = args_dict['dir']
        Path(dir_path).mkdir(parents=True, exist_ok=True)
        config = SpeciesConfigDataReader.load_config_file(
            os.path.join(os.path.dirname(__file__), 'template.properties'))
        FILE_DIRECTORY = config['OUTPUT']['FILE_DIRECTORY_OUTPUT']
        if FILE_DIRECTORY is None or len(str(FILE_DIRECTORY)) == 0:
            config['OUTPUT']['FILE_DIRECTORY_OUTPUT'] = os.path.join(dir_path, "gene_homology_data")
        with open(config['SPECIES']['SPECIES_FILE'], 'r') as species_file:
            lines = species_file.readlines()
            species_count = list(filter("".__ne__, lines)).__len__() - 1
        number_orth_files = (species_count * (species_count - 1)) / 2.0
        number_paralog_files = (species_count * (species_count + 1)) / 2.0
        total_max_files = number_orth_files + number_paralog_files
        procs = []
        chunk_size = number_paralog_files / args_dict['process_number']
        previous_max = 0
        proc_id_list = range(1, args_dict['process_number'] + 1)
        if config['DEFAULT']['remove_species_pairs'] != 'True' \
                and config['DEFAULT']['remove_species_pairs'] != 'False':
            print("here")
            raise ValueError("'remove_species_pairs' value error in the configuration file. "
                             + "It must be either True and False.")
        if config['DEFAULT']['remove_species_pairs'] == 'True':
            orth_tmp_file_path = os.path.join(dir_path, ORTHOLOGY_PREFIX + "_pairs_with_results.tmp")
            para_tmp_file_path = os.path.join(dir_path, PARALOGY_PREFIX + "_pairs_with_results.tmp")
            para_pairs_remove_list = []
            orth_pairs_remove_list = []
            if path.exists(orth_tmp_file_path):
                with open(orth_tmp_file_path, 'r') as tmp_file:
                    orth_pairs_remove_list = Util.read_tmp_file(tmp_file)
            if path.exists(para_tmp_file_path):
                with open(para_tmp_file_path, 'r') as tmp_file:
                    para_pairs_remove_list = Util.read_tmp_file(tmp_file)
            if len(para_pairs_remove_list) or len(orth_pairs_remove_list):
                total_max_files = total_max_files - len(orth_pairs_remove_list) - len(para_pairs_remove_list)
                number_orth_files = number_orth_files - len(orth_pairs_remove_list)
                number_paralog_files = number_paralog_files - len(para_pairs_remove_list)
                if number_orth_files > number_paralog_files:
                    number_paralog_files = number_orth_files
                chunk_size = number_paralog_files / args_dict['process_number']
        pbar = ProgressBar(total_max_files, counter, stop_error_msg)
        pbar.daemon = True
        pbar.start()
        for process_number in proc_id_list:
            config['DYNAMIC_VARIABLES'] = {'max_index': str(round(chunk_size * process_number))}
            if args_dict['no_ortholog']:
                config['DYNAMIC_VARIABLES']['START_INDEX_ORTHOLOG'] = str(number_orth_files)
            else:
                config['DYNAMIC_VARIABLES']['START_INDEX_ORTHOLOG'] = str(round(chunk_size * (process_number - 1)))
            if args_dict['no_paralog']:
                config['DYNAMIC_VARIABLES']['START_INDEX_PARALOG'] = str(number_paralog_files)
            else:
                config['DYNAMIC_VARIABLES']['START_INDEX_PARALOG'] = str(round(chunk_size * (process_number - 1)))
            config_process_dir = os.path.join(dir_path, "process{}_tmp".format(process_number))
            Path(config_process_dir).mkdir(parents=True, exist_ok=True)
            config_file_path = os.path.join(config_process_dir, 'config.properties')
            if process_number == proc_id_list[proc_id_list.__len__() - 1]:
                config['DYNAMIC_VARIABLES']['max_index'] = str(int(number_paralog_files))
            if not path.exists(config_file_path):
                with open(config_file_path, 'w') as configfile:
                    config.write(configfile)
            else:
                config_exist = SpeciesConfigDataReader.load_config_file(config_file_path)
                try:
                    if config['DEFAULT']['remove_species_pairs'] == 'True':
                        config_exist['DEFAULT']['remove_species_pairs'] = config['DEFAULT']['remove_species_pairs']
                        config_exist['DYNAMIC_VARIABLES']['START_INDEX_PARALOG'] = \
                            config['DYNAMIC_VARIABLES']['START_INDEX_PARALOG']
                        config_exist['DYNAMIC_VARIABLES']['START_INDEX_ORTHOLOG'] = \
                            config['DYNAMIC_VARIABLES']['START_INDEX_ORTHOLOG']
                        config_exist['DYNAMIC_VARIABLES']['MAX_INDEX'] = \
                            config['DYNAMIC_VARIABLES']['MAX_INDEX']
                        with open(config_file_path, 'w') as configfile:
                            config_exist.write(configfile)
                    else:
                        counter.value += int(config_exist['DYNAMIC_VARIABLES']['START_INDEX_ORTHOLOG']) \
                                         + int(config_exist['DYNAMIC_VARIABLES']['START_INDEX_PARALOG']) - (
                                                     2 * previous_max)
                    # Make sure that the last process defines a max index equal the estimated number of paralog files
                    if process_number == proc_id_list[len(proc_id_list) - 1]:
                        if int(config_exist['DYNAMIC_VARIABLES']['MAX_INDEX']) < number_paralog_files:
                            config_exist['DYNAMIC_VARIABLES']['MAX_INDEX'] = str(int(number_paralog_files))
                            with open(config_file_path, 'w') as configfile:
                                config_exist.write(configfile)
                except:
                    raise
            generator = GenerateHomologsParallel(config_process_dir, counter, dir_path, stop_main)
            generator.daemon = True
            generator.start()
            procs.append(generator)
            previous_max = int(chunk_size * process_number)
        for p in procs:
            p.join()
        if stop_main.is_set():
            stop_error_msg.set()
        pbar.join()
    except KeyboardInterrupt:
        print("Quiting the program and terminating all running processes...")
    except ValueError as e:
        print("Value error: " + str(e.args[0]) )
    except:
        print("The application has unexpectedly quit.")
    finally:
        for p in procs:
            if p.is_alive():
                p.terminate()
        if pbar.is_alive():
            pbar.terminate()
    return


if __name__ == '__main__':
    main(sys.argv[1:])



