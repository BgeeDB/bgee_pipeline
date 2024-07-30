## Generate Homologs App
This application generate orthologs and paralogs of a given species list. They are extracted from the OMA
hierarchical orthologous groups (HOGs) stored in a triple store compliant with the SPARQL query language. Only the genes
of the given species are considered. 

### Starting the app
First, set up the make file properties [makefile.properties](makefile.properties).
If needed, refer to sections below to know how to edit the CONFIG_TEMPLATE definition in the
[makefile.properties](makefile.properties) file. 
It is a template used to create the template.properties file (for a full example see, 
[config.properties](config.properties)). 

Second, run the make file as below:

```bash
cd makefile/directory/path
make
```

**REMARK 1:** if pipenv is not installed, you will need to run the "install_pipenv" rule from the Makefile to set up 
pipenv. To do so, add "install_pipenv" to the following line "main: config_file  running_app" in the 
[Makefile](Makefile) as shown below: 
```bash
main: config_file install_pipenv running_app
```

It requires python 3.7 or superior version and the SPARQLWrapper library. Dependencies are defined in
[Pipfile](Pipfile).

**REMARK 2:** the execution of the app can take a long time depending on the number of genomes to process. 
To avoid the overload of the public OMA SPARQL endpoint, the app only allows for executing 4 queries in parallel.
From 1 to 4 processes maximum (default is 3), the variable to set it is in the [makefile.properties](makefile.properties) 
is 
```bash
#number of processes, default is 3 maximum is 4
N_PROCESSES = 4
```

### The configuration file
* All input variable for the app are in the [makefile.properties](makefile.properties) file that will generate
  the template.properties file based on the CONFIG_TEMPLATE definition with the following parameters: 
  
* **Required parameters**
    * **file_directory_output** is *the directory where the pairwise homologous genes are saved*, for each pair
      of species a new file is created. If there is any orthology relation between two species a file is created such as
      orthologs_bgee_9606-7217.csv, where 9606 and 7217 are species NCBI ids. In addition, a separated file for paralogy
      relations is also defined such as follows paralogs_bgee_9606-7217.csv. Moreover, in-species paralogs are  
      considered too, e.g. paralogs_bgee_9606-9606.csv.
    * **oma_sparql_endpoint** is the SPARQL endpoint URI of the OMA database.
    * **ncbi_gene_crossref_property** is the property used in OMA to assign the NCBI gene cross-references such as 
      https://www.ncbi.nlm.nih.gov/gene/?term=108708356 where 108708356 is the NCBI gene identifier.
    * **ensembl_gene_crossref_property** is the property used in OMA to assign the Ensembl gene cross-references.
    * **ncbi_gene2ensembl_URL** is the tab separated value (tsv) file (a URL to download the file also works) 
    that contains NCBI to Ensembl gene id mappings. This tsv file must contain at least the following tab separated
    columns: (#tax_id	GeneID	Ensembl_gene_identifier).
    * **species_file** is the file of all species to be considered. This file is a text file where the first line is not
    taken into account by the app. Indeed, the first line is considered as the header. Each species is defined as new
      line where a species is represented with a NCBI identifier (e.g. 9606 for *homo sapiens*).
    * **taxon_namespace** is the taxonomy namespace prefix used in OMA RDF store. Currently, it is the same as UniProt
      and Bgee databases: 
        ```
        [SPECIES]
        taxon_namespace = http://purl.uniprot.org/taxonomy/
        ```
* **Optional parameters**
  * **log** is the log level mode (i.e. log=info or log=debug), if not defined the "info" mode is considered, 
      the log file (i.e. *generate_homologs_app.log*) is saved in the same directory where this app is running.
  * **ncbi_species_file** is the file of species where the genomes considered are originally from NCBI.
  * **flybase_species_file** is the file of species where the genomes considered are originally from FlyBase.
  * **species_ncbi_to_ensembl** is a list of NCBI species ids we want to map the output NCBI gene ids to Ensembl gene
    ids (e.g. \[105023\]). This will only work if the species (e.g. 105023) is in **ncbi_species_file** and there 
    is a corresponding mapping between the NCBI gene id to Ensembl id in the mapping file: **ncbi_gene2ensembl_URL**.
  * **species_to_gene_id_prefix** is a dictionary with a key-value structure where the key is a species that can be
    represented as a [NCBI taxonomy identifier](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi) and the value
    corresponds to the prefix used in the gene identifiers for this species. The syntax used is the same as the 
    [Python dictionary syntax](https://www.w3schools.com/python/python_dictionaries.asp). The default value for this
    parameter is defined below. We must use this feature for species where the Ensembl gene identifiers 
    do not start with the prefixes: 'ENS' and 'FBgn'.
    ```
    [SPECIES]
    species_to_gene_id_prefix = {6239: 'WBGene', 7740: 'BL'}
    ```
  * **start_index_ortholog** and **start_index_paralog** variables that change during the execution of the app.
      These variables set the progress state of the generation of the orthologs and paralogs (starting from default = 0)
    . Based on them, the app is able to restart the process from where it stopped.
    * To solely process paralogs, we can set **start_index_ortholog** to equal the square of the species amount.
    * To only process orthologs, we can set **start_index_paralog** to equal the square of the species amount.
  * **remove_species_pairs** by default is assigned with *False*. If this parameter is set 'True' and 
     based on the files generated by a previous execution of this application, named
    *orthologs_pairs_with_results.tmp* and *paralogs_pairs_with_results.tmp*,
    it will remove the pairs of species already processed and defined in these files. This is useful for incrementally 
    adding new species and to generate the corresponding paralog and ortholog files without needing to process all 
    species *X* species pairs from the species file. We can also define the pair of species we want to remove in 
    those files by using the JSON syntax. An example is shown below:
    ```
    [ ["http://purl.uniprot.org/taxonomy/9606", "http://purl.uniprot.org/taxonomy/10090"],
      ["http://purl.uniprot.org/taxonomy/9606", "http://purl.uniprot.org/taxonomy/7955"] ]
    ```
    Note that in this case, the **taxon_namespace** must be used as the prefix of species defined in the species files.
    The files must be saved at the same directory as the app will be executed. They must be named 
    *orthologs_pairs_with_results.tmp* and *paralogs_pairs_with_results.tmp*  for the species pairs that will not be 
    considered when generating the orthology and parlogy files, respectively . 
  * The **SPECIES_MAPPING** section is used to define mappings to be taken into account when loading the species files:
  *ncbi_species_file*, *flybase_species_file*, and *species_file*. The mappings are defined as pairs "from = to" 
    such as below. **(optional)**
    ```
    [SPECIES_MAPPING]
    7237 = 46245
    9593 = 9595
    ```
    In this example, the occurrences of the species id 7237 and 9593 are replaced with
    46245 and 9595 ids, respectively.
###  A *template.properties* file example
```
[DEFAULT]
oma_sparql_endpoint = https://sparql.omabrowser.org/sparql/

[OUTPUT]
file_directory_output = ./Bgee_data/v15/gene_homology_data/

[RDF_VOCABULARY]
ncbi_gene_crossref_property = lscr:xrefNCBIGene
ensembl_gene_crossref_property = lscr:xrefEnsemblGene

[SPECIES]
species_file = ./species_lists/all_speciesId_bgee15.txt
ncbi_species_file = ./species_lists/ncbi_speciesId_bgee15.txt
flybase_species_file = ./species_lists/flybase_speciesId_bgee15.txt
taxon_namespace = http://purl.uniprot.org/taxonomy/
```

### Generating orthology and paralogy files, incrementally
This app writes in files (*orthologs_pairs_with_results.tmp* and *paralogs_pairs_with_results.tmp*) the species already
processed in order to be used as an input for avoiding 
to recompute them in a later execution. For example, in case of adding a new species,
the user may want to solely generate the orthology and paralogy files for the missing pairs of species.
The *orthologs_pairs_with_results.tmp* and *paralogs_pairs_with_results.tmp* files are in the 
TEMPORARY_CONFIG_FILES_DIRECTORY directory as assigned in the [makefile.properties](makefile.properties).

To enable this functionality set the **REMOVE_SPECIES_PAIRS** to True in the [makefile.properties](makefile.properties)
file as show below. See further information about this feature in the **Optional parameters** list in Section
"The configuration file".
```
REMOVE_SPECIES_PAIRS = True
```

**Remark:** If it is not desirable anymore to generate orthology and paralogy incrementally, 
those *.tmp* files should be deleted. To start the homology generation from scratch, first, remove the .tmp files 
or provide a different folder for the temporary files by setting a different directory path to 
the TEMPORARY_CONFIG_FILES_DIRECTORY parameter.

### SPARQL queries 
The [QueryCatalog](QueryCatalog.py) python module contains the SPARQL queries used to extract the OMA orthology and 
paralogy pairwise relations from the hierarchical orthologous groups (HOGs). 

The RDF serialisation of the OMA Hierarchical Orthologous Groups is based on the 
[ORTH ontology specification](http://qfo.github.io/OrthologyOntology/) in the context of the Quest for Orthologs (QfO)
consortium (see GitHub repository [here](https://github.com/qfo/OrthologyOntology)).

To cross-reference other resources, this SPARQL endpoint contains annotation property assertions defined
by a first draft of the life-sciences cross-reference (LSCR) ontology that is available to download
at the Quest for Orthologs github repository [here](https://github.com/qfo/OrthologyOntology) - 
[lscr.ttl](https://github.com/qfo/OrthologyOntology/blob/master/lscr.ttl) file.

