Requirements: having successfully run the steps in [database creation](../db_creation/).

Requirements2: Check that the species are available in Ensembl or EnsemblMetazoa.

Goal: Insert the species used in Bgee, the related taxonomy from NCBI, and sex information about species.

## Details

* The file `source_files/species/bgeeSpecies.tsv` contains the species used in Bgee. This is the file to modify to add/remove a species.
* If you add/remove some species, you need to also update the files `pipeline/db_creation/insert_data_sources_to_species.sql` and `source_files/species/insert_species_sex_info.sql`.

* Details about `source_files/species/bgeeSpecies.tsv` file:
The first line is a header line, it must contain the following columns:
  * `speciesId`: The NCBI taxonomy species ID (e.g., 9606 for human).
  * `genus`: genus of the species.
  * `species`: species name of the species.
  * `speciesCommonName`: species common name
  * `genomeFilePath`: path to the genome file of the species on the Ensembl FTP. If no genome
  available for this species, this can point to the genome of a closely related species
  (e.g., use of chimpanzee genome for bonobo)
  * `genomeVersion`: genome version used
  * `dataSourceId`: ID of the data source providing the genome (currently, either Ensembl or EnsemblMetazoa)
  * `genomeSpeciesId`: the NCBI taxonomy ID of the species whose genome is being used. In most cases,
  it is the same as `speciesId`, but for some species with no genome available, it is possible
  to use the genome of a closely related species (e.g., use of chimpanzee genome for bonobo).
  * `fakeGeneIdPrefix`: If the genome of another species is used, the prefixes of the Ensembl gene IDs
  of this other species will be replaced with the provided prefix.
  * `keywords`: a list of keywords/alternative names associated to this species, separated by the character `|`.

Note that it is from this file that the information about names of species are obtained, so, no errors allowed in it.

If a line starts with `#`, it is commented and the species will not be inserted

* This pipeline step requires the NCBI taxonomy, provided as an ontology.
  * We cannot use the [http://www.obofoundry.org/cgi-bin/detail.cgi?id=ncbi_taxonomy official taxonomy ontology] because, as of Bgee 13, it does not include the last modifications that we requested to NCBI, and that were accepted (e.g., addition of a <i>Dipnotetrapodomorpha</i> term). Also, to correctly infer taxon constraints at later steps, we need this ontology to include disjoint classes axioms between sibling taxa, as explained in a Chris Mungall [http://douroucouli.wordpress.com/2012/04/24/taxon-constraints-in-owl blog post]. The default ontology does not include those.

  * This pipeline step is thus capable of generating its own version of the NCBI taxonomy ontology, in the exact same way as for the official ontology, as described [on the OBOFoundry wiki](http://www.obofoundry.org/wiki/index.php/NCBITaxon:Main_Page#Methods) (see notably the [Makefile](https://sourceforge.net/p/obo/svn/HEAD/tree/ncbitaxon/trunk/src/ontology/Makefile) generating the ontology). It is based on files available from the NCBI FTP (ftp://ftp.ebi.ac.uk/pub/databases/taxonomy/ taxonomy.dat). The code to generate disjoint classes axioms is based on the code from the [owltools](https://github.com/owlcollab/owltools) Java class `owltools.cli.TaxonCommandRunner`, in the module `OWLTools-Runner`.

  * This custom taxonomy will include the species used in Bgee and their ancestors, the taxa used in our annotations and their ancestors, the taxa used in Uberon an their ancestors. To extract taxa used in Uberon, we use the `ext` version (this is the one containing more taxa).

  * Note that the generation of the taxonomy requires about 15Go of memory.

* This pipeline steps will insert sex information about species, see `source_files/species/insert_species_sex_info.sql`.

* This pipeline steps will insert sex information about species, see `source_files/species/insert_species_sex_info.sql`.

## Data generation

* If it is the first time you execute this step in this pipeline run:
  `make clean`

* Run Makefile:
  `make`

* Modify the file `pipeline/db_creation/update_data_sources.sql`: you need to add the last modification date of the taxonomy used. This information can be found by looking at the file `taxonomy.dat` at ftp://ftp.ebi.ac.uk/pub/databases/taxonomy/.

## Data verification

* `generated_files/species/step_verification_RELEASE.txt` should contain: the total number of species inserted, the total number of taxa inserted, the number of taxa inserted that represent the least common ancestor of at least two species used in Bgee; the complete list of species ordered by their ID; the complete list of taxa least common ancestor, ordered by their position in the taxonomy (root to leaf); the complete list of taxa ordered by their position in the taxonomy (root to leaf).
  * Compare the species list between releases.
  * Check that there is no species with missing sex information. Otherwise, you need to update the file `source_files/species/insert_species_sex_info.sql`.
  * The taxa should be displayed ordered from root to leaf, and taxa of a same level should be ordered by alphabetical order of their scientific name. Verify it is correct, it is important.

* The following files should have been generated:
  * `generated_files/species/bgee_ncbitaxon.owl`, our custom taxonomy ontology
  * `generated_files/species/annotTaxIds.tsv`, a TSV file containing the IDs of the taxa used in our annotations
  * `generated_files/species/allTaxIds.tsv`, a TSV file containing the IDs of the species used in Bgee, of the taxa used in our annotations, of the taxa used in Uberon.

## Error handling

* You can have an exception thrown, saying that a specified taxon does not exist in the taxonomy ontology, for instance `java.lang.IllegalArgumentException: Taxon NCBITaxon:71164 was not found in the ontology`. It likely means that an incorrect/deprecated taxon is used in Uberon. Remove the ID of the taxon in the file `allTaxIds.tsv` (so if the exception is related to a taxon `NCBITaxon_71164`, remove the ID `71164`). Check if the taxon ID is present in the file `annotTaxIds.tsv`, if it is, the error is on us. Otherwise, report the problem on the Uberon tracker, if you identified the taxon in Uberon. You will most likely need to manually modify the ontology to remove the offending taxa in the mean time.

* If you need to re-insert the data in the database, use `make deleteSpeciesAndTaxa`, then `make clean`.

## Other notable Makefile targets
* generate the TSV files containing the taxa used in out annotations:
  `make ../../generated_files/species/annotTaxIds.tsv`
* generate the TSV files containing all taxa used in Bgee, in our annotations, in Uberon
  `make ../../generated_files/species/allTaxIds.tsv`
* generate the taxonomy ontology:
  `make ../../generated_files/species/bgee_ncbitaxon.owl`
* Remove the species and taxa from the database (this is not done when calling `clean`, to avoid wiping the database by accident)
  `make deleteSpeciesAndTaxa`
