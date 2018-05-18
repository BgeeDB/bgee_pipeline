# Insert ESTs
* **Requirements**: having successfully run the step insert species and taxonomy AND database creation AND insert genes.
* **Goal**: insert EST data used in Bgee.


## Details
* This step will insert normal EST libraries and their stage/organ.
* It will insert UniGene information from Ensembl, or from a mapping file, for EST-annotated species in Bgee. Currently only _human_, _mouse_, _zebrafish_ and _xenopus_ are concerned. _Drosophila melanogaster_ is also concerned through a mapping file done by blast in a special target of the Makefile.
  * Download of UniGene information requires about 3GB of disk space. The [Makefile.common](../Makefile.common) `$(TMPDIR)` variable points to a directory with enough free disk space.
  * Insertion of ESTs for _Drosophila melanogaster_ requires an extra mapping step because the mapping between FlyBase genes and UniGene clusters is not available in Ensembl/BioMart. So need to create a mapping file by BLAST (`dmel_mapping` step in the Makefile)
* It will insert miRNA EST information from [smirnadb](http://www.mirz.unibas.ch/smiRNAdb/). Only _human_, _mouse_, _zebrafish_ and _Drosophila melanogaster_ are currently concerned.
  * Please note that data from smirnadb are not updated anymore. So they are frozen in the svn and not re-downloaded & re-treated by the `Makefile`.
  * The file [S.xls](http://www.mirz.unibas.ch/cloningprofiles/resources/S.xls), once downloaded, has to be saved as tsv (change the name to `S.csv`: use openOffice to generate the tsv file, in UTF-8, with Tab delimiter, without quotes as field delimiter) and placed in `PIPELINEROOT/EST`. **Remove from this file the whole section concerning Rat libraries** (some have the same name than human libraries and the file is a mess to parse).
* It will fill `expression` table and update the field `estData` (quality) in `expressedSequenceTag` table.


## Data generation
* If it is the first time you execute this step in this pipeline run:
```
make clean
```

* Run Makefile:
```
make
```

**WARNING**: actually, before running `expression_est`, you should check the generated file `check_conditions`, to detect invalid conditions not supposed to exist in the related species. See 'Details' section of [pipeline/post_processing/README.md](../post_processing/README.md) for an explanation on how to fix such issues (in case the annotations were not incorrect).


## Error handling
* [insert_miRNA_est.pl](insert_miRNA_est.pl) can return warnings (mainly for human and mouse) due to MySQL uppercase/lowercase un-sensitiveness. E.g.:
  * DBD::mysql::st execute failed: Duplicate entry 'hsa-mi8-15' for key 'PRIMARY' at insert_miRNA_est.pl line 119.


## Other notable Makefile targets
* `make smiRNAdb` will download files from smirnadb. They are **NOT** updated anymore. So this target is **NOT** done by default!
  * The target will download the latest version of `S.xls`. Save it as a tsv (keep the name as `S.csv`, use openOffice to generate the tsv file, in UTF-8, with Tab delimiter, without quotes as field delimiter). The target will place it in [pipeline/ESTs](.). Remove from this file the whole section concerning Rat libraries (some have the same name than human libraries and the file is a mess to parse).
  * The target will download the latest version of the files `Report_x.csv` and place them in [pipeline/ESTs](.).
* `make dmel_mapping` will download _D. melanogaster_ cdna in fasta from Ensembl/BioMart and create a BLAST db from it. Then it will download _D. melanogaster_ UniGene clusters, and blast them against the previously created BLAST db to build the mapping file.
