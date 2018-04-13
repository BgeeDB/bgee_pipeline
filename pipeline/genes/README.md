# Insert genes & xrefs

* **Requirements**: having successfully run the step insert species and taxonomy AND database creation.
* **Goal**:         Insert the genes used in Bgee.

## Details
* This script will insert species gene information based on already inserted species (in the table species). It takes between 15 and 60 min per species.
* It uses the Ensembl API.
* It will fill the gene table with _geneId_, _geneName_ and _geneDescription_ (+ _geneBioTypeId_ + _speciesId_).
* It will also fill the tables _geneBioType_, _geneOntologyTerm_ + _geneOntologyTermAltId_ (and provides a file for obsolete GO terms), _geneNameSynonym_, _geneToGeneOntologyTerm_, and _geneXRef_. The insertion in _geneXRef_ makes use of the _dataSource_ table.
* **Important note regarding bonobo genes**: for bonobo, we take the same genome as chimpanzee (there is no bonobo genome in Ensembl, and it is debatable whether bonobo and chimp represent the same species). We use a SQL query at the end of the Makefile, to duplicate all chimpanzee genes, while providing new IDs. Consequently, it also duplicates the entries in _geneNameSynonym_, _geneToTerm_ and _geneToGeneOntologyTerm_ (with the appropriate IDs). **As a result, if you add or modify any fields in any of the tables gene, geneNameSynonym, geneToTerm, or geneToGeneOntologyTerm, you might need to modify the query used in this Makefile.**

## Data generation
* If it is the first time you execute this step in this pipeline run:
```
make clean
```
* Run Makefile:
```
make
```

## Data verification
* Before Bgee 13, mirbase Xref were not provided by Ensembl for Zebrafish. Check it at this end of the pipeline with
``` sql
SELECT x.geneId FROM gene g, geneXRef x WHERE g.geneId = x.geneId AND g.geneBioTypeId = (SELECT geneBioTypeId FROM geneBioType WHERE geneBioTypeName='miRNA') AND x.dataSourceId = (SELECT dataSourceId FROM dataSource WHERE dataSourceName='ZFIN');
```
* Before Bgee 13, Ensembl did not provide XenBase Xref mapping. Check if this is done with:
``` sql
SELECT x.geneId FROM gene g, geneXRef x WHERE g.geneId = x.geneId AND x.dataSourceId = (SELECT dataSourceId FROM dataSource WHERE dataSourceName='XenBase');
```

## Error handling
* Ensembl does NOT provide gene information the same way for each species, especially for non-Vertebrate species (_Drosophila melanogaster_ and _C. elegans_) and/or model organisms that have a non-Ensembl reference database (e.g. _Zebrafish_ or _Xenope_). Consequently Xrefs used in Bgee (linked in _dataSource_ table) are not available the same way in Ensembl. You may need to add some extra aliases in the **insert_genes.pl** script in order to catch all Xrefs you need. **You may have to do that for each new species as well as for each new Ensembl release.** E.g. not all Xrefs are available in the _zfin_ source, some others are available in _zfin_id_:
``` perl
$InsertedDataSources{'zfin_id'}          = $InsertedDataSources{'zfin'};
$InsertedDataSources{'flybasename_gene'} = $InsertedDataSources{'flybase'};
$InsertedDataSources{'flybasecgid_gene'} = $InsertedDataSources{'flybase'};
$InsertedDataSources{'flybase_symbol'}   = $InsertedDataSources{'flybase'};
$InsertedDataSources{'wormbase_gene'}    = $InsertedDataSources{'wormbase'};
$InsertedDataSources{'xenopus_jamboree'} = $InsertedDataSources{'xenbase'};
```

## Other notable Makefile targets

* Update gene table with count of genes in database with an identical Ensembl gene ID (in Bgee, for some species with no genome available, we use the genome of a closely-related species, such as chimpanzee genome for analyzing bonobo data. For this reason, a same Ensembl gene ID can be mapped to several species in Bgee.):
    `make sameIdGeneCount`

# Gene Homology
* **Requirements**: **Having contacted Adrian Altenhoff <adrian.altenhoff@inf.ethz.ch> one month in advance to request an update of the OMA HOGs based on our list of species. This is different from the data available from their download page...** Having successfully run the step insert genes.
* **Goal**: Insert the OMA Hierarchical Orthologous Groups, and mirBase families.

## Details
* NEED TO THINK ABOUT HOW TO HANDLE MIRBASE FAMILIES, OUTSIDE OF THE TABLES DESIGNED FOR OMA HOG. OR FIND A WAY TO INTEGRATE THEM INTO THE OMA HOG TABLES?

