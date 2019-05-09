# Insert genes & xrefs

* **Requirements**: having successfully run the step insert species and taxonomy AND database creation.
* **Goal**:         insert the genes used in Bgee.

## Details
### General informations
* This script will insert species gene information based on already inserted species (in the table `species`). It takes between 15 and 60 min per species.
* It will fill the gene table with `geneId`, `geneName` and `geneDescription` (+ `geneBioTypeId` + `speciesId`).
* It will also fill the tables `geneBioType`, `geneOntologyTerm` + `geneOntologyTermAltId` (and provides a file for obsolete GO terms), `geneNameSynonym`, `geneToGeneOntologyTerm`, and `geneXRef`. The insertion in `geneXRef` makes use of the `dataSource` table.

### Ensembl genes

* It uses the Ensembl API.
* **Important note regarding bonobo genes**: for bonobo, we take the same genome as chimpanzee (there is no bonobo genome in Ensembl, and it is debatable whether bonobo and chimp represent the same species). We use a SQL query at the end of the Makefile, to duplicate all chimpanzee genes, while providing new IDs. Consequently, it also duplicates the entries in `geneNameSynonym`, `geneToTerm` and `geneToGeneOntologyTerm` (with the appropriate IDs). **As a result, if you add or modify any fields in any of the tables `gene`, `geneNameSynonym`, `geneToTerm`, or `geneToGeneOntologyTerm`, you might need to modify the query used in this Makefile.**

### Non Ensembl genes

* Genome has been annotated with Maker and Blast2GO and the GFF file need to have the Maker and Blast2GO syntax. Example:

```
2_Tcm_b3v06_scaf001768	maker	mRNA	92518	119828	.	+	.	ID=TCM_04500-RA;Parent=TCM_04500;Name=TCM_04500-RA;Alias=maker-2_Tcm_b3v06_scaf001768-augustus-gene-0.2-mRNA-1;_AED=0.27;_QI=0|0|0|0.71|0.83|0.71|7|0|365;_eAED=0.27;ontology_term=GO:0016020,GO:0004767,GO:0005525,GO:0007165,GO:0003924;topblasthit_gene=;topblasthit_description=ras-related protein Ral-a;topblasthit_taxon=Halyomorpha halys;topblasthit_xref=XP_014294062.1;
```

The `description` field (last column) need to have the entries `ontology_term`, `topblasthit_gene`, `topblasthit_description`,  `topblasthit_taxon` and  `topblasthit_xref`, otherwise the tables will be filled with `undef`.

* The annotation file need to have the same base name than the genome file (with `_vbgee`) and has to be in the same directory (ie. `genomeFilePath` of the `species` table). Example:

Genome : `timema_data/1_Tdi_b3v08.fasta`
Annotation: `timema_data/1_Tdi_b3v08_vbgee.gff`

* Genome has only protein coding genes annotation. Therefore, Ensembl genes has to be inserted first in order to retrieve the `geneBioTypeId` of the `gene` table.

* Documentations:
	* [Maker](https://www.yandell-lab.org/software/maker.html)
	* [Blast2GO](https://www.blast2go.com/)
	
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
```sql
SELECT x.geneId FROM gene g, geneXRef x WHERE g.geneId = x.geneId AND g.geneBioTypeId = (SELECT geneBioTypeId FROM geneBioType WHERE geneBioTypeName='miRNA') AND x.dataSourceId = (SELECT dataSourceId FROM dataSource WHERE dataSourceName='ZFIN');
```
* Before Bgee 13, Ensembl did not provide XenBase Xref mapping. Check if this is done with:
```sql
SELECT x.geneId FROM gene g, geneXRef x WHERE g.geneId = x.geneId AND x.dataSourceId = (SELECT dataSourceId FROM dataSource WHERE dataSourceName='XenBase');
```

## Error handling
* Ensembl does NOT provide gene information the same way for each species, especially for non-Vertebrate species (_Drosophila melanogaster_ and _C. elegans_) and/or model organisms that have a non-Ensembl reference database (e.g. _Zebrafish_ or _Xenope_). Consequently Xrefs used in Bgee (linked in `dataSource` table) are not available the same way in Ensembl. You may need to add some extra aliases in the [insert_genes.pl](insert_genes.pl) script in order to catch all Xrefs you need. **You may have to do that for each new species as well as for each new Ensembl release.** E.g. not all Xrefs are available in the _zfin_ source, some others are available in `zfin_id`:
```perl
$InsertedDataSources{'zfin_id'}          = $InsertedDataSources{'zfin'};
$InsertedDataSources{'flybasename_gene'} = $InsertedDataSources{'flybase'};
$InsertedDataSources{'flybasecgid_gene'} = $InsertedDataSources{'flybase'};
$InsertedDataSources{'flybase_symbol'}   = $InsertedDataSources{'flybase'};
$InsertedDataSources{'wormbase_gene'}    = $InsertedDataSources{'wormbase'};
$InsertedDataSources{'xenopus_jamboree'} = $InsertedDataSources{'xenbase'};
```

## Other notable Makefile targets

* Update `gene` table with count of genes in database with an identical Ensembl gene ID (in Bgee, for some species with no genome available, we use the genome of a closely-related species, such as chimpanzee genome for analyzing bonobo data. For this reason, a same Ensembl gene ID can be mapped to several species in Bgee.):
    `make sameIdGeneCount`

# Gene Homology
* **Requirements**: **having contacted Adrian Altenhoff <adrian.altenhoff@inf.ethz.ch> one month in advance to request an update of the OMA HOGs based on our list of species. This is different from the data available from their download page...** Having successfully run the step insert genes.
* **Goal**: insert the OMA Hierarchical Orthologous Groups, and mirBase families.

## Details
* NEED TO THINK ABOUT HOW TO HANDLE MIRBASE FAMILIES, OUTSIDE OF THE TABLES DESIGNED FOR OMA HOG. OR FIND A WAY TO INTEGRATE THEM INTO THE OMA HOG TABLES?

## non Ensembl modifications

* Creation of the directory `generated_files/genes/`.
* Some parts are in common between `insert_genes.pl` and `insert_genes_nonEnsembl.pl` and have been factored in `Utils_insert_genes.pm`.
* For non Ensembl genes, `geneToTerm` table is not filled.
