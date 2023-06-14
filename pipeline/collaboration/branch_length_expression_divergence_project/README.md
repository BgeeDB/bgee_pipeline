Requirements: having successfully run all previous steps, this is a post-processing steps.

Goal: generate expression data for the project in collaboration with Dessimoz lab,
about the relation between branch length and expression divergence.

## Description

Participants to the project: Bastian Frederic, Dessimoz Christophe, Glover Natasha,
Mendes de Farias Tarcisio, Sima Ana-Claudia, Warwick Vesztrocy Alex.

Other contributors: Rech de Laval Valentine.

The aim of this project is to look at the relation between sequence and expression divergence.
Bgee provides expression data in tissues homologous at specific taxonomic levels.

The relations of historical homology between anatomical entities can be retrieved from the
[github anatomical similarity annotations](https://github.com/BgeeDB/anatomical-similarity-annotations)
(see `develop` branch for the latest version, notably with updated documentation).

## Description of the output files

Files can be found on our FTP: https://www.bgee.org/ftp/current/collaboration/branch_length_expression_divergence/.

### Present/absent expression calls

#### Columns

The first column contains the Ensembl IDs of genes for the species part of the requested taxon.
The following columns each contain present/absent expression calls in a specific group
of homologous anatomical entities.

Note that each group of homologous anatomical entities can contain 1 or more entities.
For instance, `retinal ganglion cell` and `rhabdomeric photoreceptor cell` both derive from
a same common ancestral structure in the ancestor of Bilateria; the data in these anatomical entities
are thus merged in a same group to compare expression data between bilaterian species.
Another example is `nervous system`, a structure homologous in all Bilateria; this homology group
thus contains only this tissue.

When a group of homologous anatomical entities contain more than 1 entity, their names
in the column header are separated with ' - '.

#### Rows

Each row represents the present/absent expression calls for each gene, in each group
of homologous anatomical entities. The value taken in each cell is one of:

* `EXPRESSED`: the gene has expression detected in this group of homologous anatomical entities
* `NOT_EXPRESSED`: there are data for this gene in this group of homologous anatomical entities,
indicating that the gene is not expressed in this group
* `NO DATA`: no data for this gene in this group of homologous anatomical entities

#### Statistics on expression calls

Each file is accompanied with a file with same name but starting with "stats_", providing various
statistics about the content of the present/absent expression call file.

### Expression calls with rank

To do/to discuss

### Expression calls with quantitative expression level categories relative to the min/max expression
level of genes

[See details](../oncomx#computation-of-the-expression-levels).
To do/to discuss

### Expression calls with quantitative expression level categories relative to the min/max expression
level in anatomical entities

[See details](../oncomx#computation-of-the-expression-levels).
To do/to discuss

## Pipeline

* `make branch_length_expression_divergence`: generates the files for this project

## Data generation

* If it is the first time you execute this step in this pipeline run:
  `make clean`

* Run Makefile:
  `make`

## Data verification

* Check the files in `generated_files/collaboration/branch_length_expression_divergence/`, notably the verification file `branch_length_expression_divergence`.
