**Requirements**: having successfully run all other expression insertion steps for all data types

**Goals**:
* make a global verification of validity of conditions used
* generate global propagated expression calls
* compute and insert gene expression ranks in database

## Details

The first step is to run a verification on expression conditions used, using the step `generated_files/post_processing/conditions_not_existing`. It detects conditions using anatomical entities and/or developmental stages not supposed to exist in the related species. If any invalid condition is detected, and the errors come from the taxon constraints, not from our annotations, see SQL file [fix_anat_taxon_constraints.sql](fix_anat_taxon_constraints.sql) in this folder for a list of SQL commands to run to fix the issues. In case you want to remap annotations, see SQL file [remap_conditions.sql](remap_conditions.sql) in this folder.

Then, global expression calls are generated

Second, ranks are computed for all data types and condition parameter combinations.

Third, ranks are normalized over all data types and species.

### Information about rank computations

The main steps are the following:

1. we compute ranks using different approaches for each data type.
2. we "normalize" ranks inside each data type (except for RNA-Seq, no normalization)
3. we "normalize" ranks over all genes, conditions and data types for each species
4. we compute a global weighted mean rank for each gene in each condition over all data types

#### Affymetrix data

[See Affymetrix rank script](ranks_affymetrix.pl)

1. For each Affymetrix chip, compute franctional ranks for each gene, based on the highest signal intensity
from all the probesets mapped to a given gene.
2. Normalize ranks between chips in a same condition and same species:
    1. First, for each chip type, compute its max rank over all conditions.
    2. Then, normalize ranks of each chip, based on the max rank of the corresponding chip type,
    as compared to the max of the max ranks of other chip types present in the same condition.
    The idea is to correct for the different genomic coverage of different chip types.
    We do not normalize simply based on the max rank in a given condition, but based on the max rank of the chip types
    represented in that condition, to not penalize conditions with a lower number of expressed genes,
    and thus higher number of ex-aequo genes, and lower fractional max ranks.
3. Compute weighted mean of normalized ranks per gene and condition, weighted by the number of distinct ranks
in each chip: we assume that chips with higher number of distinct ranks have a higher power for ranking genes.
The higher power at ranking genes of a chip is taken into account by weighting the mean 
by the number of distinct ranks in the chip, not by "normalizing" away chips with lower max ranks, 
again, not to penalize conditions with lower numbers of expressed genes.
4. Store the max of max ranks of chip types represented in each condition
(will allow to normalize ranks between conditions and data types), and sums of numbers of distinct ranks 
per gene and condition (used afterwards to compute the global weigthed mean rank over all data types in a condition,
in the application).

#### RNA-Seq data

[See RNA-Seq rank script](ranks_rnaseq.pl)

1. Identify the valid set of genes that should be considered for ranking in all libraries: 
the set of all genes that received at least one read over all libraries in Bgee. 
2. Compute gene fractional ranks for each RNA-Seq library, based on TPM values. 
3. Compute weighted mean of ranks per gene and per condition, weighted by the number of distinct ranks
in each library: we assume that libraries with a higher number of distinct ranks have a higher power for ranking genes. 
Note that we do not "normalize" ranks between samples before computing the mean, as for Affymetrix data: 
all libraries are used to produce ranking over always the same set of genes in a given species, 
so the genomic coverage is always the same, and no "normalization" is required. The higher power 
at ranking genes of a library (for instance, thanks to a higher number of mapped reads) 
is taken into account by weighting the mean by the number of distinct ranks in the library, 
not by "normalizing" away libraries with lower max ranks; this would penalize conditions 
with a lower number of expressed genes, and thus with more ex-aequo ranked genes, corresponding 
to genes receiving 0 read. 
4. Store the max ranks in each condition (will allow to normalize ranks between conditions and data types),
and the sum of distinct rank counts per gene and condition (used afterwards to compute the global weigthed mean rank
over all data types in a condition, in the application)

#### EST data

[See EST rank script](ranks_est.pl)

1. Compute "dense ranking" of the genes in each condition, based on the EST counts.
All libraries in a condition are considered together (no ranking first performed per library,
as for other data types). This is because each library usually provides information about a very low number of genes.
A fractional ranking was not appropriate, because there are many ex-aequo, leading to have an artificially
high max rank. This was giving too much weight to EST data when computing a global weighted mean rank
over all data types.
2. For each condition, retrieve the max rank, and store it (will allow to normalize ranks between conditions and data types,
and to compute the global weigthed mean rank over all data types in a condition, in the application)

#### In situ hybridization data

[See in situ hybridization rank script](ranks_in_situ.pl)

1. First, compute a score for each gene-condition, based on the detection flag
of in situ evidence (present, absent), and the quality level (high quality, poor quality).
This detection flag and quality level are retrieved from the source MOD database.
Each spot is given the following score:
    * present - high quality = 1
    * present - low quality = 0.5
    * absent - low quality = -0.5
    * absent - high quality = -1
Then we simply sum up the scores from the spots for each gene in each condition
2. Compute "dense ranking" of the genes in each condition, based on the score computed.
All experiments are all considered together (no ranking first performed per experiment,
as for other data types). This is because each in situ experiment usually studies
a very limited number of genes. A fractional ranking was not appropriate, because there are many ex-aequo,
leading to have an artificially high max rank. This was giving too much weight to in situ hybridization data
when computing a global weighted mean rank over all data types.
3. For each condition, retrieve the max rank, and store it (will allow to normalize ranks between conditions and data types,
and to compute the global weigthed mean rank over all data types in a condition, in the application)

#### "Normalization" over all data types and conditions

[Normalization script](normalize_ranks.pl)

1. Retrieve the max mean rank over all data types, conditions, and genes, for each species independently *Frederic check this*
2. Normalize mean ranks for each gene in each condition, based on the max rank for this data type in this condition,
as compared to the global max rank over all conditions and data types.

#### Computation of global weigthed mean rank for each gene and condition over all data types

This computation is done in the Java Bgee API directly, to be able to choose the data types to consider
to compute the ranks, based on the user request.

TODO: provide the formula here


## Data generation

* If it is the first time you execute this step in this pipeline run:
  `make clean`

* Run Makefile:
  `make`

## Data verification

## Error handling

## Other notable Makefile targets

* Check validity of conditions. See 'Details' section for more information:
    `make ../../generated_files/post_processing/conditions_not_existing`

* Generate global expression calls for "anat entity" combination:
    `make ../../generated_files/post_processing/globalAnatEntityExpression`

* Generate global expression calls for "anat entity - stage" combination:
    `make ../../generated_files/post_processing/globalAnatEntityStageExpression`
