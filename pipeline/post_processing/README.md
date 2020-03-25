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
2. we "normalize" ranks over all genes, conditions and data types for each species
3. we compute a global weighted mean rank for each gene in each condition over all data types

#### Affymetrix data

[See Affymetrix rank script](ranks_affymetrix.pl)

1. For each Affymetrix chip, compute franctional ranks for each gene, based on the highest signal intensity
from all the probesets mapped to a given gene.
2. Normalize ranks between chips in a same condition and same species:
    1. First, for each chip type, compute its max rank over all conditions.
    2. Then, normalize ranks of each chip, based on the max rank of the corresponding chip type, as compared to the max of the max ranks of other chip types present in the same condition.
    The idea is to correct for the different genomic coverages of different chip types.
    We do not normalize simply based on the max rank in a given condition, but based on the max rank of the chip types represented in that condition, to not penalize conditions with a lower number of expressed genes, and thus higher number of ex-aequo genes, and lower fractional max ranks.
    For gene `g`, condition `c`, chip type `t`, sample (chip) `s`
<code>
normalized_rank<sub>gs</sub> = (rank<sub>gs</sub> * (1 + max_rank<sub>c</sub> / max_rank<sub>t</sub>)) / 2
</code>

<samp>max_rank<sub>t</sub></samp>: max rank over all chips of type `t` over all conditions; <samp>max_rank<sub>c</sub></samp>: max of <samp>max_rank<sub>t</sub></samp> from the chip types `t` used in the condition `c`.
3. Compute weighted mean of normalized ranks per gene and condition, weighted by the number of distinct ranks in each chip: we assume that chips with higher number of distinct ranks have a higher power for ranking genes.
The higher power at ranking genes of a chip is taken into account by weighting the mean
by the number of distinct ranks in the chip (gave better results than by using the max rank of the chip).
4. Store the max of max ranks of chip types represented in each condition <samp>max_rank<sub>c</sub></samp> (will allow to normalize ranks between conditions and data types), and sums of numbers of distinct ranks per gene and condition (used afterwards to compute the global weigthed mean rank over all data types in a condition, in the application).

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
is taken into account by weighting the mean by the number of distinct ranks in the library
(gave better results than by using the max rank in the library).
4. Store the max ranks in each condition (will allow to normalize ranks between conditions and data types),
and the sum of distinct rank counts per gene and condition (used afterwards to compute the global weigthed mean rank over all data types in a condition, in the application). As of Bgee 14.1, the same max rank is considered in all conditions, because we assume that the same set of genes is accessible to all libraries.

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
and to compute the global weighted mean rank over all data types in a condition, in the application)

#### "Normalization" over all data types and conditions

[Normalization script](normalize_ranks.pl)

1. Retrieve the max rank for each species over all data types, conditions, and genes. 
It is the max rank over all samples and conditions, not the max mean rank per condition.
2. Retrieve the max rank for each data type and condition.
3. Normalize ranks for each gene and data type in each condition. It is based on the max rank for this data type in this condition, as compared to the global max rank over all conditions and data types:

* For gene `g`, condition `c`, data type `d`, species `s`

<code>
normalized_rank<sub>gcd</sub> = (rank<sub>gcd</sub> * (1 + max_rank<sub>s</sub> / max_rank<sub>cd</sub>)) / 2
</code>

* a condition `c` is considered to belong to a specific species `s`. For example, "brain" in human and "brain" in mouse are considered to be two different conditions; a gene `g` is considered to belong to a specific species `s` as well.
* <samp>rank<sub>gcd</sub></samp>:
  * for RNA-Seq data, the weighted mean of the fractional ranks of the gene from each library in the condition, weighted by the number of distinct ranks in the library, as described above.
  * for Affymetrix data, the weighted mean of the normalized fractional ranks of the gene from each sample in the condition, as described above (fractional ranks from each chip are first normalized between chip types and conditions, then used to compute means weighted by the number of distinct ranks on the chips).
  * for EST and in situ hybridization data, the dense rank of the gene in the condition, computed by pulling all samples in the condition together, as described above.
* <samp>max_rank<sub>s</sub></samp>, <samp>max_rank<sub>cd</sub></samp>: max rank over all samples, not an averaged rank.
  * For RNA-Seq data, the same max rank is considered in all conditions (because we assume that the same set of genes is accessible to all libraries), so <samp>max_rank<sub>cd</sub> = max_rank<sub>sd</sub></samp>.
  * For Affymetrix data, the max rank in a condition (<samp>max_rank<sub>cd</sub></samp>) is the max rank over the chip types used in the condition, not the max rank from the actual samples in the condition. The idea is to correct for the different genomic coverages of different chip types. We do not normalize simply based on the max rank in a given condition, but based on the max rank of the chip types represented in that condition, to not penalize conditions with a lower number of expressed genes, and thus higher number of ex-aequo genes, and lower fractional max ranks.

4. We blacklist some conditions and give to all genes the max rank of the species + 1, for instance, in the conditions using an anatomical term "unknown".


#### Computation of global weigthed mean rank for each gene and condition over all data types

This computation is done in the Java Bgee API directly, to be able to choose the data types to consider
to compute the ranks, based on the user request.

For a gene `g`, condition `c`:

* RNA-Seq: 
  * <samp>rna_seq_normalized_mean_rank<sub>gc</sub></samp>: weighted mean of the fractional ranks of the gene `g` over all libraries in the condition `c`, weighted by the number of distinct ranks in each library (see step RNA-Seq data); this weighted mean has then been normalized between conditions and data types (see normalization step above).
  * <samp>rna_seq_distinct_rank_sum<sub>c</sub></samp>: the sum of the number of distinct ranks in the libraries present in condition `c`
* Affymetrix:
  * <samp>affymetrix_normalized_mean_normalized_rank<sub>gc</sub></samp>: the fractional rank of the gene `g` on each chip in the condition `c` is normalized between chip types and conditions, then used to compute a weighted mean weighted by the number of distinct ranks in each chip (see step Affymetrix data); this weighted normalized rank has then been further normalized between conditions and data types (see normalization step above).
  * <samp>affymetrix_distinct_rank_sum<sub>c</sub></samp>: the sum of the number of distinct ranks in the chips present in condition `c`
* In situ hybridization:
  * <samp>in_situ_normalized_rank<sub>gc</sub></samp>: dense rank of the gene `g` in condition `c` computed by pulling together all the evidence sources for condition `c` (see step in situ hybridization data); this dense rank has then been normalized between conditions and data types (see normalization step above).
  * <samp>in_situ_max_rank<sub>c</sub></samp>: the max dense rank over all genes ranked in condition `c` by in situ data
* EST:
  * <samp>est_normalized_rank<sub>gc</sub></samp>: dense rank of the gene `g` in condition `c` computed by pulling together all the EST libraries in condition `c` (see step EST data); this dense rank has then been normalized between conditions and data types (see normalization step above).
  * <samp>est_max_rank<sub>c</sub></samp>: the max dense rank over all genes ranked in condition `c` by EST data

Depending on the data types selected, if all data types are selected:

<code>
Bgee_expression_rank<sub>gc</sub> = (rna_seq_distinct_rank_sum<sub>c</sub> * rna_seq_normalized_mean_rank<sub>gc</sub> + affymetrix_distinct_rank_sum<sub>c</sub> * affymetrix_normalized_mean_normalized_rank<sub>gc</sub> + in_situ_max_rank<sub>c</sub> * in_situ_normalized_rank<sub>gc</sub> + est_max_rank<sub>c</sub> * <samp>est_normalized_rank<sub>gc</sub>) / (rna_seq_distinct_rank_sum<sub>c</sub> + affymetrix_distinct_rank_sum<sub>c</sub> + in_situ_max_rank<sub>c</sub> + est_max_rank<sub>c</sub>)
</code>

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
