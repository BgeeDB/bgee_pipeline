Requirements: having successfully run all other expression insertion steps for all data types

Goals:

* make a global verification of validity of conditions used
* generate global propagated expression calls
* compute and insert gene expression ranks in database

## Details

The first step is to run a verification on expression conditions used, using the step `../../generated_files/post_processing/conditions_not_existing`. It detects conditions using anatomical entities and/or developmental stages not supposed to exist in the related species. If any invalid condition is detected, and the errors come from the taxon constraints, not from our annotations, see SQL file `fix_anat_taxon_constraints.sql` in this folder for a list of SQL commands to run to fix the issues. In case you want to remap annotations, see SQL file `remap_conditions.sql` in this folder.

Then, global expression calls are generated

Second, ranks are computed for all data types and condition parameter combinations.

Third, ranks are normalized over all data types and species.

## Data generation

* If it is the first time you execute this step in this pipeline run:
  `make clean`

* Run Makefile:
  `make`

## Data verification

## Error handling

## Other notable Makefile targets

* Check validity of conditions. See `Details` section for more information:
    `make ../../generated_files/post_processing/conditions_not_existing`

* Generate global expression calls for "anat entity" combination:
    `make ../../generated_files/post_processing/globalAnatEntityExpression`

* Generate global expression calls for "anat entity - stage" combination:
    `make ../../generated_files/post_processing/globalAnatEntityStageExpression`
