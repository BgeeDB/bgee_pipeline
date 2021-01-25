# TODO for Bgee 15

## Uberon

Check why a relation such as `SubClassOf(ObjectIntersectionOf(<http://purl.obolibrary.org/obo/UBERON_0002316> ObjectSomeValuesFrom(<http://purl.obolibrary.org/obo/BFO_0000050> <http://purl.obolibrary.org/obo/NCBITaxon_9443>)) <http://purl.obolibrary.org/obo/UBERON_0002437>)` has not been inserted in the database. Does it work when we redo the insertion? Check indirect relations not reached by a chain of direct relations.

## Annotations

Use the strain mapping file in the pipeline for inserting conditions

## Bgee lite

* Add expression rank and expression score info
* Add propagation information:
  * See email 19.11.19 02:04, Tom Conlin
  * We could add the values from `PropagationState` (see `org.bgee.model.expressiondata.Call.ExpressionCall#getDataPropagation()`)
  
## ExperimentExpression, experiment counts

* Get rid of all the "xxxExperimentExpression" tables and related code inserting data in them. For Bgee 15, confidence levels will be based on corrected p-values, not on number of experiments.
* Modify the globalExpression table and related code accordingly.

## Affymetrix

* Rerun Affymetrix analyses to be able to store p-values (Sara, for new FDR correction)

## RNA-Seq

* Do not produce absent calls for some gene biotypes, depending on the library type
* Same for the ranks: for now, we consider that all genes that have received
  at least one read in any library are all always accessible to rank computation in all libraries.
  
* Have different calls quality depending on the threshold intergenic/genes

* Check discarded libraries, see which one should be recovered

* Globin reduction on blood samples: we need a test to determine whether blood samples
had globin reduction or not. Let's implement the test and look at the distribution
of samples with/without reduction. Notes about that in the Bgee meeting minutes
from 2020-04-07
  * for all samples that are blood, we will run a test to check the globin depletion status
  * insert the information in the database. Maybe a specific column,
  or same information as the type of targeting of the library (miRNA, lncRNA, etc)
  * either the depletion will be known from annotation, provided by the data providers, or from the test.
  * add the result of the test in the rnaSeqInfo file already used by the pipeline.

## scRNA-Seq

### target based pipeline
  * check all files created and put their names in the Makefile.common (with variable parts like SPECIES, 
  LIBRARY_ID that will be update on the fly)
  * parallelize kallisto_bus (oher rules could be parallelize too but they take less than 2 days to run as for bgee 15.0)
  * check again all created files to be sure it is not possible to duplicate data inside of them
  * add file checking library already processed. Will allow not to reprocess it when rule is run again. Very important to do it for kallisto_bus
  * clean directory_names and and homogenize variable/path/script names to keep one naming approach (camel or snake)
  * generate target based calls using BgeeCall

## Post-processing

* post-processing to remove genes never seen expressed anywhere.
Note: this filter already exists for Affy and RNA-Seq data independently. EST only produce present calls.
Such a situation should then only happens from in situ data where only absence of expression of a gene was reported,
and with no present calls from other data types. => Do we really need a post-processing filtering step for this?
  
## Issues

Check https://github.com/BgeeDB/bgee_pipeline/issues
