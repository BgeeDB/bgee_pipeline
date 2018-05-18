# Insert _in situ_
* **Requirements**: having successfully run the step database creation AND Uberon
* **Goal**: update _in situ_ data in Bgee.


## Details
* Get _in situ_ data for different species (and different sources/db)
  * Get _in situ_ from ZFIN
  * Get _in situ_ from Xenbase
  * Get _in situ_ from BDGP & from FlyBase
  * Get _in situ_ from MGI (from MouseMine web services)
  * Get _in situ_ from WormBase
* Prepare data
  * All data are processed in the same TSV format, used in the next step
* Insert _in situ_
  * Then clean some cases of the previous insertion
* Insert expression - no expression (_can be run on a single taxon with `--taxon_id` option_)
* While running the pipeline, there will be lots of inconsistencies or synchronization issues between species-specific ontologies or MODs and Uberon. It will be needed to manually add missing mappings. See https://github.com/obophenotype/uberon/issues/664, and especially https://github.com/obophenotype/uberon/issues/1288 and https://gitlab.isb-sib.ch/Bgee/bgee_pipeline/issues/68.

## Data generation
* If it is the first time you execute this step in this pipeline run:
```
make clean
```

* Run Makefile:
```
make
```

**WARNING**: actually, before running `insert_expression`, you should check the generated file `check_conditions`, to detect invalid conditions not supposed to exist in the related species. See 'Details' section of [pipeline/post_processing/README.md](../post_processing/README.md) for an explanation on how to fix such issues (in case the annotations were not incorrect).


## Data verification
* Check in `warnings.tsv` files that the in between stages socket does not return messages that could notify problems:
  * There may have some empty start-end stages in data sources: `grep -v '\[,\]' warnings.tsv`
  * There may have some unmapped terms: `grep -v 'Could not find any OWLClass corresponding to' warnings.tsv`
  * There may have problems in our custom ontology and in annotations from data providers: `grep -v 'The start stage provided is actually a successor of the end stage provided' warnings.tsv`
  * **Other messages have to be fixed in our custom ontology !!!**
* **Xenopus**: Check on the FTP that it is still the up-to-date directory: [ftp://ftp.xenbase.org/pub/GenePageReports/](ftp://ftp.xenbase.org/pub/GenePageReports/)
  * Some XenBase warnings may be due to out-dated mapping between Ensembl genes used by XenBase and current Ensembl version. Also, Xenbase modified its developmental ontology, and improperly assigned the already used ID XAO:1000009 (previously used for "adult" stage) to the new term "frog". In bgee_12 it messed up 2 spots, that was ok.
* **MGI**: MGI stores only [normal conditions](http://www.informatics.jax.org/mgihome/homepages/tabContents/GXD_Curators.shtml) expression data (wild-type and mutants, but no treatment, etc)
* Think to check how many 'BDGP' genes remained mapped: with time their old Ensembl mapping is more and more limited (and no remapping planed from their side).
  * The Makefile downloads the BDGP database and install it.
  * The Makefile generates a mapping from BDGP terms to FBbt (in Uberon now) terms.
  * This produces an annotation file `BDGP_terms_to_FBbt_terms_bgee_vRELEASE.xls` that needs to be checked by annotators to complete the mapping.
  * **Regenerate** the mapping a last time, to check if annotators did not forget anything, and if they did not map a term with an obsolete FBbt. In that case, the BDGP term will have no mapping again. So, regenerate the mapping based on the new `BDGP_terms_to_FBbt_terms.tsv`, as already described. Then, compare the two files.
  ```
  xlscat BDGP/BDGP_terms_to_FBbt_terms_bgee_vRELEASE.xls                             >new
  xlscat new_bgee/curation/expression_data/in_situ/bdgp/BDGP_terms_to_FBbt_terms.xls >ori
  diff ori new
  ```
  * **Very important!** Check in the new mapping file that all terms that contains "no staining" are listed in the script [pipeline/In_situ/BDGP/prepare_data.pl](BDGP/prepare_data.pl) in the hash `%no_staining`. Currently, it is:
  ```
  my %no_staining = ('489' => 1,  '490' => 1,  '491' => 1,  '492' => 1,  '493' => 1,  '494' => 1,);
  ```
  * To be sure to check all terms, you can do:
  ```
  grep 'stain' BDGP_terms_to_FBbt_terms.tsv
  ```
  * Only terms also present in `source_files/In_situ/insitu_annot.csv` must be checked. If they are not present in this file, they won't be inserted. So for instance, a term like "list_9-10:no staining" does not need to be in the `%no_staining` hash.
  ```
  grep 'stain' BDGP/insitu_annot.csv | grep -i 'list'
  ```
  * Check there are no error in the BDGP stderr output (`warnings.tsv_BDGP`)
  ```
  cat warnings.tsv_BDGP | grep -v 'because no mapping found for' | grep -v 'ambiguous mapping' | grep -v 'ambiguous cgname' \
  | grep -v 'Removing main_id ' | grep -v 'Removing annot id '
  ```

  * BDGP mapping information: check in the `tsv_BDGP.log` file. E.g.

  ```
  Bgee v12                                                      Bgee v13
       8561  experiments where found                                 9207  experiments where found
      52302  evidences   where found                                56340  evidences   where found
      81570  spots       where found                                88903  spots       where found
        347  BDGP organs are mapped to FBbt organs                    348  BDGP organs are mapped to FBbt organs
  
  After filtering:                                              After filtering:
       3407  experiments are left                                    3582  experiments are left
      15690  evidences   are left                                   16664  evidences   are left
      21475  spots       are left                                   22234  spots       are left
  ```

  * Check individual warnings/cases:

  ```
  # Number of spots removed because of unmapped anatomical structures (6'547 for Bgee_v12 - 7'516 for Bgee_v13)
  cat warnings.tsv_BDGP | grep 'because no mapping found for' | wc -l
  
  # Number of ambiguous mapping to anatomical terms (0 for Bgee 12 - 0 for Bgee_v13)
  cat warnings.tsv_BDGP | grep 'ambiguous mapping' | wc -l
  
  # Number of evidences removed because all spots were removed because of absent or ambiguous mappings (1'901 for Bgee 12 - 1'326 for Bgee_v13)
  cat warnings.tsv_BDGP | grep 'because all spots from this evidence were removed' | wc -l
  
  # Number of experiments removed because all evidences were removed (616 for Bgee 12 - 663 for Bgee_v13)
  cat warnings.tsv_BDGP | grep 'because all evidences from this experiment were removed' | wc -l
  
  # EST Id used several times, weird (2 for Bgee 12 - 2 for Bgee_v13)
  cat warnings.tsv_BDGP | grep 'because of est_id used several times' | wc -l
  
  # A gene name in Bgee corresponds to a cg name in BDGP, but used for another FlBbase ID (0 for Bgee 12 - 0 for Bgee_v13)
  cat warnings.tsv_BDGP | grep 'because of inconsistent CG name for a gene retrieved by flybase ID' | wc -l
  
  # An ID in Bgee corresponds to a FlyBase ID in BDGP, but used for another CG name (3 in Bgee 12 - 4 for Bgee_v13)
  cat warnings.tsv_BDGP | grep 'because of inconsistent flybase ID for a gene retrieved by CG name' | wc -l
  
  # cgname not found in Bgee (108 in Bgee 12 - 114 for Bgee_v13)
  cat warnings.tsv_BDGP | grep 'because of flybase_id corresponding to cgname not found in bgee' | wc -l
  
  # No flybase_id nor cgname provided (29 in Bgee 12 - 29 for Bgee_v13)
  cat warnings.tsv_BDGP | grep 'because of no cgname and no flybase_id' | wc -l
  
  # Both the flybase_id and the cgname are provided in BDGP, but a gene in Bgee retrieved by the BDGP flybase_id does not match the Bgee cgname, or   vice-versa
    (3'655 in Bgee 12 - 3'984 for Bgee_v13)
  cat warnings.tsv_BDGP | grep 'because of ambiguous cgname and flybase_id' | wc -l
  
  # Gene not present in Bgee (741 in Bgee 12 - 829 for Bgee_v13)
  cat warnings.tsv_BDGP | grep 'because gene not present in Bgee' | wc -l
  ```

  * Check unmapped structures: (in Bgee 12, 31 structures - in Bgee 13, 33 structures (all weird terms) should not be used ("gap", "segmentally repeated", ...). Check in `BDGP_terms_to_FBbt_terms.tsv` for their name)
  ```
  cat warnings.tsv_BDGP | grep 'because no mapping found for' | egrep -o " [0-9]+$" | sort | uniq
  ```


* **Flybase** does not provide _in situ_ quality information (available only for RNA Seq). It also provides very few absent annotation, in a non-structured way. So we don't catch them for now and let them mixed with "normal/regular" present annotation.

* **WormBase**: try to use WormMine next time. It is still in beta so try to ask for missing things before next release.
  * Currently WormBase (file from their acedb) provides very few absent annotation, in a non-structured way. So we don't catch them for now.


## Other notable Makefile targets
* `map_source` are targets that get the original mapping
* `tsv_source` are targets that prepare data in a single TSV format
* `insert_insitu` inserts all tsv_source files found in Bgee
* `insert_expression` inserts expression / no expression in Bgee
* `all` will recapitulate all targets and checks/creates statistics about _in situ_ insertions


## Data Quality Information
* **MGI**

| code | quality        | percentage in source file | quality in Bgee |
|------|----------------|---------------------------|-----------------|
| -2   | Not Applicable | 0                         | _Not included_  |
| -1   | Not Specified  | 0                         | High            |
| 1    | Absent         | 7.17                      | Absent          |
| 2    | Present        | 83.21                     | High            |
| 3    | Ambiguous      | 0                         | Poor            |
| 4    | Trace          | 0.31                      | Poor            |
| 5    | Weak           | 4.05                      | Poor            |
| 6    | Moderate       | 1.16                      | High            |
| 7    | Strong         | 4.08                      | High            |
| 8    | Very strong    | 0.03                      | High            |

* **WormBase**

| quality   | percentage in source file | quality in Bgee |
|-----------|---------------------------|-----------------|
| _empty_   | 33.64                     | High            |
| certain   | 62.01                     | High            |
| partial   | 3.17                      | High            |
| uncertain | 1.18                      | Poor            |

* **ZFIN**

| Thisse stars | quality                                                                                                         | percentage in source file | quality in Bgee |
|--------------|-----------------------------------------------------------------------------------------------------------------|---------------------------|-----------------|
| 1 star       | Probe is difficult to use. General basal level of expression with more intense labeling in particular structure | 0.34                      | Poor            |
| 2 stars      | Weak expression pattern                                                                                         | 1.20                      | Poor            |
| 3 stars      | Moderate expression pattern                                                                                     | 1.14                      | High            |
| 4 stars      | Nice strong expression pattern                                                                                  | 0.85                      | High            |
| 5 stars      | Simple to use, intense expression pattern restricted to a few structures                                        | 0.63                      | High            |
| NO star      | (Experiments not made by Thisse)                                                                                | 95.83                     | High            |


## WormBase data info
With the file coming from WormBase AceDB ! Should be fixed with WormMine.




Include data from WormBase. They provide their _in situ_ data to us here: [ftp://caltech.wormbase.org/pub/wormbase/expr_dump/](ftp://caltech.wormbase.org/pub/wormbase/expr_dump/). Note about their _in situ_ hybridizations, from Daniela Raciti:


"Generally, we capture the expression patterns without keeping a stage specific identity. I know it's not ideal but it was done historically like that and change it would mean revise all the previous annotations (I guess at some point we might as well do it). For example if the authors say: 'gene xyz is expressed in the pharynx at L1 stage and in neurons from embryo to adult' then the annotation would be pharynx and neuron in embryo, larva and adult. We loose stage specific info. Even when there is only one stage, the authors could still have said: 'gene xyz is expressed in the pharynx, vulva, and intestine. At stage L1 expression is observed in the nerve ring.' We might capture L1 for nerve ring but then you have all the other stuff. How about including the ones that have only one anatomical term and one -or multiple- life stages. That should be safe" (file on the git, `extra/pipeline/In_situ/WormBase/expr_pattern.ace.20140526`)


Note about their strains, from Daniela Raciti: "in the file that you have is normally just based on the N2 strain (wild type). Expression in mutant background is captured via gene regulation and is not included here. There are some cases in which the strain used was different -as the example you picked- but in this case contains an extra-chromosomal array to express cye-1. The genotype is still WT. It is generally safe to say that all the data in the file are good to be displayed on the WT. I can provide you with a list of strains though -attached" (file on the git, `extra/pipeline/In_situ/WormBase/Strains.ace`)
