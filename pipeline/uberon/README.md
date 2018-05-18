**Requirements**: having successfully run the steps in [database creation/](../db_creation/) and [species/](../species/).

**Goal**: insert Uberon into the database. Requires to generate our own taxonomy ontology (already done for the `species` step), to generate taxon constraints for each terms in Uberon and each species in Bgee, to produce our own dev. stage ontology, and to extract a custom anatomical ontology for Bgee from Uberon.

## Details

### Species

* The file `source_files/species/bgeeSpecies.tsv` contains the species used in Bgee. This is the file to modify to add/remove a species.

### Makefile.taxon_info

Overall these pipeline steps, you will edit the file [pipeline/Makefile.taxon_info](../Makefile.taxon_info) for overriding taxon constraints. If new species were added to Bgee:

* you need to be extra careful in redefining the species members of the various taxa defined in this file. E.g., `BOREOEUTHERIA_IDS`, `GLIRES_IDS`.
* You need to update the variables `OVERRIDING_TAXA` and `TAXA_VAR_PREFIXES` approriately.
* You might need to update the variable `SIMPLIFICATION_STEPS`, used to generate taxon constraints in several steps (otherwise, they cannot be computed for some taxa, due to implementation issues in Elk reasoner)

### Taxon constraints (TODOs before pipeline run)

* You need to manually define taxon constraints, that will override some taxon constraints automatically generated using Uberon. See the variable `OVERRIDE_TAXON_CONSTRAINTS` in [pipeline/Makefile.taxon_info](../Makefile.taxon_info). This is because some are incorrect. And because, as of bgee 14, it is not possible to extract taxon constraints from the `composite-metazoan` version, so we need to manually define taxon constraints for species-specific ontologies.
  * First, define taxon constraints for each species-specific dev. stage ontologies, check if new species were added, see `DEV_STAGE_ONT_TAXON_CONSTRAINTS`.
  * Then you might need to edit `DEV_STAGE_FIX_TAXON_CONSTRAINTS` to fix incorrect taxon constraints in the dev. stage ontology.
  * Then, define taxon constraints for each taxon-specific anatomical ontologies, see `ANAT_TAXON_CONSTRAINTS`. To do so:
    * Check the declaration `treat-xrefs-as-reverse-genus-differentia` at the beginning of Uberon (e.g., visible here http://www.ontobee.org/ontology/UBERON), check if new ontologies were added.
    * Extract all species-specific terms from the composite-metazoan ontology and check if an ontology was not declared, e.g. `grep "\-\->" source_files/uberon/composite-metazoan.owl | grep -v UBERON_ | grep -v NCBITaxon_ | sort | uniq`. Check if new ontologies were added.
  * Then, potentially fix incorrect taxon constraints in the taxonomy by editing `ANAT_FIX_TAXON_CONSTRAINTS`.
  * Of note, when working on dev. stage ontology, it is more convenient to generate dev. stage constraints only on the ontology provided by Chris Mungall for this work, as only "errors" will needs to be overridden in the command line, rather than to use the constraints generated in `generated_files/uberon/taxonConstraints.tsv`.
  * Then it is convenient to merge this dev. stage taxon constraints with our taxon constraints generated from Uberon, to have one unique taxonConstraints file for our pipeline.

### Dev. stage ontology (TODOs before pipeline run)

* The composite-metazoan dev. stage ontology is now built separately on the related tracker, see https://github.com/obophenotype/developmental-stage-ontologies/tree/master/external

* Chris Mungall generates a "dev. stage only" composite-metazoan ontology once we are happy with our species-specific dev. stage ontologies.

* Then, from this ontology, we extract our own version to be later manually tuned. From the `Makefile` in this folder, the command used was:
    ```
    @$(JAVA) UberonDevStage generateStageOntology $(UBERON_COMPOSITE_FILE_PATH) $(UBERON_OUTPUT_PATH)$(DEV_STAGE_ONT_PREFIX) UBERON:0000067$(LIST_SEP)UBERON:0000071$(LIST_SEP)UBERON:0000105$(LIST_SEP)UBERON:0000000$(LIST_SEP)FBdv:00007008$(LIST_SEP)MmusDv:0000041$(LIST_SEP)ZFS:0000000$(LIST_SEP)UBERON:0035944$(LIST_SEP)UBERON:0035945$(LIST_SEP)BFO:0000015 - UBERON:0000481$(KEY_VAL_SEP)NCBITaxon:6072 BFO:0000050$(LIST_SEP)BFO:0000062$(LIST_SEP)RO:0002087 UBERON:0000104$(LIST_SEP)WBls:0000075$(LIST_SEP)FBdv:00005259$(LIST_SEP)ZFS:0100000$(LIST_SEP)XAO:1000094$(LIST_SEP)NCBITaxon:1
    ```

* Manually add all xrefs of the term `UBERON:0000105` in the original ontology to the term `UBERON:0000104` in our custom ontology. **TODO**: make `generateStageOntology` to do this automatically.

* Then, the extracted dev. stage ontology needs to be manually finely tuned, using a command that needs to be run multiple times to detect errors, until there is no error left, e.g.:
    ```
    java -Xmx{RAM_GB_VALUE}g -Dbgee.dao.jdbc.username={ROOT_MYSQL_LOGIN} -Dbgee.dao.jdbc.password={ROOT_MYSQL_PASS} -Dbgee.dao.jdbc.driver.names=com.mysql.jdbc.Driver -Dbgee.dao.jdbc.url=jdbc:mysql://127.0.0.1:3306/mysql -jar bgee-pipeline-14-with-dependencies.jar InsertUberon insertStages dev_stage_ont.obo taxonConstraints.tsv FBdv://7227,XAO://8364,WBls://6239,ZFS://7955,UBERON:0000069://7227--7955 NCBITaxon:1 species.tsv
    ```
    (where `species.tsv` can correspond to `source_files/species/bgeeSpecies.tsv` and `taxonConstraints.tsv` to `generated_files/uberon/taxonConstraints.tsv`)

  * To do so, it is necessary to add `immediately_preceded_by` and `preceded_by` relations.
  * And to override taxon constraints in the command line when they are incorrect. For an explanation about the constraints overridden, see https://github.com/obophenotype/developmental-stage-ontologies/tree/master/external/bgee/taxonConstraintsOverridden.md. You need to edit this file whenever you change the constraints overriden.
  * Of note, it is more convenient to generate the taxon constraints file `taxonConstraints.tsv` on the ontology provided by Chris Mungall for this work, as only "errors" will needs to be overridden in the command line, rather than to use the constraints generated in `generated_files/uberon/taxonConstraints.tsv`.
  * Then it is convenient to merge this dev. stage taxon constraints with our taxon constraints generated from Uberon, to have one unique taxonConstraints file for our pipeline.

### Anatomical ontology (TODOs before pipeline run)

You need to generate a custom version of the `composite-metazoan` ontology, for fixing errors, for removing terms we're not interested in, and for improving the display of the graph (although as of Bgee 13 we never display the complete ontology, so this is not a priority for now)

* First, you might want to generate a version with only is_a/part_of-like relations. This will allow you to immediately spot terms with missing relations, roots of subgraphs to remove, and roots of subgraphs to keep. Can be done by editing the Makefile, or by running a command like:
    ```
    java -Xmx32g -jar bgee-pipeline-VERSION-with-dependencies.jar Uberon simplifyUberon composite-metazoan.owl custom_uberon simplification_composite_info.tsv - BFO:0000050 - - - - -
    ```

* Then open the file in OBO-Edit (in Protege, part_of relations do not appear in the graph, so a lots of roots will appear that simply miss a is_a relation)

    1. Custom modifications to source ontology. As of Bgee 14:
        * delete `EHDAA2:0000000` and add it as a xref to `UBERON:0001062` anatomical entity

    2. Identify terms that you simply want to go away, with their incoming edges merged with their outgoing edges. Their subgraph is conserved. Be careful about important XRefs that might disappear! As of Bgee 14, they are:
        * `CARO:0000006 material anatomical entity` (redundant with `UBERON:0000465`)
        * `UBERON:0010000 multicellular anatomical structure` (what does it add for us as compared to "anatomical structure"?)
        * `UBERON:0008979 carcass` (we don't annotate such terms)
        * `GO:0016265 obsolete death` (seems that some of its children are not obsolete)
        * `BFO:0000030 object` (maybe there are still some useful terms in the subgraph, so we don't remove the whole subgraph)
        * `BFO:0000020 specifically dependent continuant`
        * `BFO:0000004 independent continuant`
        * `BFO:0000002 continuant` (all these continuants will make "anatomical entity" to become a root)
        * `UBERON:0001062 anatomical entity`
        * `BFO:0000040 material entity`

    3. Identify relations that you want to remove. This is especially important if you want a nice display, a bit less as long as we don't display the graph. As of Bgee 14:
        * `UBERON:0010707 appendage girdle complex` in_lateral_side_of `UBERON:0000468 multi-cellular organism` (what does it mean exactly?)
        * `UBERON:0000463 organism substance` part_of `UBERON:0000468 multi-cellular organism` (is annotated as to be maybe removed)
        * `UBERON:0000481 multi-tissue structure` part_of `NCBITaxon:6072 Eumetazoa` (weird, should be "only_in_taxon" Metazoa I guess)
        * `UBERON:0005212 Leydig cell region of testis` part_of `NCBITaxon:7742 Vertebrata` (check whether an "only_in_taxon" relations was added, and then stop removing this relation)

    4. Identify the relation types you want to keep in the ontology. All their sub-relations will be mapped to them. As of Bgee 14:
        * `BFO:0000050 part_of`
        * `RO:0002202 develops_from`
        * `RO:0002494 transformation_of`

    5. Identify subgraphs you absolutely don't want to be part of the ontology (terms related by is_a/part_of/develops_from/transformation_of relations, as these are the only relations we keep before filtering). As of Bgee 14:
        * `NBO:0000313 behavior process` (not anatomy)
        * `CEPH:0000300 cephalopod trait` (not anatomy)
        * `CHEBI:24431 chemical entity` (not anatomy)
        * `PATO:0000001 quality` (not anatomy)
        * `BFO:0000019 quality` (not anatomy)
        * `ENVO:01000254 environmental system` (not anatomy)
        * `GO:0003674 molecular_function` (not anatomy)
        * `CHEBI:23367 molecular entity` (not anatomy)
        * `CL:0000063 obsolete cell by histology` (regroups obsolete cell terms)
        * `UBERON:0000952 obsolete adductor muscle` (regroups obsolete cell terms)
        * `WBls:0000075 Nematoda Life Stage` (life stages are extracted in a different ontology)
        * `BFO:0000003 occurrent` (life stages are extracted in a different ontology)
        * `FBdv:00007008 occurrent` (life stages are extracted in a different ontology)
        * `CS:0 Carnegie stage` (life stages are extracted in a different ontology)
        * `ZFS:0100000 stages` (life stages are extracted in a different ontology)
        * `OlatDv:0000010 developmental stage` (life stages are extracted in a different ontology)
        * `UBERON:0000466 immaterial anatomical entity` (should never be used for gene expression annotation...)
        * `PORO:0000019 immaterial anatomical entity` (should never be used for gene expression annotation...)
        * `BFO:0000141 immaterial anatomical entity` (should never be used for gene expression annotation...)
        * `PORO:0000923 cell component` (redundant and almost empty)
        * `PORO:0000963 spiculogenesis` (seems to be a process?)
        * `SO:0000704 gene`
        * Do **NOT** select: `NCBITaxon:1 root` (otherwise, all the GCI relations will go away)

    6. Identify subgraphs you want to conserve in the ontology (terms related by is_a/part_of/develops_from/transformation_of relations, as these are the only relations we keep before filtering). This allows to conserve classes that are part of a subgraph to remove, but also of subgraph to keep. So this selection is extremely important. Parents of these terms will also be kept. As of Bgee 14, we select `UBERON:0000465 material anatomical entity`, as well as the Cell Ontology:
        * `UBERON:0000465 material anatomical entity`
        * `GO:0005575 cellular_component`
        * Do **NOT** select: `UBERON:0001062 anatomical entity` (otherwise "immaterial anatomical entity" will be kept)
        * Do **NOT** select: `NCBITaxon:1 root` (otherwise all terms with a GCI relation will be kept)

    7. Define subsets for which you want all their incoming edges to be removed, unless a term would be left orphan. This is to allow a nicer display of the graph. As of Bgee 14:
        * upper_level "abstract upper-level terms not directly useful for analysis", e.g., `processual entity`, `anatomical surface`, `anatomical wall`.
        * You can specify classes that should be excluded from this edge removal. As of Bgee 14:
          * `UBERON:0013701 main body axis`
          * `UBERON:0000026 appendage`
          * `UBERON:0000480 anatomical group`
          * `FBbt:00007276 anatomical group (Drosophila)`
          * `UBERON:0000479 tissue`
          * `UBERON:0011676 subdivision of organism along main body axis`
          * `UBERON:0005423 developing anatomical structure`
          * `UBERON:0000463 organism substance`
          * `GO:0005575 cellular_component`

    8. While you run the pipeline, you will have to do a lot of remapping manually: there are lots of inconsistencies or synchronization issues between species-specific ontologies or MODs and Uberon. See https://github.com/obophenotype/uberon/issues/664, and especially https://github.com/obophenotype/uberon/issues/1288 and https://gitlab.isb-sib.ch/Bgee/bgee_pipeline/issues/68.
    
    * The overall command line as of Bgee 14 is:
    ```
    @$(JAVA) Uberon simplifyUberon $(UBERON_COMPOSITE_FILE_PATH) $(UBERON_OUTPUT_PATH)$(CUSTOM_UBERON_PREFIX).tmp $(UBERON_OUTPUT_PATH)simplification_composite_info.tsv CARO:0000006$(LIST_SEP)UBERON:0010000$(LIST_SEP)UBERON:0008979$(LIST_SEP)GO:0016265$(LIST_SEP)BFO:0000030$(LIST_SEP)BFO:0000020$(LIST_SEP)BFO:0000004$(LIST_SEP)BFO:0000002$(LIST_SEP)UBERON:0001062$(LIST_SEP)BFO:0000040 BFO:0000050$(LIST_SEP)RO:0002202$(LIST_SEP)RO:0002494 UBERON:0010707$(KEY_VAL_SEP)UBERON:0000468$(ENTRY_SEP)UBERON:0000463$(KEY_VAL_SEP)UBERON:0000468$(ENTRY_SEP)UBERON:0000481$(KEY_VAL_SEP)NCBITaxon:6072$(ENTRY_SEP)UBERON:0005212$(KEY_VAL_SEP)NCBITaxon:7742 NBO:0000313$(LIST_SEP)CEPH:0000300$(LIST_SEP)CHEBI:24431$(LIST_SEP)PATO:0000001$(LIST_SEP)BFO:0000019$(LIST_SEP)ENVO:01000254$(LIST_SEP)GO:0003674$(LIST_SEP)CHEBI:23367$(LIST_SEP)CL:0000063$(LIST_SEP)UBERON:0000952$(LIST_SEP)WBls:0000075$(LIST_SEP)BFO:0000003$(LIST_SEP)FBdv:00007008$(LIST_SEP)CS:0$(LIST_SEP)ZFS:0100000$(LIST_SEP)OlatDv:0000010$(LIST_SEP)UBERON:0000466$(LIST_SEP)PORO:0000019$(LIST_SEP)BFO:0000141$(LIST_SEP)PORO:0000923$(LIST_SEP)PORO:0000963$(LIST_SEP)SO:0000704 UBERON:0000465$(LIST_SEP)GO:0005575 upper_level UBERON:0013701$(LIST_SEP)UBERON:0000026$(LIST_SEP)UBERON:0000480$(LIST_SEP)FBbt:00007276$(LIST_SEP)UBERON:0000479$(LIST_SEP)UBERON:0011676$(LIST_SEP)UBERON:0005423$(LIST_SEP)UBERON:0000463$(LIST_SEP)GO:0005575
    ```

#### After generation

* Then, it's a good idea to look at our custom ontology with only is_a/part_of relations, and by removing GCI relations to see roots of all species. With a command such as:
    `java -Xmx32g -jar bgee-pipeline-VERSION-with-dependencies.jar Uberon simplifyUberon custom_uberon.owl custom_uberon_isapartof simplification_composite_info.tsv - BFO:0000050 - NCBITaxon:1 - - -`

* And also to look precisely at the file generated `generated_files/uberon/simplification_composite_info.tsv` (see "Data verification" below).

* Add to the root of the ontology xrefs to the "unspecified" and "unknown" terms (they should have been removed by graph filtering). Perform the modifications in `generated_files/uberon/custom_composite.obo` and `generated_files/uberon/custom_composite.owl`, e.g., add xrefs to `XAO:0003003` and `ZFA:0001093` to `UBERON:0000465 material anatomical entity`.

* If the query allowing to check inconsistent taxon constraints between anatomical entities and anatomical entity relations returns some results (it should be none), it's possible to use this command:
	```
	java -Dbgee.dao.jdbc.username={ROOT_MYSQL_USERNAME} -Dbgee.dao.jdbc.password={ROOT_MYSQL_PASS} -Dbgee.dao.jdbc.driver.names=com.mysql.jdbc.Driver -Dbgee.dao.jdbc.url=jdbc:mysql://{SERVER_NAME}:3306/bgee_v{BGEE_VERSION} -jar bgee-pipeline-{BGEE_VERSION}-with-dependencies.jar CorrectTaxonConstraints fromRDB`
	```
	It will correct all inconsistencies. Please check WARNING output messages. It's sometime necessary to execute some mysql commands manually.
	This tool has been developed by J. Wollbrett in order to correct inconsistencies present in the RDB of Bgee 14 Beta. Good practices would be to directly correct the ontology...
	
### Cycle removal

It is needed to remove some cycles between terms of Uberon. The script [pipeline/uberon/check_cycles.pl](check_cycles.pl)
is very useful for this task. It will identify cycles, and the paths causing the cycles,
as well as the most represented paths among all the cycles.

**Important**: Note that, as of Bgee 14, all the cycles are caused for this reason: https://github.com/obophenotype/uberon/issues/651 (look also at comments of this issue).

* Run the step, e.g., from [pipeline/uberon/](.), `make ../../generated_files/uberon/detect_cycles`,
and check carefully the file `generated_files/uberon/detect_cycles`, notably the end of the files,
listing the most represented cycles and most represented subgraphs in cycles. Note that currently
the script only checks 'is_a part_of' relations, but it is easy to change in the script.

* Check whether you identified new subgraphs as compared to the cycles identified in https://github.com/obophenotype/uberon/issues/651
(look also at the comments). You need to remove each of these already-identified relations causing cycles.
To do so, you can:
  * Find the appropriate relation with a command such as (the last column allows you the see
  the taxa in which the relation is valid in. A `NULL` value means that it is valid in all taxa):
  ```
  SELECT t1.anatEntityRelationId, t1.anatEntitySourceId, t3.anatEntityName, t1.relationType, t1.anatEntityTargetId, t4.anatEntityName, GROUP_CONCAT(t2.speciesId ORDER BY t2.speciesId) FROM anatEntityRelation AS t1 INNER JOIN anatEntityRelationTaxonConstraint AS t2 ON t1.anatEntityRelationId = t2.anatEntityRelationId INNER JOIN anatEntity AS t3 ON t1.anatEntitySourceId = t3.anatEntityId INNER JOIN anatEntity AS t4 ON t1.anatEntityTargetId = t4.anatEntityId WHERE t1.relationStatus != 'reflexive' and t1.anatEntitySourceId = 'UBERON:0000074' and t1.anatEntityTargetId = 'UBERON:0005325' GROUP BY t1.anatEntityRelationId;
  ```
  * Delete the relation with a command such as:
  ```
  delete from anatEntityRelation where anatEntityRelationId = XXX
  ```

* Once you will have remove all direct relations causing cycles, there will still be cycles
detected thanks to the remaining indirect relations in the database. But there will be no subgraph
identified to explain them. So, you just need to now identify which indirect relations
to delete from the database. To do so, you can:
  * Find the appropriate relation with a command such as (the last column allows you the see
  the taxa in which the relation is valid in. A `NULL` value means that it is valid in all taxa;
  it is very informative since most likely the relation to delete is valid in only one species):
  ```
  SELECT t1.anatEntityRelationId, t1.anatEntitySourceId, t3.anatEntityName, t1.relationType, t1.anatEntityTargetId, t4.anatEntityName, GROUP_CONCAT(t2.speciesId ORDER BY t2.speciesId) FROM anatEntityRelation AS t1 INNER JOIN anatEntityRelationTaxonConstraint AS t2 ON t1.anatEntityRelationId = t2.anatEntityRelationId INNER JOIN anatEntity AS t3 ON t1.anatEntitySourceId = t3.anatEntityId INNER JOIN anatEntity AS t4 ON t1.anatEntityTargetId = t4.anatEntityId WHERE t1.relationStatus != 'reflexive' and t1.anatEntitySourceId IN ('UBERON:0011299', 'UBERON:0019261') and t1.anatEntityTargetId IN ('UBERON:0011299', 'UBERON:0019261') GROUP BY t1.anatEntityRelationId;
  ```
  * Then again, delete the relation with a command such as:
  ```
  delete from anatEntityRelation where anatEntityRelationId = XXX
  ```

* Once you are happy with the deletion, re-run the step to identify cycles. There should be none.

### Taxonomy

* This pipeline step requires the NCBI taxonomy, provided as an ontology.
  * We cannot use the [official taxonomy ontology](http://www.obofoundry.org/cgi-bin/detail.cgi?id=ncbi_taxonomy) because to correctly infer taxon constraints at later steps, we need this ontology to include disjoint classes axioms between sibling taxa, as explained in a Chris Mungall [blog post](http://douroucouli.wordpress.com/2012/04/24/taxon-constraints-in-owl). The default ontology does not include those.
  * This custom taxonomy ontology was generated at a previous step, see [species step](../species/).
  * This custom taxonomy should include the species used in Bgee and their ancestors, the taxa used in our annotations and their ancestors, the taxa used in Uberon an their ancestors. To extract taxa used in Uberon, we use the `ext` version (this is the one containing more taxa).

## Data generation

* If it is the first time you execute this step in this pipeline run:
  `make clean`

* Run Makefile:
  `make`

* Modify the file [pipeline/db_creation/update_data_sources.sql](../db_creation/update_data_sources.sql): you need to add the last modification date of the uberon version used, and update the URL.

* Release our updated dev. stage ontology, along with formal reports, see https://github.com/obophenotype/developmental-stage-ontologies/tree/master/external/bgee/dev_stage_ontology.obo, https://github.com/obophenotype/developmental-stage-ontologies/tree/master/external/bgee/report.md, and https://github.com/obophenotype/developmental-stage-ontologies/tree/master/external/bgee/known_issues.md

## Data verification

* Check very carefully the file `generated_files/uberon/simplification_composite_info.tsv`. This contains information about all terms that were deleted from Uberon following our simplification, with information about the reason for deletion. You need to absolutely make sure that all classes were intentionally deleted.

* Check very carefully the file `generated_files/uberon/insert_stages`.
  * In this file you notably have a nicely-formatted report of dev. stages structure for each species
  * You also have a table listing for each UBERON term the species the stage exists in. For each of this term, you need to check that we correctly have an xref corresponding to each of these species. So, if a term is supposed to exist in all Bgee species, we should have 29 xrefs to other dev. stage ontologies (as of Bgee 14). Be careful that not all xrefs point to the proper ontologies, so it is not enough to simply counting xrefs. Check that each xref maps to the appropriate term.

* `generated_files/uberon/step_verification_RELEASE.txt` should contain:

* The following files should have been downloaded/generated:
  * `source_files/uberon/ext.owl`
  * `source_files/uberon/composite-metazoan.owl`
  * `source_files/uberon/dev_stage_ontology.obo`
  * `generated_files/uberon/taxonConstraints.tsv`
  * `generated_files/uberon/simplification_composite_info.tsv`
  * `generated_files/uberon/custom_composite.obo`
  * `generated_files/uberon/custom_composite.owl`
  * `generated_files/uberon/insert_stages`
  * `generated_files/uberon/insert_anatomy`
  * `generated_files/uberon/uberon_sex_info.tsv`

## Following step
  *  See post_processing step

## Error handling

* You can have an exception thrown, saying that a specified taxon does not exist in the taxonomy ontology, for instance `java.lang.IllegalArgumentException: Taxon NCBITaxon:71164 was not found in the ontology`. It likely means that an incorrect/deprecated taxon is used in Uberon. Remove the ID of the taxon in the file `generated_files/species/allTaxIds.tsv` (so if the exception is related to a taxon `NCBITaxon_71164`, remove the ID `71164`). Check if the taxon ID is present in the file `generated_files/species/annotTaxIds.tsv`, if it is, the error is on us. Otherwise, report the problem on the Uberon tracker, if you identified the taxon in Uberon. You will most likely need to manually modify the ontology to remove the offending taxa in the mean time.

## Other notable Makefile targets

* To delete the anatomical ontology inserted into database:
    `make delete_anatomy_from_db`

* To delete the dev. stage ontology inserted into database:
    `make delete_stage_from_db`

* Generate taxon constraints:
    `make ../../generated_files/uberon/taxonConstraints.tsv`

* Generate custom Uberon version:
    `make ../../generated_files/uberon/custom_composite.owl`

* Generate sex-related info about anatomical terms:
    `make ../../generated_files/uberon/uberon_sex_info.tsv`
