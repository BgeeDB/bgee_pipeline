# Operations on annotations/conditions

1. [Deletion of unused conditions](#deletion-of-unused-conditions)
2. [Condition remapping](#condition-remapping)

## Deletion of unused conditions

Use the SQL script `delete_unused_conditions.sql` present in this folder.


## Condition remapping

Sometimes, it is needed to update the database after data insertion, to update a condition.
For instance, in case we spot an incorrect annotation, not supposed to exist
in the targeted species.

Procedure to update conditions:

### 1) Update the mapping file

The mapping file to edit is [source_files/annotations/condition_remapping.tsv](../../source_files/annotations/condition_remapping.tsv). Warning, the script using this file are based on
the column order for now, not the column names. Don't change the order.

* The first column is the `conditionId` of the condition to be remapped.
* The second column is the anat. entity ID to be remapped to. Should be the same
as the original one if you don't want to change this value.
* Third column is the dev. stage ID to be remapped to. Should be the same
as the original one if you don't want to change this value.
* Fourth column is the sex value to be remapped to. Should be the same
as the original one if you don't want to change this value.
* Fifth column is the strain value to be remapped to. Should be the same
as the original one if you don't want to change this value.

### 2) Run the script to identify or create the conditions to remap to

In this folder, run the Makefile target `make ../../generated_files/annotations/remap_cond`.

This script with only populate the table `remapCond` in the database. Column `incorrectConditionId`
is the condition ID that needed to be remapped. `remappedConditionId` is the ID of the condition
to map to.

### 3) Check that the mapping is what you expected

Run the following SQL commands.

* One for checking that you didn't unintentionally change some values you did not want to.
For instance, if you only wanted to update the dev. stage ID, check that the anat. entity ID,
sex, strain, sexInferrence, are the same between the condition to remap and the condition to map to
(should return no result).
For instance:

```
SELECT t1.*, t2.anatEntityId AS oldAnatEntityId, t3.anatEntityId AS remappedAnatEntityId,
       t2.sex AS oldSex, t3.sex AS remappedSex,
       t2.strain AS oldStrain, t3.strain AS remappedStrain,
       t2.sexInferred AS oldSexInferred, t3.sexInferred AS remappedSexInferred
FROM remapCond AS t1
INNER JOIN cond AS t2 ON t1.incorrectConditionId = t2.conditionId
INNER JOIN cond AS t3 ON t1.remappedConditionId = t3.conditionId
WHERE t2.anatEntityId != t3.anatEntityId OR
      t2.sex != t3.sex OR
      t2.strain != t3.strain OR
      t2.sexInferred != t3.sexInferred;
```

* One for checking that you didn't miss updating a field of a condition to remap, so that
the condition to remap is remapped to itself (should return no result).

```
SELECT * FROM remapCond WHERE incorrectConditionId = remappedConditionId;
```

### 4) Update the conditions you wanted to modify

For instance, if you wanted to change annotation of RNA-Seq libraries, and only RNA-Seq libraries:

```
UPDATE t1 FROM rnaSeqLibrary AS t1
INNER JOIN remapCond AS t2 ON t1.conditionId = t2.incorrectConditionId
SET t1.conditionId = t2.remappedConditionId;
```

Do this operation for every table containing annotations you want to update
(`affymetrixChip`, `estLibrary`, `inSituSpot`, ...).

### 5) Delete unused conditions

Follow [the procedure described in this README](#deletion-of-unused-conditions)
to delete unused conditions.

### 6) Empty the table `remapCond`

Do that when you are completely done, it's not a big deal to not empty it,
but you might run into surprises if you do another update later.