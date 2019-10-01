# Insert Affymetrix

* **Requirements**: having successfully run the step database creation AND the step initialization.
* **Goal**:         insert Affymetrix data used in Bgee.

## Details
* This step will download newly annotated Affymetrix experiences (new cel and new mas5) from annotator's machine
```
make new_annotation
```

  * Once downloaded, **remove them from the annotators' computers to not re-process them next time! (Currently Anne's machine `annotbioinfo`)**

  * Currently the rsync command used at this step uses SSH and the step runner **has to** cat its ssh public key on annotator's machine at the end of the `$(ANNOTATORLOGIN)HOME/.ssh/authorized\_keys` file. See http://troy.jdmz.net/rsync/ about how to do this. Possible solutions could be to use rsync through a rsync daemon on annotator's machine, but it may cause file-system problems on Mac; or to log on annotator's machine through a guest account, without password. **Investigating** ...

* It will compute the quality scores, and check for duplicated chips.
  * Check for absence of conflicts between new files, and already present files, and move new files that are OK to the directory where files are permanently stored
  ```
  make maintenance
  ```

  * MAS5 files need to be standardized
  ```
  make check_mas5
  ```

* It will generate information about experiments, chip, ...
```
make gen_info
```

* It will check annotations
```
make check_annot
```
  * **Problem! Experiment not inserted in `microarrayExperiment` or experiment annotated twice in `affymetrixChip`**. Warning can most of the time be skipped by sorting the `source_files/Affymetrix/affymetrixChip.tsv` file: sort based on column B then column A.
  * If you need to remove an experiment (e.g. gone from ArrayExpress), use the command
  ```
  perl Maintenance/delete_affy_experiment.pl $(AFFYDATAPATH) $(AFFY_CHIPINFO_FILEPATH) <ExpId_to_delete>
  ```

* It will normalize Affymetrix analyses
```
make normalization
```

* It will do final checks before db insertion
```
make check_after
```

* It will insert Affymetrix data in the database
```
make insert_affy
```

* It will check validity of inserted conditions: before running `insert_expression`, you should generate the file `check_conditions`, to detect invalid conditions not supposed to exist in the related species. See 'Details' section of [pipeline/post_processing/README.md](../post_processing/README.md) for an explanation on how to fix such issues (in case the annotations were not incorrect).
```
make check_conditions
```

* It will insert expression data in the database
```
make insert_expression
```


## Data generation
* If it is the first time you execute this step in this pipeline run:
```
make clean
```

* Run Makefile:
```
make
```


## Error handling
* If there are discrepancies between cel chipId extension case between `source_files/Affymetrix/affymetrixChip.tsv` and [pipeline/Affymetrix/affymetrixChipInformation](affymetrixChipInformation) files, there is the script [Maintenance/fix_cel_extension_name.pl](Maintenance/fix_cel_extension_name.pl) to help fixing that.

* If you see the following warning in the `check_annot` make step "**Problem! Experiment not inserted in microarrayExperiment or experiment annotated twice in affymetrixChip**", most of the time you can fix it by sorting the `source_files/Affymetrix/affymetrixChip.tsv` file: sort based on column B (`Experiment ID`) then column A (`Chip ID`).


## General remarks
* Because many folders (for cel files) downloaded in ArrayExpress finish by ".raw" or ".raw.1", you can rename them using:
```
for i in *.raw;       do mv $i ${i/.raw/};       done
for i in *.raw.[0-9]; do mv $i ${i/.raw.[0-9]/}; done
```

* Data files (".cel") and MAS5 are stored on `devbioinfo` and are not on the git because too BIG -> find a better solution? They are mirrored on `annotbioinfo` and Fred's machines.


## Check new chipType(s) supported by Ensembl
[Affymetrix/Maintenance/test_Ensembl_Chiptype_Support.pl](Maintenance/test_Ensembl_Chiptype_Support.pl) allows to check chipType(s) support for a particular species in Ensembl.


Ensembl chipType name is very often different from ArrayExpress name. So double check!!


## Clean up MAS5 files
```
make check_mas5
```
MAS5 files can have various formats, number of columns, calls correspondences, etc... [clean_mas5_files.pl](MAS5/clean_mas5_files.pl) is a script to extract the relevant information from the files, and generate clean and standard new MAS5 files. It will read MAS5 files from `$(MAS5ORIPATH)`, and put the cleaned files in `$(MAS5PATH)`, with exactly the same names, directories tree, etc.

You need to be in [pipeline/Affymetrix/](.) to run the script.


### Clean all files
The first usage of this script is to clean all files. Already-cleaned files will not be deleted, nor overwritten. Only files not already cleaned will be. If you want to clean again all files, remove all files and directories in `$(MAS5PATH)` (it is harmless, as all original files are kept anyway).

* Usage:
    ```
    perl MAS5/clean_mas5_files.pl $(MAS5ORIPATH) $(MAS5PATH)
    ```
    * `$(MAS5ORIPATH)`: path where not-already cleaned mas5 files are stored
    * `$(MAS5PATH)`: path to the directory where to store cleaned mas5 files

* Classical error messages will be for instance:
    ```
    Error, invalid mas5 file format: $(MAS5ORIPATH)E-MEXP-10/NT4WT3C_Norm
    ```
    This means that the script couldn't detect proper columns in the file: probeset_id, MAS5 call, signal intensity. You need to manually relaunch the script for this file (see below)

* Classical warning messages, for instance:
    ```
    Warning, unrecognized line, # 500, in E-MEXP-10/NT4WT3C_Norm: [...line displayed here]
    ```
    It means that the line was not recognized. Might simply be an additional header line for instance. Check what went wrong, and if needed, add manually the line into the new generated file in `$(MAS5PATH)`

* And also:
    ```
    Warning, unrecognized header, in E-MEXP-10/NT4WT3C_Norm: [...header displayed here]
    ```
    It means that it was not possible to identify proper columns from the header of the file. If it is not followed by the "Error, invalid mas5 file format", it means that it was possible to guess the proper columns from their content. It is then harmless, and the message is displayed for your information, and eventually to modify the function `get_mas5_columns_from_header` in [mas5_utils.pl](mas5_utils.pl) (but do not be too permissive, maybe it is better to keep a few warning messages).


### Check that everything went fine
**Run a diff command on directories to ensure that all files were cleaned.**

From [pipeline/Affymetrix/](.):
```
diff -rq $(MAS5ORIPATH) $(MAS5PATH) | grep -v ' diff' | grep -v '\.svn' | grep -v 'not_separated'
perl MAS5/check_mas5_filtered.pl -affyChipFilesDir=$(AFFYDATAPATH)  -affymetrixChip=$(PIPELINEROOT)$(AFFY_CHIP_FILEPATH)
```
It will tell you if invalid files are present.


### If you need to rerun the script manually for one file
You can run the script for one file only, and specifying the relevant column indexes (rank of the column containing probeset IDs, MAS5 calls, signal intensities). The clean file will be overwritten even if it already exists.

Usage:
```
perl MAS5/clean_mas5_files.pl $(MAS5ORIPATH) $(MAS5PATH) <file path> <probeset_id column> <call column> <signal column>
```
* `$(MAS5ORIPATH)`: path where not-already cleaned mas5 files are stored
* `$(MAS5PATH)`: path to the directory where to store cleaned mas5 files
* `<file path>`: path to an original mas5 file, RELATIVE TO `<path to mas5 original files>` (basically, 'expId/chipId')
* `<probeset_id column>`: index of the column containing probeset IDs in the file. Usually, it is 0.
* `<call column>`: index of the column containing the mas5 calls in the file. Usually, it is 1.
* `<signal column>`: index of the column containing signal intensities in the file. Usually, it is 2.

For instance, from [pipeline/Affymetrix/](.):
```
perl MAS5/clean_mas5_files.pl $(MAS5ORIPATH) $(MAS5PATH)  E-MEXP-1027/H_WT1_Norm 1 5 4
```


## Generate information
```
make gen_info
```

### Generation
This part of the pipeline generates checksums, quality score, percent present, looks into cel files to get cdf name, scan dates, etc. This is used to identified duplicated chips, and to remove low quality chips.

From [pipeline/Affymetrix/](.):
```
perl Generate_information/generate_affy_chips_information.pl -affyChipFilesDir=$(AFFYDATAPATH) -affymetrixChip=$(PIPELINEROOT)$(AFFY_CHIP_FILEPATH) \
       -affymetrixChipInformation=$(PIPELINEROOT)$(AFFY_CHIPINFO_FILEPATH) -chipTypeQual=$(PIPELINEROOT)$(AFFY_CHIPTYPEQUAL_FILEPATH)
```
This script will try to compute new information for chips present in [pipeline/Affymetrix/affymetrixChip](affymetrixChip) and not present in [pipeline/Affymetrix/affymetrixChipInformation](affymetrixChipInformation). It will then display detailed information about this new generated information: potentially duplicated chips, low quality or incompatible chips. This can be used by the annotators to define which chips they should actually annotated.

Look at all the info provided by the script very carefully.

Check that everything went fine:
```
ls -l $(AFFYDATAPATH)chip_information/logs/*.out.PROB
grep -H -c 'proc.time' $(AFFYDATAPATH)chip_information/logs/*.out | grep ':0'
```

Also check warnings:
```
grep -i "warning" -A 5 $(AFFYDATAPATH)chip_information/logs/*.out
```

Also check manually the entries added to [pipeline/Affymetrix/affymetrixChipInformation](affymetrixChipInformation):
```
git diff $(PIPELINEROOT)$(AFFY_CHIPINFO_FILEPATH) $(PIPELINEROOT)$(AFFY_CHIPTYPEQUAL_FILEPATH)
```

Also check that you have as many new information generated, as new chips moved to the server at the previous step.


### If you need to regenerate information for some chips
Remove the corresponding lines in [pipeline/Affymetrix/affymetrixChipInformation](affymetrixChipInformation). If you don't have as many new information generated, as new chips moved to the server at the previous step, it might be because the information was already present in the file. You should remove the lines to recompute the information.


### Remarks
If you've added a new chip type (incompatible chip type detected), you might need to modify the file [pipeline/Affymetrix/chipTypeCorrespondencesAndQualityThresholds](chipTypeCorrespondencesAndQualityThresholds). If some duplicated chips were detected, you need to comment in [pipeline/Affymetrix/affymetrixChip](affymetrixChip) the ones you want to remove. If there are some low quality chips, don't bother, they won't be taken into account by the pipeline.

Output of R when generating the info are in `$(AFFYDATAPATH)chip_information/logs/`. File names follow the pattern **expId_chipId.out**. Actual results for each chip are in `$(AFFYDATAPATH)chip_information/results/`. File names follow the pattern **expId_chipId.tsv**.


## Annotations checking
```
make check_annot
```
* In [pipeline/Affymetrix/](.) check that there is no problem left with the annotation:
```
perl Annotation_checking/check_affy_curation.pl  -bgee=$(BGEECMD) -normalizationType=$(PIPELINEROOT)$(AFFY_NORMTYPE_FILEPATH) \
       -detectionType=$(PIPELINEROOT)$(AFFY_DETCTYPE_FILEPATH) -chipType=$(PIPELINEROOT)$(AFFY_CHIPTYPE_FILEPATH) \
       -microarrayExperiment=$(PIPELINEROOT)$(MICROARRAY_EXPERIMENT_FILEPATH) -cel_data=$(CELPATH) -processed_mas5=$(MAS5PATH) \
       -affyChipInformation=$(PIPELINEROOT)$(AFFY_CHIPINFO_FILEPATH) -chipTypeQual=$(PIPELINEROOT)$(AFFY_CHIPTYPEQUAL_FILEPATH) \
       -affymetrixChip=$(PIPELINEROOT)$(AFFY_CHIP_FILEPATH) -processed_schuster=$(SCHUSTERPATH)   before
```
(Before normalization, to detect problems in the annotation and files. Be careful, the script should be able to connect to the database, and it should be run on the computer having all the data. It checks a lot of small common mistakes by the annotators)


* Check that there is no problem with the generation of the information.
```
perl Annotation_checking/check_affy_info.pl  -affyChipInformation=$(PIPELINEROOT)$(AFFY_CHIPINFO_FILEPATH) \
       -chipTypeQual=$(PIPELINEROOT)$(AFFY_CHIPTYPEQUAL_FILEPATH) -affymetrixChip=$(PIPELINEROOT)$(AFFY_CHIP_FILEPATH)
```
This will check for duplicated files, incompatible chip types, and will inform you about low quality chips. You need to take actions for duplicated chips (comment the chips you want to remove in [pipeline/Affymetrix/affymetrixChip](affymetrixChip)), and incompatible chip types (if the chip type is not already present in [pipeline/Affymetrix/chipTypeCorrespondencesAndQualityThresholds](chipTypeCorrespondencesAndQualityThresholds); if it is already present, you don't need to do anything). For low quality chips, don't bother.


### Potential problems
With processed data (processed_mas5):
* Annotation problems: the name in the annotation file does not correspond to the name of the files in the experiment folder.
* Corrupted files: sometimes the exported files on ArrayExpress (AE) do not include the probeIds. The first column if the downloaded files is empty -> find the list of probes from another experiment using the same chip, or contact AE so that they correct their file.
* Corrupted files: sometimes files don't have the correct number of lines (probes). This is mysterious and is probably due to a bad submission to AE.


## How to remove a chip?
Deprecate the chip/experiment from the annotation files (comment it?)

If the chip/experiment has already been used in the pipeline, you have to remove it from `/var/bgee/extra/pipeline/Affymetrix/cel_data/` and `/var/bgee/extra/pipeline/Affymetrix/processed_schuster/` or `processed_mas5/`. As well as from the different backup servers.

Better to keep the chip/experiment anyway, so move the chip/experiment cell file to `/var/bgee/extra/pipeline/Affymetrix/cel_data/retired/` keeping the chip/experiment structure.

Do the same for result file in `/var/bgee/extra/pipeline/Affymetrix/processed_schuster/retired/` or `processed_mas5/retired/`.


## Normalization
```
make normalization
```

### Raw data normalization (using gcRMA, when several ".cel" files are available in an experiment)
```
perl Normalization/launch_affy_analysis.pl  -affyChipInformation=$(PIPELINEROOT)$(AFFY_CHIPINFO_FILEPATH) \
       -chipTypeQual=$(PIPELINEROOT)$(AFFY_CHIPTYPEQUAL_FILEPATH) -affymetrixChip=$(PIPELINEROOT)$(AFFY_CHIP_FILEPATH) \
       -cel_data=$(CELPATH) -processed_schuster=$(SCHUSTERPATH) -bioconductorout=$(BIOCONDUCTOROUT) -bioconductoraffin=$(BIOCONDUCTORAFFIN)
(or perl Normalization/launch_affy_analysis.pl ... EXP_ID if you want to analyse only one experiment)
```

Check normalization process by looking if some "*.out.PROB" files exist in `$(BIOCONDUCTOROUT)` AND if all "*.out" files contains a line with `proc.time()` at the end:
```
ls -l ls -l $(BIOCONDUCTOROUT)*.out.PROB
grep -H -c 'proc.time' $(BIOCONDUCTOROUT)*.out | grep ':0'
```

Also check warnings:
```
grep -i "warning" -A 5 $(BIOCONDUCTOROUT)*.out
```


### If there is a problem (files ".out.PROB")
It might be that there is only one ".cel" file used from this experiment. In that case, it's not possible to use gcRMA. Check in the annotation file directly that this is indeed the problem (you can have several ".cel" files in the directory of the experiment, but still having only one chip annotated). In that case, see below. Otherwise, see the **potential problems** section. Examples of error messages:
```
> data.gcrma <- gcrma(data, affinity.info=ai, type="affinities")
Adjusting for optical effect.Done.
Error in model.frame(formula, rownames, variables, varnames, extras, extranames,  :
   variable lengths differ (found for 'x')
> data.gcrma <- gcrma(data, affinity.info=ai, type="affinities")
Adjusting for optical effect.Done.
Error in model.frame.default(formula = y ~ x, drop.unused.levels = TRUE) :
variable lengths differ (found for 'x')
Calls: gcrma ... lm -> eval -> eval -> model.frame -> model.frame.default
```

It is now possible to normalize using only one CEL file, this part should not be needed anymore. If it should become needed again in the future, one should check that the chip passed quality controls (percent present and compatible chip type) before launching the analysis.


### Note for large ".cel" files
If large cel files have to be normalized, you can use cluster computers with a lot of memory. Copy required files, data structure and scripts there:
```
rsync -Wav -essh --exclude '*.gz' $(CELPATH)       YYYYYY@dee-serv02.vital-it.ch:/scratch/temporary/YYYYYY/bgee/Affymetrix/cel_data/
rsync -Wav -essh                  $(BIOCONDUCTOR)  YYYYYY@dee-serv02.vital-it.ch:/scratch/temporary/YYYYYY/bgee/Affymetrix/bioconductor/
rsync -Wav -essh                  $(AFFYPATH)      YYYYYY@dee-serv02.vital-it.ch:/scratch/temporary/YYYYYY/bgee/Affymetrix/
ssh YYYYYY@dee-serv02.vital-it.ch
module add R/3.2.2
cd /scratch/temporary/YYYYYY/bgee/Affymetrix/
RUN normalization scripts
```
* `YYYYYY` being your login account on cluster, and the repository in `/scratch/temporary/` you need to ask for.
* `dee-serv03.vital-it.ch` could also be used !
* You may have to change some paths in order to reach `/scratch/temporary/YYYYYY/bgee/Affymetrix/` properly!


### Note if you need to re-run some normalizations
If you need to re-run normalization for some experiments or all experiments:
* remove the corresponding ".out" file(s) in `$(BIOCONDUCTOROUT)`. File names use the pattern **experimentId_arrayId.out**. You should remove every files corresponding to a given experiment, just in case a previously annotated chip was removed from current analyses.
* remove the corresponding file(s) in `$(SCHUSTERPATH)`. File names use the pattern **experimentId/arrayId**. You should remove the whole path corresponding to an experiment, just in case a previously annotated chip was removed from current analyses.
* uncompress the corresponding ".cel" files in `$(CELPATH)` (do not forget to re-compress them at the end).


### Potential problems
* If the annotation package for the chipType is missing check the list of all annotations packages [here](http://www.bioconductor.org/packages/2.2/ChipName.html). You may need to install it on the devbioinfo server. If it's on cluster, ask them to install it [mailto:projects@vital-it.ch](mailto:projects@vital-it.ch).
* Annotations errors. A problem happens if annotators put the wrong chip type. This is often the case when multiple chip types are used in the same experiment. Sometimes the mistake is present in ArrayExpress! -> You can have a look at the chip type directly in the header of the CEL file.
* Corrupted files. Usually there is nothing to do... You can remove the problematic chip from the experiment in the annotation file.
* Custom chips: the annotation package does not exist in bioconductor (encode chips for example)
* Memory problems: shouldn't occur on cluster (256Go memory)
* Other problems may be solved using the [bioconductor mailing list](http://dir.gmane.org/gmane.science.biology.informatics.conductor) usually. If the problem is strange, you may want to avoid losing too much time and be pragmatic: throw away the experiment or chip (remove it from `affymetrixChip` and `microarrayExperiment`).


## After normalization
```
make check_after
```
* In [pipeline/Affymetrix/](.) check that there is no problem with the normalized files:
```
perl Annotation_checking/check_affy_curation.pl  -bgee=$(BGEECMD) -normalizationType=$(PIPELINEROOT)$(AFFY_NORMTYPE_FILEPATH) \
       -detectionType=$(PIPELINEROOT)$(AFFY_DETCTYPE_FILEPATH) -chipType=$(PIPELINEROOT)$(AFFY_CHIPTYPE_FILEPATH) \
       -microarrayExperiment=$(PIPELINEROOT)$(MICROARRAY_EXPERIMENT_FILEPATH) -cel_data=$(CELPATH) -processed_mas5=$(MAS5PATH) \
       -affyChipInformation=$(PIPELINEROOT)$(AFFY_CHIPINFO_FILEPATH) -chipTypeQual=$(PIPELINEROOT)$(AFFY_CHIPTYPEQUAL_FILEPATH) \
       -affymetrixChip=$(PIPELINEROOT)$(AFFY_CHIP_FILEPATH) -processed_schuster=$(SCHUSTERPATH)   after
```

* You have then to compress the cel files not yet compressed (new files).
```
find $(CELPATH) ! -name \*.gz -type f -exec gzip -9 {} \;
```

* **Think to back-up these processed cel files at the end.** Backup the whole `$(AFFYDATAPATH)` repository on another machine for safety !
```
rsync -av --del -f "- /lost+found" --exclude=.svn --exclude=.git $(AFFYDATAPATH)  $(DATALOGIN)@$(DATAHOST):$(DATAPATH)Affymetrix/
```


## Data insertion


### Insert the data
* Before insertion, list all unique strains present in our annotations, to merge inconsistent strains, e.g., 'Wild Type' or 'wild type' instead of 'wild-type' (see, e.g., https://gitlab.isb-sib.ch/Bgee/bgee_pipeline/issues/67#note_4654)
* Useful commands to have a list of major strains potentially redundant

#### C57BL/6 129... not mixed not and
* C57BL/6 129/Sv
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep -i "sv" | grep -iv "ev" | grep -iv "and" | grep -iv "mix" | grep -iv "6J"
```
* C57BL/6 129/SvEv
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep -i "sv" | grep -i "ev" | grep -iv "and" | grep -iv "mix" | grep -iv "6J"
```
* C57BL/6 129
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep "129" | grep -vi "sv" | grep -vi "ev" | grep -iv "and" | grep -iv "mix" | grep -iv "6J"
```
* C57BL/6J 129/Sv
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep -i "sv" | grep -iv "ev" | grep -iv "and" | grep -iv "mix" | grep -i "6J"
```
* C57BL/6J 129/SvEv
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep -i "sv" | grep -i "ev" | grep -iv "and" | grep -iv "mix" | grep -i "6J"
```
* C57BL/6 129
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep "129" | grep -vi "sv" | grep -vi "ev" | grep -iv "and" | grep -iv "mix" | grep -i "6J"
```


#### C57BL/6 129... not mixed, and
* C57BL/6 129/Sv
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep -i "sv" | grep -iv "ev" | grep -i "and" | grep -iv "mix" | grep -iv "6J"
```

* C57BL/6 129/SvEv
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep -i "sv" | grep -i "ev" | grep -i "and" | grep -iv "mix" | grep -iv "6J"
```

* C57BL/6 129
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep "129" | grep -vi "sv" | grep -vi "ev" | grep -i "and" | grep -iv "mix" | grep -iv "6J"
```

* C57BL/6J 129/Sv
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep -i "sv" | grep -iv "ev" | grep -i "and" | grep -iv "mix" | grep -i "6J"
```

* C57BL/6J 129/SvEv
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep -i "sv" | grep -i "ev" | grep -i "and" | grep -iv "mix" | grep -i "6J"
```

* C57BL/6 129
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep "129" | grep -vi "sv" | grep -vi "ev" | grep -i "and" | grep -iv "mix" | grep -i "6J"
```

#### C57BL/6 129... mixed, not and
* C57BL/6 129/Sv
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep -i "sv" | grep -iv "ev" | grep -iv "and" | grep -i "mix" | grep -iv "6J"
```

* C57BL/6 129/SvEv
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep -i "sv" | grep -i "ev" | grep -iv "and" | grep -i "mix" | grep -iv "6J"
```

* C57BL/6 129
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep "129" | grep -vi "sv" | grep -vi "ev" | grep -iv "and" | grep -i "mix" | grep -iv "6J"
```

* C57BL/6J 129/Sv
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep -i "sv" | grep -iv "ev" | grep -iv "and" | grep -i "mix" | grep -i "6J"
```

* C57BL/6J 129/SvEv
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep -i "sv" | grep -i "ev" | grep -iv "and" | grep -i "mix" | grep -i "6J"
```

* C57BL/6 129
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep "129" | grep -vi "sv" | grep -vi "ev" | grep -iv "and" | grep -i "mix" | grep -i "6J"
```

#### C57BL/6 129... mixed, and
* C57BL/6 129/Sv
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep -i "sv" | grep -iv "ev" | grep -i "and" | grep -i "mix" | grep -iv "6J"
```

* C57BL/6 129/SvEv
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep -i "sv" | grep -i "ev" | grep -i "and" | grep -i "mix" | grep -iv "6J"
```

* C57BL/6 129
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep "129" | grep -vi "sv" | grep -vi "ev" | grep -i "and" | grep -i "mix" | grep -iv "6J"
```

* C57BL/6J 129/Sv
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep -i "sv" | grep -iv "ev" | grep -i "and" | grep -i "mix" | grep -i "6J"
```

* C57BL/6J 129/SvEv
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep -i "sv" | grep -i "ev" | grep -i "and" | grep -i "mix" | grep -i "6J"
```

* C57BL/6 129
```
cut -f23 source_files/Affymetrix/affymetrixChip.tsv | sort | uniq | grep -i "C57BL" | grep "129" | grep -vi "sv" | grep -vi "ev" | grep -i "and" | grep -i "mix" | grep -i "6J"
```

* `make insert_affy`

  * Fill the rnaSeqLibrary, affymetrixProbeset, ... tables.
  ```
  perl Data_insertion/insert_affy.pl -bgee=$(BGEECMD) -ensembl=$(ENSCMD)  -chipType=$(PIPELINEROOT)$(AFFY_CHIPTYPE_FILEPATH) \
     -chipTypeQual=$(PIPELINEROOT)$(AFFY_CHIPTYPEQUAL_FILEPATH) -microarrayExperiment=$(PIPELINEROOT)$(MICROARRAY_EXPERIMENT_FILEPATH) \
     -affymetrixChipInfo=$(PIPELINEROOT)$(AFFY_CHIPINFO_FILEPATH) -affymetrixChip=$(PIPELINEROOT)$(AFFY_CHIP_FILEPATH) \
     -annotations=$(ANNOTATIONPATH) -processed_mas5=$(MAS5PATH) -processed_schuster=$(SCHUSTERPATH)  -exp=both
  ```

  * See remarks and warnings to see details about the error and warning messages!
  * Note: If you need to rerun the script because something went wrong, delete information from **expression** and **noExpression** tables without removing information inserted by EST or RNASeq steps: `make deleteAffy`

* Delete `microarrayExperiment` with no `affymetrixChip`
```
DELETE t1 FROM microarrayExperiment AS t1 LEFT OUTER JOIN affymetrixChip AS t2 ON t1.microarrayExperimentId = t2.microarrayExperimentId
WHERE t2.microarrayExperimentId IS NULL
```

* Check carefully the statistics and strains inserted in file `insert_affy`.

### Remarks/Warnings
* [Data_insertion/insert_affy.pl](Data_insertion/insert_affy.pl) messages:
  * _Warning! Incorrectly formated file: [...]_

  It means that, in a `processed_mas5 file`, the script could not find the column for mas5 calls and/or the column for signal intensities. Look at the files to see if this is an error on our side. Most likely, it is a problem with ArrayExpress that you should report to them.

  * _Warning, unrecognized header..._

  is harmless as long as you don't have another error message following it, like:
  ```
  Warning: unrecognized line:
  ```

  * _Warning! The mapping in the ".out" file is outdated (XXX probesets are no longer mapped)_

  It means custom mapping for some arrays (Su et al.) for which some gene id have changed in Ensembl. So it's OK and harmless.

  * BUT, if you have a message such as: **Warning! There seems to be a bad mapping for this chip: 0 probesets**

  It means you encountered a real issue. Check carefully the chipType file that you generated at a previous step. In the 4th column, you supposed to have, either the name of the table to retrieve the mapping in the Ensembl database, or the name of the mapping file we generated, stored in `$(ANNOTATIONPATH)`. If it is a Ensembl table name, then maybe this table has disappeared, check in the Ensembl database. If it is a mapping file, then either this mapping file is absent, or obsolete... Check in the wiki of the annotators to get more details about how these mapping files are generated. Or maybe it is just the wrong mapping file that is indicated, etc.

  It means that in a mas5 file or a processed schuster file, instead of getting a P (Present), M (Marginal), or A (Absent), you get an unknown signal. There is something wrong with the file and you have to check it, and maybe add this signal in the list of known signals.
  ```
  Warning! Some expression calls are not standard (X) for ...
  ```

  It means there are more than 2% of "null" signal in the file. Maybe this file should be discarded.
  ```
  Warning! Some expression calls have too many null values (X%) for ...
  ```

  Such warnings can then be followed by
  ```
  Warning! Problem with gene ENSMUSG00000005232: not possible to infer summary of presence and quality
  ```

  It means that, when trying to defined the expression state and the quality for this gene, we got no data (presence high qual = 0; presence low qual = 0; absence high qual = 0; absence low qual = 0) If you have this warning without the previous one for the same chip, then you get a problem...


* The script [Annotation_checking/check_affy_curation.pl](Annotation_checking/check_affy_curation.pl) could be modified to signal when an experiment is both in the annotation file and the "not_included" or "not_included_for_now" files.
* Consensus concerning probeset quality (when multiple probesets are present for the same genes):

| Pst/High | Pst/Low | Abs/High | Abs/Low | Consensus                 |
|:--------:|:-------:|:--------:|:-------:|---------------------------|
| 1        | 1       | 0        | 0       | Pst/High                  |
| 1        | 0       | 0        | 0       | Pst/High                  |
| 1        | 1       | 1        | 1       | Pst/Low                   |
| 1        | 1       | 1        | 0       | Pst/Low                   |
| 1        | 1       | 0        | 1       | Pst/Low                   |
| 1        | 0       | 1        | 1       | Pst/Low                   |
| 1        | 0       | 1        | 0       | Pst/Low                   |
| 1        | 0       | 0        | 1       | Pst/Low                   |
| 0        | 1       | 1        | 1       | Pst/Bronze (not inserted) |
| 0        | 1       | 1        | 0       | Pst/Bronze (not inserted) |
| 0        | 1       | 0        | 1       | Pst/Bronze (not inserted) |
| 0        | 1       | 0        | 0       | Pst/Low                   |
| 0        | 0       | 1        | 0       | Abs/High                  |
| 0        | 0       | 1        | 1       | Abs/High                  |
| 0        | 0       | 0        | 1       | Abs/Low (not inserted)    |
| 0        | 0       | 0        | 0       | Not possible              |

* Note that in some cases the probesets are inserted in `affymetrixProbeset` table, but nothing is inserted into the `expression`/`noExpression` tables. This is the case for:
  * Consensus Abs/Low (Abs seen only by mas5)
  * Consensus Pst/Bronze (no probeset Pst/High is seen for this gene/stage/organ)
  * Pre-filtering: probesets that are never Pst/High on the whole dataset with gcRMA, and never Pst(/Low) with mas5. Detection flags are always absent or marginal.

* Note: the same consensus table is used for _in situ_ data


## Expression data insertion

### Insert the data
```
make insert_expression
```

* Fill the `expression` and `noExpression` tables and update `affymetrixProbeset` table.
```
perl Data_insertion/insert_affy_expression.pl -bgee=$(BGEECMD) -type=both
```

* Check carefully the statistics and strains inserted in file `insert_expression`.


## Explanations about the Affymetrix files and directories
* `$(AFFYPATH)`
  * `$(AFFYPATH)affymetrixChip.tsv`: the file from annotators
  * `$(AFFYPATH)affymetrixChipInformation`: for each chip, the cdf name, quality score, percent present, unique_id, cel data checksum, processed mas5 original file checksum, processed mas5 filtered file checksum
  * `$(AFFYPATH)chipTypeCorrespondencesAndQualityThresholds`: for each cdf name, correspondence to chipTypeId, status (compatible/incompatible), quality score threshold, percent present threshold
  * `$(AFFYPATH)Generate_information/`: R scripts used to generate the info. These scripts are launched by `$(AFFYPATH)Generate_information/generate_affy_chips_information.pl`
* `$(AFFYDATAPATH)chip_information/`: folder to store scripts and results of generating information on affy chips
  * `$(AFFYDATAPATH)chip_information/logs/`: R output of the R scripts, with the name 'experimentId_chipId.out'. Renamed with ".out.PROB" when a problem occurs.
  * `$(AFFYDATAPATH)chip_information/results/`: results of the analyses of the R scripts, with names experimentId_chipId.tsv. Data are extracted from these files using `$(AFFYPATH)Generate_information/generate_affy_chips_information.pl`
* `$(AFFYDATAPATH)bioconductor/`: R scripts used for affy analyses.
  * `$(AFFYDATAPATH)bioconductor/affinities/`: R data generated by `$(AFFYPATH)Normalization/affy_analysis.R`, launched by `$(AFFYPATH)Normalization/launch_affy_analysis.pl`. Names follow the pattern 'chipTypeId.RData'
  * `$(AFFYDATAPATH)bioconductor/out/`: outputs of `$(AFFYPATH)Normalization/affy_analysis.R`, launched by `$(AFFYPATH)Normalization/launch_affy_analysis.pl`. Names follow the pattern 'experimentId_chipTypeId.out'. These files are read by `launch_affy_analysis.pl` to determine whether a normalization has already been performed for this chipType in this experiment. Should be removed to re-normalized.
  * `$(AFFYDATAPATH)bioconductor/differential/`: outputs of `extra/pipeline/Affymetrix/bioconductor/diff_analysis.R`, launched by `extra/pipeline/pipeline/launch_diff_analysis-Affy.pl`. Names follow the pattern 'experimentId_chipTypeId.out'. These files are used by `launch_diff_analysis-Affy.pl` to determine whether a diff analysis has already been performed for this chipType in this experiment. Should be removed to re-launch diff analysis.
  * `$(AFFYDATAPATH)bioconductor/targets/`: ".target" files, created by `extra/pipeline/pipeline/launch_diff_analysis-Affy.pl`, used by `extra/pipeline/Affymetrix/bioconductor/diff_analysis.R` to retrieve which files to used containing signal intensities of probesets (either in `$(AFFYDATAPATH)processed_schuster/`, or `$(AFFYDATAPATH)processed_mas5/`), and by `extra/pipeline/pipeline/insert_diff_affy.pl` to retrieve the chipIds used in each analysis (either by assuming that files names follow the pattern 'chipId.cel.out' case insensitive or 'chipId.out', or the pattern 'chipId'). They contain information about the organs/stages to analyze, and the raw data files names (i.e., 'chipId.cel.out'), e.g., GSM341502.CEL.out `../Affymetrix/processed_schuster/`. Names of these files follow the pattern 'experimentId___chipTypeId.target'. They should be removed if annotations or chips included have changed for an experiment. If a chipId has changed, the ".target" file should be updated.
* `$(AFFYDATAPATH)annotations/`: files used to create custom mapping from probesets to Ensembl genes, usually when it is a custom chip with no mappings found in Ensembl. Just kept for the records basically.
* `$(AFFYDATAPATH)cel_data/`: CEL files, stored experiment by experiment (e.g., `cel_data/expId1/`, `cel_data/expId2/`). Note that files that are gunziped have already been normalized. Files not gunziped are supposed to be normalized during the next pipeline run. Files recently annotated to be added in the next pipeline run are first stored on an annotator computer (currently Anne's computer, Anne's account), then moved to `devbioinfo` for normalization and analysis, and removed from the annotator computer. All files are also backuped on an annotator computer (currently Anne's computer, bgee's account). Names of the CEL files follow the pattern 'chipId.cel' or 'chipId.CEL' (potentially ".gz" if already normalized). You should gunzip all CEL files of an experiment you want to re-normalize.
* `$(AFFYDATAPATH)processed_differential/`: results of diff analyses, stored experiment by experiment (e.g., `processed_differential/expId1/`, `processed_differential/expId2/`). These directories are created by `extra/pipeline/pipeline/launch_diff_analysis-Affy.pl` before launching `extra/pipeline/Affymetrix/bioconductor/diff_analysis.R`, that writes the ".out" files in the created directory, containing results of the analysis. Names of the files follow the pattern 'chipTypeId_factors.out', i.e. 'chipTypeId_organId_stageId.out'. They contains for each probeset the foldchange and the associated p-value from the multiple comparisons to the mean procedure. The directory of an experiment should be removed if you want to re-analyse it. These files are read by `extra/pipeline/pipeline/insert_diff_affy.pl` to insert results into the Bgee database.
* `$(AFFYDATAPATH)processed_schuster/`: results of the Schuster analyses, generated by `$(AFFYPATH)Normalization/affy_analysis.R`. Files stored experiment by experiment (script launched and directories created by `$(AFFYPATH)Normalization/launch_affy_analysis.pl`). Names of the files follow the pattern 'CELFile.out', i.e. 'chipId.cel.out' or 'chipId.CEL.out'. These files contain normalized signal intensities and detection flag for each probeset, based on the Schuster procedure. The content of these files are used by `$(AFFYPATH)Data_insertion/insert_affy.pl` to generate files inserted into the Bgee database, and by `extra/pipeline/Affymetrix/bioconductor/diff_analysis.R` to perform the differential analyses. The files names are used in `extra/pipeline/pipeline/launch_diff_analysis-Affy.pl` to generate the ".target" files in `$(AFFYDATAPATH)bioconductor/targets/`. These ".target" files will be used by `extra/pipeline/Affymetrix/bioconductor/diff_analysis.R` to identify which files to used to perform the differential analyses, and by `extra/pipeline/pipeline/insert_diff_affy.pl` to retrieve the chips used in each differential analysis. If a chipId has changed, name of the ".out" file should be modified.
* `$(AFFYDATAPATH)processed_mas5_original_files/`: results of the MAS5 analyses, generated by `extra/pipeline/Affymetrix/bioconductor/affy_analysis_mas5.R` when we have only one CEL file for an experiment (not possible to normalize using gcRMA), or MAS5 outputs directly downloaded from ArrayExpress or GEO repositories when no CEL files are available. Files stored experiment by experiment (directories created manually by curators). Names of the files follow the pattern 'CELFile.out', or just 'fileName', i.e. 'chipId.cel.out' or 'chipId.CEL.out', or just 'chipId'. These files contain normalized signal intensities and detection flag for each probeset, based on the MAS5 procedure. These files are used by "clean_mas5_files.pl" to generate "clean" files, stored in `$(AFFYDATAPATH)processed_mas5/` (see below).
* `$(AFFYDATAPATH)processed_mas5/`: same as `$(AFFYDATAPATH)processed_mas5_original_files/`, but the files are converted to a proper format: sometimes the columns of the mas5 files are not in the proper order, or there are additional header lines, etc... All the files are converted to a proper format, where columns are always in the same order, and where there is only one header line. This is achieved by the script `$(AFFYPATH)MAS5/clean_mas5_files.pl`. All cleaned files are stored in this directory, and all original files are kept in `$(AFFYDATAPATH)processed_mas5_original_files/`. All scripts then used the cleaned file in `processed_mas5/`: the content of these files are used by `$(AFFYPATH)Data_insertion/insert_affy.pl` to generate files inserted into the Bgee database, and by `extra/pipeline/Affymetrix/bioconductor/diff_analysis.R` to perform the differential analyses. The files names are used in `extra/pipeline/pipeline/launch_diff_analysis-Affy.pl` to generate the ".target" files in `$(AFFYDATAPATH)bioconductor/targets/`. These ".target" files will be used by `extra/pipeline/Affymetrix/bioconductor/diff_analysis.R` to identify which files to used to perform the differential analyses, and by `extra/pipeline/pipeline/insert_diff_affy.pl` to retrieve the chips used in each differential analysis. If a chipId has changed, name of the ".out" file should be modified.

