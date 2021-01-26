PIPELINEROOT := ../
DIR_NAME := Affymetrix/
include $(PIPELINEROOT)Makefile.common


all: $(VERIFICATIONFILE)


#######################
# Maintenance
#NOTE  /!\  need to put makefile runner ssh key in ANNOTATORHOST $(ANNOTATORLOGIN)_HOME/.ssh/authorized_keys end of file
# see http://troy.jdmz.net/rsync/
# OR to use a anonymous/guest login without password on ANNOTATORHOST
new_annotation:
	# Download new 'cel' annotated data
#TODO re-do rsync -a to fix permission problems directly!
	@$(RSYNC) $(ANNOTATORLOGIN)@$(ANNOTATORHOST):$(CELANNPATH)/      $(CELNEWPATH)/     >  $@.tmp
	# Download new 'MAS5 only' annotated data
	@$(RSYNC) $(ANNOTATORLOGIN)@$(ANNOTATORHOST):$(MAS5ANNPATH)/     $(MAS5NEWPATH)/    >> $@.tmp
	@$(RSYNC) $(ANNOTATORLOGIN)@$(ANNOTATORHOST):$(MAS5ORIANNPATH)/  $(MAS5ORINEWPATH)/ >> $@.tmp
# Clean mac useless files
	@find $(AFFYNEWDATAPATH)/ -type f -name .DS_Store    | sudo xargs $(RM)
	@find $(AFFYDATAPATH)/    -type f -name .DS_Store    | sudo xargs $(RM)
	@find $(AFFYDATAPATH)/    -type f -name .DS_Store.gz | sudo xargs $(RM)
# Reminder
	@echo -e "\t\e[1m\e[37m\e[41mRemove them from the annotators computers to not re-process them next time !\e[m\n"
	@$(MV) $@.tmp $@

maintenance: new_annotation
	# Check file presence, without moving them
	@perl Maintenance/check_and_move_new_files.pl $(AFFYNEWDATAPATH)/  $(AFFYDATAPATH)/
	# Move new Affy annotated files, if any
	@perl Maintenance/check_and_move_new_files.pl $(AFFYNEWDATAPATH)/  $(AFFYDATAPATH)/  1 > $@.tmp
	@$(MV) $@.tmp $@

get_annot: maintenance
	# Retrieve up-to-date annotation files from https://gitlab.sib.swiss/Bgee/expression-annotations
	@$(WGET) $(ANNOTATION_GIT_URL)/Strains/StrainMapping.tsv                && $(MV) StrainMapping.tsv             $(STRAIN_MAPPING_FILE)
	@$(WGET) $(ANNOTATION_GIT_URL)/Affymetrix/affymetrixChip.tsv            && $(MV) affymetrixChip.tsv            $(AFFY_CHIP_FILEPATH)                  2>/dev/null  || rm -f affymetrixChip.tsv
	@$(WGET) $(ANNOTATION_GIT_URL)/Affymetrix/chipType.tsv                  && $(MV) chipType.tsv                  $(AFFY_CHIPTYPE_FILEPATH)              2>/dev/null  || rm -f chipType.tsv
	@$(WGET) $(ANNOTATION_GIT_URL)/Affymetrix/microarrayExperiment.tsv      && $(MV) microarrayExperiment.tsv      $(MICROARRAY_EXPERIMENT_FILEPATH)      2>/dev/null  || rm -f microarrayExperiment.tsv
# Download wormbase annotation as well. Will be merged with our annotation
	@$(WGET) $(ANNOTATION_GIT_URL)/Affymetrix/affymetrixChip_worm.tsv       && $(MV) affymetrixChip_worm.tsv       $(AFFY_CHIP_FILEPATH_WORM)             2>/dev/null  || rm -f affymetrixChip_worm.tsv
	@$(WGET) $(ANNOTATION_GIT_URL)/Affymetrix/chipType_worm.tsv             && $(MV) chipType_worm.tsv             $(AFFY_CHIPTYPE_FILEPATH_WORM)         2>/dev/null  || rm -f chipType_worm.tsv
	@$(WGET) $(ANNOTATION_GIT_URL)/Affymetrix/microarrayExperiment_worm.tsv && $(MV) microarrayExperiment_worm.tsv $(MICROARRAY_EXPERIMENT_FILEPATH_WORM) 2>/dev/null  || rm -f microarrayExperiment_worm.tsv
	@touch $@

prepare_annot: get_annot
	# Concatenate our annotation file with wormbase annotations (without header, and sorted by experiments and libraries)
	@tail -n+2 $(AFFY_CHIP_FILEPATH_WORM)             | sort -k2,2 -k1,1 | cat $(AFFY_CHIP_FILEPATH)             - > $(AFFY_CHIP_FILEPATH_FULL).ori
	@tail -n+2 $(AFFY_CHIPTYPE_FILEPATH_WORM)         | sort -k1,1       | cat $(AFFY_CHIPTYPE_FILEPATH)         - > $(AFFY_CHIPTYPE_FILEPATH_FULL)
	@tail -n+2 $(MICROARRAY_EXPERIMENT_FILEPATH_WORM) | sort -k1,1       | cat $(MICROARRAY_EXPERIMENT_FILEPATH) - > $(MICROARRAY_EXPERIMENT_FILEPATH_FULL)
	# Map strain names
	@perl -e 'use lib ".."; use Utils; Utils::map_strain_names("$(AFFY_CHIP_FILEPATH_FULL).ori", "$(STRAIN_MAPPING_FILE)")' >$(AFFY_CHIP_FILEPATH_FULL)
	@$(RM) $(AFFY_CHIP_FILEPATH_FULL).ori
	@touch $@

#FIXME Need to properly set directory rights: Some are owned by bgee users, others by yourself, others by root (e.g. R library folder)


#######################
# MAS5
clean_mas5: prepare_annot
	# Clean up MAS5 files
	@perl MAS5/clean_mas5_files.pl $(MAS5ORIPATH) $(MAS5PATH) > $@.out 2> $@.err
# remove $@.err if empty
	@find . -type f -empty -name $@.err -exec $(RM) {} \;
	@find . -type f        -name $@.err -exec cat   {} \;
	@$(MV) $@.out $@

check_mas5: clean_mas5 $(AFFY_CHIP_FILEPATH_FULL)
	# Global diff
	-@diff -rq $(MAS5ORIPATH)/ $(MAS5PATH)/ | grep -v ' diff' | grep -v '\.svn' | grep -v '\.git' | grep -v 'not_separated' | cat > $@.tmp
	# Check if invalid files are present
	@perl MAS5/check_mas5_filtered.pl  -affyChipFilesDir=$(AFFYDATAPATH)  -affymetrixChip=$(AFFY_CHIP_FILEPATH_FULL) >> $@.tmp
	@$(MV) $@.tmp $@

adj_call: check_mas5 $(AFFY_CHIP_FILEPATH) $(MAS5PATH)
	# Correction of the calls
	@echo --- starting the calls correction for all MAS5 files ---
	@$(CLUSTER_R_CMD) R CMD BATCH --no-save --no-restore '--args affy_anno="$(AFFY_CHIP_FILEPATH)" mas5_path="$(MAS5PATH)" cut_off="0.01"' MAS5/adj_calls_mas5.R MAS5/adj_calls_mas5.Rout > $@.tmp
	@$(MV) MAS5/adj_calls_mas5.Rout $@.log
	@echo --- DONE ---
	@$(MV) $@.tmp $@

#######################
# Generate information
gen_info: prepare_annot $(AFFY_CHIP_FILEPATH_FULL) $(AFFY_CHIPINFO_FILEPATH) $(AFFY_CHIPTYPEQUAL_FILEPATH) adj_call
	@which R >/dev/null
	# Generate information
	@perl Generate_information/generate_affy_chips_information.pl  -affyChipFilesDir=$(AFFYDATAPATH)  -affymetrixChip=$(AFFY_CHIP_FILEPATH_FULL) -affymetrixChipInformation=$(AFFY_CHIPINFO_FILEPATH) -chipTypeQual=$(AFFY_CHIPTYPEQUAL_FILEPATH) > $@.tmp
	# Checks
	-@ls -l $(AFFYDATAPATH)chip_information/logs/*.out.PROB                          >> $@.tmp
	-@grep -H -c 'proc.time' $(AFFYDATAPATH)chip_information/logs/*.out | grep ':0'  >> $@.tmp
	-@grep -i "warning" -A 5 $(AFFYDATAPATH)chip_information/logs/*.out              >> $@.tmp
	@$(GIT) diff $(AFFY_CHIPINFO_FILEPATH)                                           >> $@.tmp
	@$(GIT) diff $(AFFY_CHIPTYPEQUAL_FILEPATH)                                       >> $@.tmp
	@$(MV) $@.tmp $@


#######################
# Annotations checking
check_annot: gen_info $(AFFY_NORMTYPE_FILEPATH) $(AFFY_DETCTYPE_FILEPATH) $(AFFY_CHIPTYPE_FILEPATH_FULL) $(MICROARRAY_EXPERIMENT_FILEPATH_FULL) $(AFFY_CHIPINFO_FILEPATH) $(AFFY_CHIPTYPEQUAL_FILEPATH) $(AFFY_CHIP_FILEPATH_FULL)
	# Check curation before
	@perl Annotation_checking/check_affy_curation.pl  -bgee=$(BGEECMD) -normalizationType=$(AFFY_NORMTYPE_FILEPATH) -detectionType=$(AFFY_DETCTYPE_FILEPATH) -chipType=$(AFFY_CHIPTYPE_FILEPATH_FULL) -microarrayExperiment=$(MICROARRAY_EXPERIMENT_FILEPATH_FULL) -cel_data=$(CELPATH) -processed_mas5=$(MAS5PATH) -affyChipInformation=$(AFFY_CHIPINFO_FILEPATH) -chipTypeQual=$(AFFY_CHIPTYPEQUAL_FILEPATH) -affymetrixChip=$(AFFY_CHIP_FILEPATH_FULL) -processed_schuster=$(SCHUSTERPATH)   before                > $@.tmp
	@echo >> $@.tmp
	# Check info
	@perl Annotation_checking/check_affy_info.pl  -affyChipInformation=$(AFFY_CHIPINFO_FILEPATH) -chipTypeQual=$(AFFY_CHIPTYPEQUAL_FILEPATH) -affymetrixChip=$(AFFY_CHIP_FILEPATH_FULL) >> $@.tmp
	@$(MV) $@.tmp $@


#######################
# Normalization
normalization: check_annot $(AFFY_CHIPINFO_FILEPATH) $(AFFY_CHIPTYPEQUAL_FILEPATH) $(AFFY_CHIP_FILEPATH_FULL)
	@which R >/dev/null
	# Normalization
	@perl Normalization/launch_affy_analysis.pl  -affyChipInformation=$(AFFY_CHIPINFO_FILEPATH) -chipTypeQual=$(AFFY_CHIPTYPEQUAL_FILEPATH) -affymetrixChip=$(AFFY_CHIP_FILEPATH_FULL) -cel_data=$(CELPATH) -processed_schuster=$(SCHUSTERPATH) -bioconductorout=$(BIOCONDUCTOROUT) -bioconductoraffin=$(BIOCONDUCTORAFFIN) >$@.tmp
	# Checks
	-@find $(BIOCONDUCTOROUT) -empty -name \*.out.PROB | xargs $(RM)
	-@ls -l $(BIOCONDUCTOROUT)*.out.PROB                          >> $@.tmp
	-@grep -H -c 'proc.time' $(BIOCONDUCTOROUT)*.out | grep ':0'  >> $@.tmp
	-@grep -i "warning" -A 5 $(BIOCONDUCTOROUT)*.out              >> $@.tmp
	@$(MV) $@.tmp $@


#######################
check_after: normalization $(AFFY_NORMTYPE_FILEPATH) $(AFFY_DETCTYPE_FILEPATH) $(AFFY_CHIPTYPE_FILEPATH_FULL) $(MICROARRAY_EXPERIMENT_FILEPATH_FULL) $(AFFY_CHIPINFO_FILEPATH) $(AFFY_CHIPTYPEQUAL_FILEPATH) $(AFFY_CHIP_FILEPATH_FULL)
	# Check curation after
	@perl Annotation_checking/check_affy_curation.pl  -bgee=$(BGEECMD) -normalizationType=$(AFFY_NORMTYPE_FILEPATH) -detectionType=$(AFFY_DETCTYPE_FILEPATH) -chipType=$(AFFY_CHIPTYPE_FILEPATH_FULL) -microarrayExperiment=$(MICROARRAY_EXPERIMENT_FILEPATH_FULL) -cel_data=$(CELPATH) -processed_mas5=$(MAS5PATH) -affyChipInformation=$(AFFY_CHIPINFO_FILEPATH) -chipTypeQual=$(AFFY_CHIPTYPEQUAL_FILEPATH) -affymetrixChip=$(AFFY_CHIP_FILEPATH_FULL) -processed_schuster=$(SCHUSTERPATH)   after       > $@.tmp 2> warnings.$@
	@echo >> $@.tmp
# Gzip new cel files normalized, if any
	@find $(CELPATH) ! -name \*.gz -type f -exec gzip -9 {} \; >> $@.tmp
	# Rsync files for backup on sensitive cluster (see wiki for password)
	@echo -e "\t\e[1m\e[37m\e[41mRsync files in '$(AFFYDATAPATH)' to sensitive machine!\e[m\n"
	@rsync -av --del -f "- /lost+found" --exclude=.svn --exclude=.git $(AFFYDATAPATH)  $(DATALOGIN)@$(DATAHOST):$(DATAPATH)/Affymetrix/  >>  $@.tmp
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@


#######################
insert_affy: check_after $(AFFY_CHIPTYPE_FILEPATH_FULL) $(AFFY_CHIPTYPEQUAL_FILEPATH) $(MICROARRAY_EXPERIMENT_FILEPATH_FULL) $(AFFY_CHIPINFO_FILEPATH) $(AFFY_CHIP_FILEPATH_FULL) $(UBERON_SEX_INFO_FILE_PATH)
	# Start anatomy id mapping socket
	@$(IDMAPPING)  $(IDMAPPINGPORT)  &
	@$(STGMAPPING) $(STGMAPPINGPORT) &
	@sleep 50
	# Inserting affymetrix data. Debug mode is more verbose
	@perl Data_insertion/insert_affy.pl -bgee=$(BGEECMD) -ensembl=$(ENSCMD)  -chipType=$(AFFY_CHIPTYPE_FILEPATH_FULL) -chipTypeQual=$(AFFY_CHIPTYPEQUAL_FILEPATH) -microarrayExperiment=$(MICROARRAY_EXPERIMENT_FILEPATH_FULL) -affymetrixChipInfo=$(AFFY_CHIPINFO_FILEPATH) -affymetrixChip=$(AFFY_CHIP_FILEPATH_FULL) -annotations=$(ANNOTATIONPATH) -processed_mas5=$(MAS5PATH) -processed_schuster=$(SCHUSTERPATH) -sex_info=$(UBERON_SEX_INFO_FILE_PATH) -Aport=$(IDMAPPINGPORT) -Sport=$(STGMAPPINGPORT) -exp=both -debug > $@.tmp 2>warnings.$@
	# Delete microarrayExperiment with no affymetrixChip
	@$(MYSQL) -e "DELETE t1 FROM microarrayExperiment AS t1 LEFT OUTER JOIN affymetrixChip AS t2 ON t1.microarrayExperimentId = t2.microarrayExperimentId WHERE t2.microarrayExperimentId IS NULL" >> $@.tmp 2>> warnings.$@
	@echo "Check inconsistencies between condition species and gene species (there should be none): " >> $@.tmp
	@$(MYSQL) -e "SELECT t1.* FROM affymetrixProbeset AS t1 INNER JOIN gene AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId INNER JOIN affymetrixChip AS t3 ON t1.bgeeAffymetrixChipId = t3.bgeeAffymetrixChipId INNER JOIN cond AS t4 on t3.conditionId = t4.conditionId WHERE t2.speciesId != t4.speciesId" >> $@.tmp
	@echo >> $@.tmp
	@echo "Distinct strains in Affymetrix conditions, check that they are correct (e.g., no 'wild type' instead of 'wild-type')" >> $@.tmp
	@$(MYSQL) -e "SELECT DISTINCT t1.strain FROM cond AS t1 INNER JOIN affymetrixChip AS t2 ON t1.conditionId = t2.conditionId ORDER BY t1.strain" >> $@.tmp
	@echo >> $@.tmp
	@echo "Statistics on chips" >> $@.tmp
	@$(MYSQL) -e "SELECT t1.speciesId, COUNT(DISTINCT t2.bgeeAffymetrixChipId), MIN(t2.percentPresent), MAX(t2.percentPresent), AVG(t2.percentPresent) FROM cond AS t1 INNER JOIN affymetrixChip AS t2 ON t1.conditionId = t2.conditionId GROUP BY t1.speciesId ORDER BY COUNT(DISTINCT t2.bgeeAffymetrixChipId) DESC" >> $@.tmp
	@echo "Statistics on conditions" >> $@.tmp
	@$(MYSQL) -e "SELECT t1.speciesId, COUNT(DISTINCT t2.bgeeAffymetrixChipId), COUNT(DISTINCT t2.conditionId), COUNT(DISTINCT t1.exprMappedConditionId), COUNT(DISTINCT t1.anatEntityId), COUNT(DISTINCT t1.stageId), COUNT(DISTINCT t1.anatEntityId, t1.stageId, t1.sex), GROUP_CONCAT(DISTINCT t1.sex ORDER BY t1.sex SEPARATOR ', '), GROUP_CONCAT(DISTINCT t1.strain ORDER BY t1.strain SEPARATOR ', ') FROM cond AS t1 INNER JOIN affymetrixChip AS t2 ON t1.conditionId = t2.conditionId GROUP BY t1.speciesId ORDER BY COUNT(DISTINCT t2.bgeeAffymetrixChipId) DESC" >> $@.tmp
	@echo "Same condition information, but for mapped conditions of expression tables" >> $@.tmp
	@$(MYSQL) -e "SELECT t3.speciesId, COUNT(DISTINCT t2.bgeeAffymetrixChipId), COUNT(DISTINCT t3.anatEntityId), COUNT(DISTINCT t3.stageId), COUNT(DISTINCT t3.anatEntityId, t3.stageId, t3.sex), GROUP_CONCAT(DISTINCT t3.sex ORDER BY t3.sex SEPARATOR ', '), GROUP_CONCAT(DISTINCT t3.strain ORDER BY t3.strain SEPARATOR ', ') FROM cond AS t1 INNER JOIN affymetrixChip AS t2 ON t1.conditionId = t2.conditionId INNER JOIN cond AS t3 ON t1.exprMappedConditionId = t3.conditionId GROUP BY t3.speciesId ORDER BY COUNT(DISTINCT t2.bgeeAffymetrixChipId) DESC" >> $@.tmp
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@

TMP_TABLE_CLAUSE   := SELECT t10.expressionId, IF((SELECT COUNT(microarrayExperimentId) FROM microarrayExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.microarrayExperimentCallDirection = 'present' AND t11.microarrayExperimentCallQuality = 'high quality') >= 2, 'pstGold', IF((SELECT COUNT(microarrayExperimentId) FROM microarrayExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.microarrayExperimentCallDirection = 'present' AND t11.microarrayExperimentCallQuality = 'high quality') = 1, 'pstSilverFromHigh', IF((SELECT COUNT(microarrayExperimentId) FROM microarrayExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.microarrayExperimentCallDirection = 'present' AND t11.microarrayExperimentCallQuality = 'poor quality') >= 2, 'pstSilverFromLow', IF((SELECT COUNT(microarrayExperimentId) FROM microarrayExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.microarrayExperimentCallDirection = 'present' AND t11.microarrayExperimentCallQuality = 'poor quality') = 1, 'pstBronze', IF((SELECT COUNT(microarrayExperimentId) FROM microarrayExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.microarrayExperimentCallDirection = 'absent' AND t11.microarrayExperimentCallQuality = 'high quality') >= 2, 'absGold', IF((SELECT COUNT(microarrayExperimentId) FROM microarrayExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.microarrayExperimentCallDirection = 'absent' AND t11.microarrayExperimentCallQuality = 'high quality') = 1, 'absSilverFromHigh', IF((SELECT COUNT(microarrayExperimentId) FROM microarrayExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.microarrayExperimentCallDirection = 'absent' AND t11.microarrayExperimentCallQuality = 'poor quality') >= 2, 'absSilverFromLow', IF((SELECT COUNT(microarrayExperimentId) FROM microarrayExperimentExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.microarrayExperimentCallDirection = 'absent' AND t11.microarrayExperimentCallQuality = 'poor quality') = 1, 'absBronze', 'IMPOSSIBLE')))))))) AS expressionCall FROM microarrayExperimentExpression AS t10

check_conditions: insert_affy
	@echo "Conditions with anat. entity not existing in related species: " > $@.tmp
	@$(MYSQL) -e "select DISTINCT t1.speciesId, t1.conditionId, t1.exprMappedConditionId, t1.anatEntityId, t3.anatEntityName, t1.stageId, t4.stageName, t1.sex, t1.strain from cond as t1 left outer join anatEntityTaxonConstraint as t2 on t1.anatEntityId = t2.anatEntityId and (t2.speciesId is null or t1.speciesId = t2.speciesId) left outer join anatEntity AS t3 ON t3.anatEntityId = t1.anatEntityId LEFT OUTER JOIN stage AS t4 ON t1.stageId = t4.stageId left outer join affymetrixChip as t10 on t1.conditionId = t10.conditionId where t2.anatEntityId is null AND t10.conditionId IS NOT NULL order by t1.speciesId;" >> $@.tmp 2> $@.warn
	@echo "Conditions with dev. stage not existing in related species: " >> $@.tmp
	@$(MYSQL) -e "select DISTINCT t1.speciesId, t1.conditionId, t1.exprMappedConditionId, t1.anatEntityId, t3.anatEntityName, t1.stageId, t4.stageName, t1.sex, t1.strain from cond as t1 left outer join stageTaxonConstraint as t2 on t1.stageId = t2.stageId and (t2.speciesId is null or t1.speciesId = t2.speciesId) left outer join anatEntity AS t3 ON t3.anatEntityId = t1.anatEntityId LEFT OUTER JOIN stage AS t4 ON t1.stageId = t4.stageId left outer join affymetrixChip as t10 on t1.conditionId = t10.conditionId where t2.stageId is null AND t10.conditionId IS NOT NULL order by t1.speciesId;" >> $@.tmp 2>> $@.warn
	@$(MV) $@.tmp $@

insert_expression: check_conditions insert_affy
	# Insert the expression summary
	@perl Data_insertion/insert_affy_expression.pl -bgee=$(BGEECMD) > $@.tmp 2>warnings.$@
	@echo "Searching for incorrect updates of affymetrixProbeset/expression tables (there should be none)" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM affymetrixProbeset WHERE reasonForExclusion = 'not excluded' AND expressionId IS NULL" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM affymetrixProbeset WHERE reasonForExclusion != 'not excluded' AND expressionId IS NOT NULL" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM affymetrixProbeset AS t1 WHERE expressionId IS NOT NULL AND NOT EXISTS(SELECT 1 FROM microarrayExperimentExpression AS t2 WHERE t2.expressionId = t1.expressionId)" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM microarrayExperimentExpression AS t1 WHERE NOT EXISTS (SELECT 1 FROM affymetrixProbeset AS t2 WHERE t2.expressionId = t1.expressionId)" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM microarrayExperimentExpression WHERE (presentHighMicroarrayChipCount + presentLowMicroarrayChipCount + absentHighMicroarrayChipCount + absentLowMicroarrayChipCount) = 0" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM microarrayExperimentExpression WHERE (presentHighMicroarrayChipCount > 0 AND (microarrayExperimentCallDirection != 'present' OR microarrayExperimentCallQuality != 'high quality')) OR (presentHighMicroarrayChipCount = 0 AND presentLowMicroarrayChipCount > 0 AND (microarrayExperimentCallDirection != 'present' OR microarrayExperimentCallQuality != 'poor quality')) OR (presentHighMicroarrayChipCount = 0 AND presentLowMicroarrayChipCount = 0 AND absentHighMicroarrayChipCount > 0 AND (microarrayExperimentCallDirection != 'absent' OR microarrayExperimentCallQuality != 'high quality')) OR (presentHighMicroarrayChipCount = 0 AND presentLowMicroarrayChipCount = 0 AND absentHighMicroarrayChipCount = 0 AND (microarrayExperimentCallDirection != 'absent' OR microarrayExperimentCallQuality != 'poor quality'))" >> $@.tmp
	@$(MYSQL) -e "SELECT t1.* FROM microarrayExperiment AS t1 WHERE NOT EXISTS (SELECT 1 FROM microarrayExperimentExpression AS t2 WHERE t2.microarrayExperimentId = t1.microarrayExperimentId)" >> $@.tmp
	@echo >> $@.tmp
	@echo "Statistics for affymetrixProbeset" >> $@.tmp
	@$(MYSQL) -e "SELECT speciesId, expressionCall, detectionFlag, affymetrixData, COUNT(*) AS affymetrixProbesetCount FROM affymetrixProbeset AS t1 INNER JOIN gene AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId INNER JOIN ($(TMP_TABLE_CLAUSE)) AS tmpTable ON t1.expressionId = tmpTable.expressionId GROUP BY speciesId, expressionCall, detectionFlag, affymetrixData ORDER BY speciesId, expressionCall, detectionFlag, affymetrixData" >> $@.tmp
	@echo "Statistics for expression table" >> $@.tmp
	@$(MYSQL) -e "SELECT speciesId, expressionCall, COUNT(*), COUNT(DISTINCT t1.bgeeGeneId), COUNT(DISTINCT t1.conditionId), COUNT(DISTINCT t3.anatEntityId, t3.stageId), COUNT(DISTINCT t3.anatEntityId, t3.stageId, t3.sex) FROM expression AS t1 INNER JOIN cond AS t3 ON t1.conditionId = t3.conditionId INNER JOIN ($(TMP_TABLE_CLAUSE)) AS tmpTable ON t1.expressionId = tmpTable.expressionId GROUP BY speciesId, expressionCall ORDER BY speciesId, expressionCall" >> $@.tmp
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@

#######################
# External annotation
# They have to do done in Other_species/ where scripts are able to create fake annotation.xls sheets.


$(VERIFICATIONFILE): get_annot maintenance prepare_annot check_mas5 gen_info check_annot normalization check_after check_conditions insert_expression
	# Get MySQL Affymetrix statistics
	#TODO
	@$(MYSQL) -e "SELECT COUNT(*) FROM affymetrixProbeset WHERE 1=1" > $@.tmp
	@echo "TODO statistics per species"
	@$(MV) $@.tmp $@
	# Add the verification file to GIT
	-@$(GIT) add $@
	-@$(GIT) commit $@ -m "Commit verification file for Affymetrix"
	@$(GIT) push

clean:
	-@$(RM) $(VERIFICATIONFILE) new_annotation maintenance get_annot prepare_annot clean_mas5 check_mas5 gen_info check_annot normalization check_after insertion *.tmp *.err *.out check_conditions insert_expression
	-@$(GIT) rm $(VERIFICATIONFILE)
	-@$(GIT) commit $(VERIFICATIONFILE) -m "Removing verification file for Affymetrix"
	@$(GIT) push

# special target to remove species and taxa from database. We make the name long on purpose, and does not
# incude it in the clean target, to not wipe the database accidentally
.PHONY := $(.PHONY) deleteAffy
deleteAffy:
# expression has already been started filling in EST step
	@$(MYSQL) -e "DELETE FROM affymetrixProbeset; DELETE FROM affymetrixChip; DELETE FROM microarrayExperiment; DELETE FROM chipType;"

