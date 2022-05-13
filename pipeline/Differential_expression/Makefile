PIPELINEROOT := ../
include $(PIPELINEROOT)Makefile.common

all: $(VERIFICATIONFILE)


backup_previous_data:
	# Backup previous release data
	@cd $(DIFFEXPRPATH);                        let PREV=$(RELEASE)-1; if [[ ! -e bgee_v$$PREV.tar.gz ]]; then tar cf bgee_v$$PREV.tar --exclude=bgee_v$$PREV.tar --exclude=*.gz *; gzip -9 bgee_v$$PREV.tar; rm -Rf *[^gz]; fi
	@cd $(BIOCONDUCTORDIFF);                    let PREV=$(RELEASE)-1; if [[ ! -e bgee_v$$PREV.tar.gz ]]; then tar cf bgee_v$$PREV.tar --exclude=bgee_v$$PREV.tar --exclude=*.gz *; gzip -9 bgee_v$$PREV.tar; rm -Rf *[^gz]; fi
	@cd $(BIOCONDUCTORTARG);                    let PREV=$(RELEASE)-1; if [[ ! -e bgee_v$$PREV.tar.gz ]]; then tar cf bgee_v$$PREV.tar --exclude=bgee_v$$PREV.tar --exclude=*.gz *; gzip -9 bgee_v$$PREV.tar; rm -Rf *[^gz]; fi
	@cd $(RNASEQDIFFEXPRPATH_DEVANDANAT);       let PREV=$(RELEASE)-1; if [[ ! -e bgee_v$$PREV.tar.gz ]]; then tar cf bgee_v$$PREV.tar --exclude=bgee_v$$PREV.tar --exclude=*.gz *; gzip -9 bgee_v$$PREV.tar; rm -Rf *[^gz]; fi
	#XXX Is this directory used anymore?
	@cd $(RNASEQBIOCONDUCTORDIFF);              let PREV=$(RELEASE)-1; if [[ ! -e bgee_v$$PREV.tar.gz ]]; then tar cf bgee_v$$PREV.tar --exclude=bgee_v$$PREV.tar --exclude=*.gz *; gzip -9 bgee_v$$PREV.tar; rm -Rf *[^gz]; fi
	@cd $(RNASEQBIOCONDUCTORTARG_DEVANDANAT);   let PREV=$(RELEASE)-1; if [[ ! -e bgee_v$$PREV.tar.gz ]]; then tar cf bgee_v$$PREV.tar --exclude=bgee_v$$PREV.tar --exclude=*.gz *; gzip -9 bgee_v$$PREV.tar; rm -Rf *[^gz]; fi
	# Rsync files for backup on sensitive cluster (see wiki for password)
	@echo -e "\t\e[1m\e[37m\e[41mRsync files in '$(AFFYDATAPATH) | $(RNASEQDATAPATH) '\n sensitive machine!\e[m\n"
	@rsync -av --del -f "- /lost+found" --exclude=.svn --exclude=.git $(AFFYDATAPATH)    $(DATALOGIN)@$(DATAHOST):$(DATAPATH)/Affymetrix/   >  $@.tmp
	@rsync -av --del -f "- /lost+found" --exclude=.svn --exclude=.git $(RNASEQDATAPATH)  $(DATALOGIN)@$(DATAHOST):$(DATAPATH)/rna_seq/     >>  $@.tmp
	@$(MV) $@.tmp $@


launch_diff_expression_affy: backup_previous_data
	# Launch differential analysis Affymetrix
	@perl launch_diff_analysis_affy.pl   -bgee=$(BGEECMD) -path_mas5=$(MAS5PATH)   -path_schuster=$(SCHUSTERPATH)  -path_target=$(BIOCONDUCTORTARG)       -path_processed=$(DIFFEXPRPATH)           >$@.tmp 2>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@

launch_diff_expression_rna_seq: backup_previous_data
	# Launch differential analysis RNA-seq
	@perl launch_diff_analysis_rna_seq.pl -bgee=$(BGEECMD) -path_generes=$(RNASEQGENERES) -ensRelease=$(ENSRELEASE) -path_target=$(RNASEQBIOCONDUCTORTARG_DEVANDANAT) -path_processed=$(RNASEQDIFFEXPRPATH_DEVANDANAT)   >>$@.tmp 2>>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@
	
launch_sex_diff_expression_rna_seq:
	# Launch differential analysis RNA-seq
	@perl launch_sex_diff_analysis_rna_seq.pl -bgee=$(BGEECMD) -path_generes=$(RNASEQALLRES) -abundance_file=$(RNASEQABUNDANCEFILE) -path_library=$(RNASEQ_LIB_FILEPATH) -path_library_worm=$(RNASEQ_LIB_FILEPATH_WORM) -path_target=$(RNASEQBIOCONDUCTORTARG_SEX) -path_processed=$(RNASEQDIFFEXPRPATH_SEX)   >>$@.tmp 2>>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@

generate_files_sex_diff_expression_rna_seq: # launch_sex_diff_expression_rna_seq
	# Generate download files for sex differential expression from RNA-Seq
	@perl generate_files_sex_diff_rna_seq.pl -bgee=$(BGEECMD) -path_target=$(RNASEQBIOCONDUCTORTARG_SEX) -path_processed=$(RNASEQDIFFEXPRPATH_SEX) -path_download_file=$(RNASEQDIFFEXPRPATH_SEX)   >>$@.tmp 2>>warnings.$@
	@find $(RNASEQDIFFEXPRPATH_SEX) -maxdepth 1 -type f  -name '*sex-single.tsv' -exec gzip --verbose --best --force {} \; 2>$@.tmp
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@

check_diff_expression: launch_diff_expression
	# Check results
	@-ls -l $(DIFFEXPRPATH)/*/*.out.PROB 2>/dev/null                       >$@.tmp
	@-grep -H -c 'proc.time' $(DIFFEXPRPATH)/*/*.log       | grep ':0$$'  >>$@.tmp
	@-ls -l $(RNASEQDIFFEXPRPATH_DEVANDANAT)/*/*.out.PROB 2>/dev/null                >>$@.tmp
	@-grep -H -c 'proc.time' $(RNASEQDIFFEXPRPATH_DEVANDANAT)/*/*.log | grep ':0$$'  >>$@.tmp
	# Check warnings
	@-grep -i 'warning' -A 5 $(DIFFEXPRPATH)/*/*.log                      >>$@.tmp
	@-grep -i 'warning' -A 5 $(RNASEQDIFFEXPRPATH_DEVANDANAT)/*/*.log                >>$@.tmp
	@$(MV) $@.tmp $@


fill_summary_tables: check_diff_expression
	# Fill the deaAffymetrixProbesetSummary table
	@perl insert_diff_affy.pl   -bgee=$(BGEECMD) -path_target=$(BIOCONDUCTORTARG)       -path_processed=$(DIFFEXPRPATH)        >$@.tmp 2>warnings.$@
	# Check AffymetrixP p-values
	@$(MYSQL) -e "SELECT 'At least one deaAffymetrixProbesetSummary with incorrect p-value' AS errorMessage, deaAffymetrixProbesetSummary.* FROM deaAffymetrixProbesetSummary WHERE deaRawPValue < 0 OR deaRawPValue > 1 LIMIT 1" 2>>warnings.$@
# The sub-query identified probesetSummaryIds that should be excluded for a given deaId, because the probeset was always seen as 'not expressed'
# in all conditions of the related DEA. This is done by using a GROUP_CONCAT: if the probeset is always seen 'not expressed' in the analysis,
# the GROUP_CONCAT will have the value '0'; if it is always seen as 'expressed' in the analysis, the GROUP_CONCAT will have the value '1';
# if there is a mixture, GROUP_CONCAT will have the value '0,1' or '1,0'; So, we simply identify probeset/DEA with GROUP_CONCAT having value '0'.
	@$(MYSQL) -e "UPDATE deaAffymetrixProbesetSummary AS t1 INNER JOIN deaSampleGroup AS t2 ON t1.deaSampleGroupId = t2.deaSampleGroupId INNER JOIN (SELECT t4.deaAffymetrixProbesetSummaryId, t3.deaId FROM deaSampleGroup AS t3 INNER JOIN deaAffymetrixProbesetSummary AS t4 ON t3.deaSampleGroupId = t4.deaSampleGroupId GROUP BY t4.deaAffymetrixProbesetSummaryId, t3.deaId HAVING GROUP_CONCAT(DISTINCT (IF (differentialExpressionAffymetrixData = 'no data' OR differentialExpressionAffymetrixData = 'not expressed', 0, 1))) = '0') AS t5 ON t1.deaAffymetrixProbesetSummaryId = t5.deaAffymetrixProbesetSummaryId AND t2.deaId = t5.deaId SET t1.reasonForExclusion = 'not expressed'" >>$@.tmp
	# Fill the deaRNASeqSummary table
	@echo >>$@.tmp
	@perl insert_diff_rna_seq.pl -bgee=$(BGEECMD) -path_target=$(RNASEQBIOCONDUCTORTARG_DEVANDANAT) -path_processed=$(RNASEQDIFFEXPRPATH_DEVANDANAT) >>$@.tmp 2>>warnings.$@
	# Check RNASeq p-values
	@$(MYSQL) -e "SELECT 'At least one deaRNASeqSummary with incorrect p-value' AS errorMessage, deaRNASeqSummary.* FROM deaRNASeqSummary WHERE deaRawPValue < 0 OR deaRawPValue > 1 LIMIT 1" 2>>warnings.$@
# The sub-query identified geneSummaryIds that should be excluded for a given deaId, because the gene was always seen as 'not expressed'
# in all conditions of the related DEA. This is done by using a GROUP_CONCAT: if the gene is always seen 'not expressed' in the analysis,
# the GROUP_CONCAT will have the value '0'; if it is always seen as 'expressed' in the analysis, the GROUP_CONCAT will have the value '1';
# if there is a mixture, GROUP_CONCAT will have the value '0,1' or '1,0'; So, we simply identify gene/DEA with GROUP_CONCAT having value '0'.
	@$(MYSQL) -e "UPDATE deaRNASeqSummary AS t1 INNER JOIN deaSampleGroup AS t2 ON t1.deaSampleGroupId = t2.deaSampleGroupId INNER JOIN (SELECT t4.geneSummaryId, t3.deaId FROM deaSampleGroup AS t3 INNER JOIN deaRNASeqSummary AS t4 ON t3.deaSampleGroupId = t4.deaSampleGroupId GROUP BY t4.geneSummaryId, t3.deaId HAVING GROUP_CONCAT(DISTINCT (IF (differentialExpressionRNASeqData = 'no data' OR differentialExpressionRNASeqData = 'not expressed', 0, 1))) = '0') AS t5 ON t1.geneSummaryId = t5.geneSummaryId AND t2.deaId = t5.deaId SET t1.reasonForExclusion = 'not expressed'" >>$@.tmp
	#check analyses with no data inserted
	@$(MYSQL) -e "SELECT 'DIFF ANALYSIS WITH NO DATA INSERTED', differentialExpressionAnalysis.* FROM differentialExpressionAnalysis WHERE deaId IN (SELECT DISTINCT t1.deaId FROM differentialExpressionAnalysis AS t1 INNER JOIN deaSampleGroup AS t2 ON t1.deaId = t2.deaId WHERE IF(t1.microarrayExperimentId IS NOT NULL, (SELECT t4.speciesId FROM gene AS t4 WHERE t4.geneId = (SELECT t3.geneId FROM deaAffymetrixProbesetSummary AS t3 WHERE t2.deaSampleGroupId = t3.deaSampleGroupId LIMIT 1)), (SELECT t4.speciesId FROM gene AS t4 WHERE t4.geneId = (SELECT t3.geneSummaryId FROM deaRNASeqSummary AS t3 WHERE t2.deaSampleGroupId = t3.deaSampleGroupId LIMIT 1))) is null)" 2>>warnings.$@
	# Check analyses with not enough conditions or wrong total replicate number
	@$(MYSQL) -e "SELECT 'INCORRECT CONDITION/SAMPLE COUNT', deaId, COUNT(DISTINCT deaSampleGroupId) AS 'sample group count', COUNT(DISTINCT anatEntityId, stageId) AS 'condition count', ((SELECT COUNT(DISTINCT bgeeAffymetrixChipId) FROM deaSampleGroupToAffymetrixChip AS t2 INNER JOIN deaSampleGroup AS t3 ON t2.deaSampleGroupId = t3.deaSampleGroupId WHERE t3.deaId = t1.deaId) + (SELECT COUNT(DISTINCT rnaSeqLibraryId) FROM deaSampleGroupToRnaSeqLibrary AS t2 INNER JOIN deaSampleGroup AS t3 ON t2.deaSampleGroupId = t3.deaSampleGroupId WHERE t3.deaId = t1.deaId)) AS 'chip/library count' FROM deaSampleGroup AS t1 GROUP BY t1.deaId HAVING `sample group count` < 3 OR `sample group count` != `condition count` OR `chip/library count` < (`sample group count` * 2)" 2>>warnings.$@
	# Check conditions with not enough replicates
	@$(MYSQL) -e "SELECT 'INCORRECT CONDITION/SAMPLE COUNT', deaId, deaSampleGroupId, ((SELECT COUNT(DISTINCT bgeeAffymetrixChipId) FROM deaSampleGroupToAffymetrixChip AS t2 WHERE t2.deaSampleGroupId = t1.deaSampleGroupId) + (SELECT COUNT(DISTINCT rnaSeqLibraryId) FROM deaSampleGroupToRnaSeqLibrary AS t2 WHERE t2.deaSampleGroupId = t1.deaSampleGroupId)) AS 'chip/library count' FROM deaSampleGroup AS t1 HAVING `chip/library count` < 2" 2>>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@


fill_dea_table: fill_summary_tables
	# Fill the differentialExpression table
	@perl insert_diff_expression.pl -bgee=$(BGEECMD) -type=affy    >$@.tmp 2>warnings.$@
	@echo >>$@.tmp
	@perl insert_diff_expression.pl -bgee=$(BGEECMD) -type=rnaseq >>$@.tmp 2>>warnings.$@
	@echo >>$@.tmp
	# UPDATE SUMMARY DATA INTO DIFF EXPRESSION TABLE
	# Consistent counts
	@$(MYSQL) -e "UPDATE differentialExpression AS t0 INNER JOIN (SELECT t1.differentialExpressionId, COUNT(DISTINCT t3.deaId) AS consistentDEACount FROM differentialExpression AS t1 INNER JOIN deaAffymetrixProbesetSummary AS t2 ON t1.differentialExpressionId = t2.differentialExpressionId AND IF(t2.differentialExpressionAffymetrixData = 'no diff expression', 'no diff expression', IF(t2.foldChange >= 0, 'over-expression', 'under-expression')) = t1.diffExprCallAffymetrix INNER JOIN deaSampleGroup AS t3 ON t2.deaSampleGroupId = t3.deaSampleGroupId WHERE t2.reasonForExclusion = 'not excluded' AND t2.differentialExpressionAffymetrixData NOT IN ('no data', 'not expressed') AND t1.diffExprCallAffymetrix NOT IN ('no data', 'not expressed') GROUP BY t1.differentialExpressionId) AS tempTable ON t0.differentialExpressionId = tempTable.differentialExpressionId SET t0.consistentDEACountAffymetrix = tempTable.consistentDEACount" 2>>warnings.$@
	@$(MYSQL) -e "UPDATE differentialExpression AS t0 INNER JOIN (SELECT t1.differentialExpressionId, COUNT(DISTINCT t3.deaId) AS consistentDEACount FROM differentialExpression AS t1 INNER JOIN deaRNASeqSummary AS t2 ON t1.differentialExpressionId = t2.differentialExpressionId AND IF(t2.differentialExpressionRNASeqData = 'no diff expression', 'no diff expression', IF(t2.foldChange >= 0, 'over-expression', 'under-expression')) = t1.diffExprCallRNASeq INNER JOIN deaSampleGroup AS t3 ON t2.deaSampleGroupId = t3.deaSampleGroupId WHERE t2.reasonForExclusion = 'not excluded' AND t2.differentialExpressionRNASeqData NOT IN ('no data', 'not expressed') AND t1.diffExprCallRNASeq NOT IN ('no data', 'not expressed') GROUP BY t1.differentialExpressionId) AS tempTable ON t0.differentialExpressionId = tempTable.differentialExpressionId SET t0.consistentDEACountRNASeq = tempTable.consistentDEACount" 2>>warnings.$@
	# Inconsistent counts
	@$(MYSQL) -e "UPDATE differentialExpression AS t0 INNER JOIN (SELECT t1.differentialExpressionId, COUNT(DISTINCT t3.deaId) AS inconsistentDEACount FROM differentialExpression AS t1 INNER JOIN deaAffymetrixProbesetSummary AS t2 ON t1.differentialExpressionId = t2.differentialExpressionId AND IF(t2.differentialExpressionAffymetrixData = 'no diff expression', 'no diff expression', IF(t2.foldChange >= 0, 'over-expression', 'under-expression')) != t1.diffExprCallAffymetrix INNER JOIN deaSampleGroup AS t3 ON t2.deaSampleGroupId = t3.deaSampleGroupId WHERE t2.reasonForExclusion = 'not excluded' AND t2.differentialExpressionAffymetrixData NOT IN ('no data', 'not expressed') AND t1.diffExprCallAffymetrix NOT IN ('no data', 'not expressed') GROUP BY t1.differentialExpressionId) AS tempTable ON t0.differentialExpressionId = tempTable.differentialExpressionId SET t0.inconsistentDEACountAffymetrix = tempTable.inconsistentDEACount" 2>>warnings.$@
	@$(MYSQL) -e "UPDATE differentialExpression AS t0 INNER JOIN (SELECT t1.differentialExpressionId, COUNT(DISTINCT t3.deaId) AS inconsistentDEACount FROM differentialExpression AS t1 INNER JOIN deaRNASeqSummary AS t2 ON t1.differentialExpressionId = t2.differentialExpressionId AND IF(t2.differentialExpressionRNASeqData = 'no diff expression', 'no diff expression', IF(t2.foldChange >= 0, 'over-expression', 'under-expression')) != t1.diffExprCallRNASeq INNER JOIN deaSampleGroup AS t3 ON t2.deaSampleGroupId = t3.deaSampleGroupId WHERE t2.reasonForExclusion = 'not excluded' AND t2.differentialExpressionRNASeqData NOT IN ('no data', 'not expressed') AND t1.diffExprCallRNASeq NOT IN ('no data', 'not expressed') GROUP BY t1.differentialExpressionId) AS tempTable ON t0.differentialExpressionId = tempTable.differentialExpressionId SET t0.inconsistentDEACountRNASeq = tempTable.inconsistentDEACount" 2>>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@


$(VERIFICATIONFILE): check_diff_expression fill_summary_tables fill_dea_table
	# Check Differential Expression data insertion and report statistics per species and geneBioTypeName
	@$(RM) $@.tmp
	# Overview
	@echo "Overview" > $@.tmp
	@$(MYSQL) -e "SELECT speciesId, COUNT(deaId) FROM (SELECT DISTINCT t1.deaId, IF(t1.microarrayExperimentId IS NOT NULL, (SELECT t4.speciesId FROM gene AS t4 WHERE t4.geneId = (SELECT t3.geneId FROM deaAffymetrixProbesetSummary AS t3 WHERE t2.deaSampleGroupId = t3.deaSampleGroupId LIMIT 1)), (SELECT t4.speciesId FROM gene AS t4 WHERE t4.geneId = (SELECT t3.geneSummaryId FROM deaRNASeqSummary AS t3 WHERE t2.deaSampleGroupId = t3.deaSampleGroupId LIMIT 1))) AS speciesId FROM differentialExpressionAnalysis AS t1 INNER JOIN deaSampleGroup AS t2 ON t1.deaId = t2.deaId) AS tempTable GROUP BY speciesId" >> $@.tmp
	@$(MYSQL) -e "SELECT deaId, COUNT(DISTINCT deaSampleGroupId) AS 'sample count', COUNT(DISTINCT anatEntityId, stageId) AS 'condition count', ((SELECT COUNT(DISTINCT bgeeAffymetrixChipId) FROM deaSampleGroupToAffymetrixChip AS t2 INNER JOIN deaSampleGroup AS t3 ON t2.deaSampleGroupId = t3.deaSampleGroupId WHERE t3.deaId = t1.deaId) + (SELECT COUNT(DISTINCT rnaSeqLibraryId) FROM deaSampleGroupToRnaSeqLibrary AS t2 INNER JOIN deaSampleGroup AS t3 ON t2.deaSampleGroupId = t3.deaSampleGroupId WHERE t3.deaId = t1.deaId)) AS 'chip/library count' FROM deaSampleGroup AS t1 GROUP BY t1.deaId" >> $@.tmp
	# Check that all results not excluded have a differentialExpressionId corresponding to their gene/organ/stage/comparison factor
	@echo "Check that all results not excluded have a differentialExpressionId corresponding to their gene/organ/stage/comparison factor" >> $@.tmp
	@$(MYSQL) -e "SELECT 'At least one non-excluded deaAffymetrixProbesetSummary not inserted in table differentialExpression' AS errorMessage, t1.* FROM deaAffymetrixProbesetSummary AS t1 WHERE t1.reasonForExclusion = 'not excluded' AND t1.differentialExpressionId IS NULL LIMIT 1" >> $@.tmp
	@echo >>$@.tmp
	@$(MYSQL) -e "SELECT 'At least one non-excluded deaRNASeqSummary not inserted in table differentialExpression' AS errorMessage, t1.* FROM deaRNASeqSummary AS t1 WHERE t1.reasonForExclusion = 'not excluded' AND t1.differentialExpressionId IS NULL LIMIT 1" >> $@.tmp
	@echo >>$@.tmp
	@$(MYSQL) -e "SELECT 'Incorrect differentialExpressionId assigned to deaAffymetrixProbesetSummary' AS errorMessage, t1.*, t4.* FROM deaAffymetrixProbesetSummary AS t1 INNER JOIN deaSampleGroup AS t2 ON t1.deaSampleGroupId = t2.deaSampleGroupId INNER JOIN differentialExpressionAnalysis AS t3 ON t2.deaId = t3.deaId INNER JOIN differentialExpression AS t4 ON t4.comparisonFactor = t3.comparisonFactor AND t4.geneId = t1.geneId AND t4.anatEntityId = t2.anatEntityId AND t4.stageId = t2.stageId WHERE t1.reasonForExclusion = 'not excluded' AND t1.differentialExpressionId != t4.differentialExpressionId LIMIT 1" >> $@.tmp
	@echo >>$@.tmp
	@$(MYSQL) -e "SELECT 'Incorrect differentialExpressionId assigned to deaRNASeqSummary' AS errorMessage, t1.*, t4.* FROM deaRNASeqSummary AS t1 INNER JOIN deaSampleGroup AS t2 ON t1.deaSampleGroupId = t2.deaSampleGroupId INNER JOIN differentialExpressionAnalysis AS t3 ON t2.deaId = t3.deaId INNER JOIN differentialExpression AS t4 ON t4.comparisonFactor = t3.comparisonFactor AND t4.geneId = t1.geneSummaryId AND t4.anatEntityId = t2.anatEntityId AND t4.stageId = t2.stageId WHERE t1.reasonForExclusion = 'not excluded' AND t1.differentialExpressionId != t4.differentialExpressionId LIMIT 1" >> $@.tmp
	@echo >>$@.tmp
	# Check data that should not have been inserted, because no data for any data type
	@echo "Check data that should not have been inserted, because no data for any data type" >>$@.tmp
	@$(MYSQL) -e "SELECT 'At least one line in table differentialExpression with no data for all data types' AS errorMessage, differentialExpression.* FROM differentialExpression WHERE (diffExprCallAffymetrix = 'no data' OR diffExprCallAffymetrix = 'not expressed') AND (diffExprCallRNASeq = 'no data' OR diffExprCallRNASeq = 'not expressed') LIMIT 1" >> $@.tmp
	@echo >>$@.tmp
	# Check error of comparison factor for insertion into diff expression table
	@echo "Check error of comparison factor for insertion into diff expression table" >>$@.tmp
	@$(MYSQL) -e "SELECT 'At least one line in table differentialExpression with incorrect comparison factor' AS errorMessage, t1.* FROM differentialExpression AS t1 INNER JOIN deaAffymetrixProbesetSummary AS t2 ON t1.differentialExpressionId = t2.differentialExpressionId INNER JOIN deaSampleGroup AS t3 ON t2.deaSampleGroupId = t3.deaSampleGroupId INNER JOIN differentialExpressionAnalysis AS t4 ON t3.deaId = t4.deaId WHERE t4.comparisonFactor != t1.comparisonFactor LIMIT 1" >> $@.tmp
	@echo >>$@.tmp
	@$(MYSQL) -e "SELECT 'At least one line in table differentialExpression with incorrect comparison factor' AS errorMessage, t1.* FROM differentialExpression AS t1 INNER JOIN deaRNASeqSummary AS t2 ON t1.differentialExpressionId = t2.differentialExpressionId INNER JOIN deaSampleGroup AS t3 ON t2.deaSampleGroupId = t3.deaSampleGroupId INNER JOIN differentialExpressionAnalysis AS t4 ON t3.deaId = t4.deaId WHERE t4.comparisonFactor != t1.comparisonFactor LIMIT 1" >> $@.tmp
	@echo >>$@.tmp
	# Check correct insertion into diff expression table
	@echo "Check correct insertion into diff expression table" >> $@.tmp
	@$(MYSQL) -e "SELECT t5.speciesCommonName, t1.deaId, t1.comparisonFactor, t1.microarrayExperimentId, IF(t1.microarrayExperimentId IS NOT NULL, (SELECT GROUP_CONCAT(DISTINCT t10.chipTypeId) FROM affymetrixChip AS t10 INNER JOIN deaSampleGroupToAffymetrixChip AS t20 ON t10.bgeeAffymetrixChipId = t20.bgeeAffymetrixChipId INNER JOIN deaSampleGroup AS t30 ON t20.deaSampleGroupId = t30.deaSampleGroupId WHERE t30.deaId = t1.deaId GROUP BY t30.deaId), NULL) AS chipTypeId, t1.rnaSeqExperimentId, COUNT(DISTINCT t2.anatEntityId, t2.stageId) numberOfConditions, (SELECT t6.anatEntityName FROM anatEntity AS t6 WHERE t6.anatEntityId = t2.anatEntityId LIMIT 1) AS anatEntity, GROUP_CONCAT(t7.stageName ORDER BY t7.stageLeftBound) AS developmentalStage FROM differentialExpressionAnalysis AS t1 INNER JOIN deaSampleGroup AS t2 ON t1.deaId = t2.deaId INNER JOIN species AS t5 ON t5.speciesId = IF(t1.microarrayExperimentId IS NOT NULL, (SELECT t4.speciesId FROM gene AS t4 WHERE t4.geneId = (SELECT t3.geneId FROM deaAffymetrixProbesetSummary AS t3 WHERE t2.deaSampleGroupId = t3.deaSampleGroupId LIMIT 1)), (SELECT t4.speciesId FROM gene AS t4 WHERE t4.geneId = (SELECT t3.geneSummaryId FROM deaRNASeqSummary AS t3 WHERE t2.deaSampleGroupId = t3.deaSampleGroupId LIMIT 1))) INNER JOIN stage AS t7 ON t7.stageId = t2.stageId WHERE t1.comparisonFactor = 'development' GROUP BY t1.deaId ORDER BY t5.speciesCommonName, numberOfConditions DESC" >>$@.tmp
	@$(MYSQL) -e "SELECT t5.speciesCommonName, t1.deaId, t1.comparisonFactor, t1.microarrayExperimentId, IF(t1.microarrayExperimentId IS NOT NULL, (SELECT GROUP_CONCAT(DISTINCT t10.chipTypeId) FROM affymetrixChip AS t10 INNER JOIN deaSampleGroupToAffymetrixChip AS t20 ON t10.bgeeAffymetrixChipId = t20.bgeeAffymetrixChipId INNER JOIN deaSampleGroup AS t30 ON t20.deaSampleGroupId = t30.deaSampleGroupId WHERE t30.deaId = t1.deaId GROUP BY t30.deaId), NULL) AS chipTypeId, t1.rnaSeqExperimentId, COUNT(DISTINCT t2.anatEntityId, t2.stageId) numberOfConditions, (SELECT t6.stageName FROM stage AS t6 WHERE t6.stageId = t2.stageId LIMIT 1) AS stageId, GROUP_CONCAT(t7.anatEntityName ORDER BY t7.anatEntityId) AS anatEntity FROM differentialExpressionAnalysis AS t1 INNER JOIN deaSampleGroup AS t2 ON t1.deaId = t2.deaId INNER JOIN species AS t5 ON t5.speciesId = IF(t1.microarrayExperimentId IS NOT NULL, (SELECT t4.speciesId FROM gene AS t4 WHERE t4.geneId = (SELECT t3.geneId FROM deaAffymetrixProbesetSummary AS t3 WHERE t2.deaSampleGroupId = t3.deaSampleGroupId LIMIT 1)), (SELECT t4.speciesId FROM gene AS t4 WHERE t4.geneId = (SELECT t3.geneSummaryId FROM deaRNASeqSummary AS t3 WHERE t2.deaSampleGroupId = t3.deaSampleGroupId LIMIT 1))) INNER JOIN anatEntity AS t7 ON t7.anatEntityId = t2.anatEntityId WHERE t1.comparisonFactor = 'anatomy' GROUP BY t1.deaId ORDER BY t5.speciesCommonName, numberOfConditions DESC" >>$@.tmp
	@$(MYSQL) -e "SELECT t3.speciesCommonName, comparisonFactor, diffExprCallAffymetrix, diffExprAffymetrixData, diffExprCallRNASeq, diffExprRNASeqData, COUNT(*) FROM differentialExpression AS t1 INNER JOIN gene AS t2 ON t1.geneId = t2.geneId INNER JOIN species AS t3 ON t2.speciesId = t3.speciesId GROUP BY t3.speciesId, comparisonFactor, diffExprCallAffymetrix, diffExprAffymetrixData, diffExprCallRNASeq, diffExprRNASeqData ORDER BY t3.speciesCommonName, comparisonFactor, diffExprCallAffymetrix, diffExprCallRNASeq, diffExprAffymetrixData, diffExprRNASeqData" >> $@.tmp
	@echo >>$@.tmp
	@$(MYSQL) -e "SELECT t3.speciesCommonName, t2.geneBioTypeName, t10.diffExprCallAffymetrix, t10.diffExprCallRNASeq, COUNT(DISTINCT t1.geneId) AS geneCount, COUNT(DISTINCT t10.differentialExpressionId) diffCallCount FROM gene AS t1 LEFT OUTER JOIN geneBioType AS t2 ON t1.geneBioTypeId = t2.geneBioTypeId LEFT OUTER JOIN differentialExpression AS t10 ON t1.geneId = t10.geneId INNER JOIN species AS t3 ON t1.speciesId = t3.speciesId GROUP BY t1.speciesId, t2.geneBioTypeName, t10.diffExprCallAffymetrix, t10.diffExprCallRNASeq ORDER BY t1.speciesId, diffCallCount DESC, geneCount DESC" >> $@.tmp
	@echo >>$@.tmp
	# Check if the number of successful analysis (*.log files in $(DIFFEXPRPATH)/*/ & $(RNASEQDIFFEXPRPATH_DEVANDANAT)/*/ with proc.time()) is the same than the one in the db
	@echo "Check if the number of successful analysis (*.log files in $(DIFFEXPRPATH)/*/ & $(RNASEQDIFFEXPRPATH_DEVANDANAT)/*/ with proc.time()) is the same than the one in the db" >>$@.tmp
	@grep '^> proc.time()$$' $(DIFFEXPRPATH)/*/*.log $(RNASEQDIFFEXPRPATH_DEVANDANAT)/*/*.log | wc -l >> $@.tmp
	@tail -1 $@.tmp
	@$(MYSQL) -e "SELECT COUNT(*) FROM differentialExpressionAnalysis" | tail -1 >>$@.tmp
	@tail -1 $@.tmp
	@if [[ `tail -1 $@.tmp` -ne `tail -2 $@.tmp | head -1` ]]; then echo "Number of $(DIFFEXPRPATH)/*/*.out & $(RNASEQDIFFEXPRPATH_DEVANDANAT)/*/*.out files is NOT the same than the number of rows in the differentialExpressionAnalysis table" >>$@.tmp | tail -1 $@.tmp; fi
	@$(MV) $@.tmp $@

clean:
	-@$(RM) -R $(VERIFICATIONFILE) backup_previous_data launch_diff_expression check_diff_expression  fill_summary_tables   fill_dea_table  *.tmp  warnings.*

# special target to remove Differential Expression data from database. We make the name long on purpose, and does not
# incude it in the clean target, to not wipe the database accidentally
.PHONY := $(.PHONY) deleteDiffExpr
deleteDiffExpr:
	@$(MYSQL) -e "DELETE FROM deaAffymetrixProbesetSummary; DELETE FROM deaRNASeqSummary; DELETE FROM deaSampleGroupToAffymetrixChip; DELETE FROM deaSampleGroupToRnaSeqLibrary; DELETE FROM deaSampleGroup; DELETE FROM differentialExpressionAnalysis; DELETE FROM differentialExpression;"
