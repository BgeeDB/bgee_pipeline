PIPELINEROOT := ../
DIR_NAME := In_situ/
include $(PIPELINEROOT)Makefile.common

all: $(VERIFICATIONFILE)

get_annot:
	@$(GIT) submodule init
	@$(GIT) submodule update
	# Retrieve up-to-date annotation files from https://gitlab.sib.swiss/Bgee/expression-annotations submodule
	@$(CP) $(ANNOT_SUBMODULE_DIR)/Strains/StrainMapping.tsv  $(STRAIN_MAPPING_FILE)
	@touch $@


#TODO it looks there is no strain info
#TODO get Ensembl mapping directly from ZebrafishMine???
map_zfin: get_annot
	# Get data from ZebrafishMine
	@python3 ZFIN/zebra_query.py  > $(SOURCE_FILES_DIR)$(DIR_NAME)ZebrafishMine.data 2>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@touch $@

tsv_zfin: map_zfin $(SOURCE_FILES_DIR)$(DIR_NAME)ZebrafishMine.data
	# Start in between stages socket
	@$(INBETWEENSTAGES) 7955 $(INBETWEENSTAGESPORT) &
	# Start anatomy id mapping socket
	@$(IDMAPPING) $(IDMAPPINGPORT) &
	@sleep 50
	# Prepare ZFIN tsv for insertion
	@perl ZFIN/prepare_data.pl -bgee=$(BGEECMD) -data=$(SOURCE_FILES_DIR)$(DIR_NAME)ZebrafishMine.data -Sport=$(INBETWEENSTAGESPORT) -Aport=$(IDMAPPINGPORT) >$@.tmp 2>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@


#TODO it looks there is no strain info
map_mgi: get_annot
	# Get MGI main RNA in situ data for wild-type mice with Ensembl gene mapped
	@perl MGI/get_mousemine_data.pl > $(SOURCE_FILES_DIR)$(DIR_NAME)MouseMine.data 2>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@touch $@

tsv_mgi: map_mgi $(SOURCE_FILES_DIR)$(DIR_NAME)MouseMine.data
	# Start in between stages socket
	@$(INBETWEENSTAGES) 10090 $(INBETWEENSTAGESPORT) &
	# Start anatomy id mapping socket
	@$(IDMAPPING) $(IDMAPPINGPORT) &
	@sleep 50
	# Prepare MGI tsv for insertion
	@perl MGI/prepare_data.pl  -bgee=$(BGEECMD) -data=$(SOURCE_FILES_DIR)$(DIR_NAME)MouseMine.data -Sport=$(INBETWEENSTAGESPORT) -Aport=$(IDMAPPINGPORT) >$@.tmp 2>warnings.$@
	@gzip -f -9 $(SOURCE_FILES_DIR)$(DIR_NAME)MouseMine.data
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@


#NOTE it looks there is no strain info
#TODO It looks to have in situ data in XenMine now but it is extremely slow
#FIXME issue with FTP passive mode to get GeneExpression_tropicalis.txt
tsv_xenbase: get_annot
	# Start in between stages socket
	@$(INBETWEENSTAGES) 8364 $(INBETWEENSTAGESPORT) &
	# Start anatomy id mapping socket
	@$(IDMAPPING) $(IDMAPPINGPORT) &
	@sleep 50
	# Prepare XenBase tsv for insertion
	@perl Xenbase/get_in_situ_xenbase.pl -bgee=$(BGEECMD) -Sport=$(INBETWEENSTAGESPORT) -Aport=$(IDMAPPINGPORT) -src_dir=$(SOURCE_FILES_DIR)$(DIR_NAME)  >$@.tmp  2>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@


#NOTE it looks there is no strain info, only 'wt' or 'transgenic line'
data_BDGP: get_annot
	# Download up-to-date BDGP in situ data
	@$(WGET) 'http://insitu.fruitfly.org/insitu-mysql-dump/insitu.sql.gz'       -O $(SOURCE_FILES_DIR)$(DIR_NAME)insitu.sql.gz.tmp       && $(MV) $(SOURCE_FILES_DIR)$(DIR_NAME)insitu.sql.gz.tmp       $(SOURCE_FILES_DIR)$(DIR_NAME)insitu.sql.gz 2>/dev/null  || rm -f $(SOURCE_FILES_DIR)$(DIR_NAME)insitu.sql.gz.tmp
	@$(GUNZIP) $(SOURCE_FILES_DIR)$(DIR_NAME)insitu.sql.gz
	@$(WGET) 'http://insitu.fruitfly.org/insitu-mysql-dump/insitu_annot.csv.gz' -O $(SOURCE_FILES_DIR)$(DIR_NAME)insitu_annot.csv.gz.tmp && $(MV) $(SOURCE_FILES_DIR)$(DIR_NAME)insitu_annot.csv.gz.tmp $(SOURCE_FILES_DIR)$(DIR_NAME)insitu_annot.csv.gz 2>/dev/null  || rm -f $(SOURCE_FILES_DIR)$(DIR_NAME)insitu_annot.csv.gz.tmp
	@$(GUNZIP) $(SOURCE_FILES_DIR)$(DIR_NAME)insitu_annot.csv.gz
	# Create BDGP local db
	@$(MYSQL) -e "DROP DATABASE IF EXISTS bdgp; CREATE DATABASE bdgp CHARACTER SET utf8; SHOW WARNINGS;" > $@.tmp
	@$(MYSQLNODBNAME) bdgp < $(SOURCE_FILES_DIR)$(DIR_NAME)insitu.sql >> $@.tmp 2>&1
	@gzip -9 $(SOURCE_FILES_DIR)$(DIR_NAME)insitu.sql
	# Stage info
	@echo -e "#count\tstage_id\torgan name" > BDGP/annot_info
	@cut -f4,5 -d, ../../source_files/In_situ/insitu_annot.csv | sort | uniq -c | sort -n -k1 | sed -e 's/,"/\t/; s/"$$//; s/^ *\([0-9][0-9]*\) /\1\t/' >> BDGP/annot_info
	@$(MV) $@.tmp $@

map_BDGP: data_BDGP $(SOURCE_FILES_DIR)$(DIR_NAME)$(BDGP2FBBT_MAPPING_FILE) BDGP/annot_info
	# Put the new mapping back in the annotator file
	@perl BDGP/complete_mapping_from_bdgp_to_bgee.pl -bgee=$(BGEECMD) -bdgp=$(BDGPCMD) -annot_info=BDGP/annot_info -old=$(SOURCE_FILES_DIR)$(DIR_NAME)$(BDGP2FBBT_MAPPING_FILE) -new=BDGP/$(BDGP2FBBT_MAPPING_NEW)  >$@.tmp 2>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@echo -e "After checking with annotators, replace $(SOURCE_FILES_DIR)$(DIR_NAME)$(BDGP2FBBT_MAPPING_FILE) by BDGP/$(BDGP2FBBT_MAPPING_NEW)"
	@$(CP) BDGP/$(BDGP2FBBT_MAPPING_NEW) $(SOURCE_FILES_DIR)$(DIR_NAME)$(BDGP2FBBT_MAPPING_FILE)
	@$(MV) $@.tmp $@

tsv_BDGP: map_BDGP $(SOURCE_FILES_DIR)$(DIR_NAME)$(BDGP2FBBT_MAPPING_FILE) $(SOURCE_FILES_DIR)$(DIR_NAME)$(STAGECORRESP_FILE) $(SOURCE_FILES_DIR)$(DIR_NAME)insitu_annot.csv
	@echo "Be sure $(SOURCE_FILES_DIR)$(DIR_NAME)$(BDGP2FBBT_MAPPING_FILE) has been checked and completed by annotators"
	# Start in between stages socket
	@$(INBETWEENSTAGES) 7227 $(INBETWEENSTAGESPORT) &
	# Start anatomy id mapping socket
	@$(IDMAPPING) $(IDMAPPINGPORT) &
	@sleep 50
	# Prepare BDGP tsv for insertion
	@perl BDGP/prepare_data.pl  -bgee=$(BGEECMD) -bdgp=$(BDGPCMD) -Sport=$(INBETWEENSTAGESPORT) -Aport=$(IDMAPPINGPORT) -mapping=$(SOURCE_FILES_DIR)$(DIR_NAME)$(BDGP2FBBT_MAPPING_FILE) -annotation=$(SOURCE_FILES_DIR)$(DIR_NAME)insitu_annot.csv -stagecorresp=$(SOURCE_FILES_DIR)$(DIR_NAME)$(STAGECORRESP_FILE)  >$@.out 2>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
#TODO Data should be extracted from FlyMine in the future!
	@$(MV) $@.out log.$@


#TODO it looks there is no strain info
map_flybase: get_annot
	# Get FlyBase mapping
	@export PERL5LIB=FlyBase:$$PERL5LIB ; perl FlyBase/gene_expression_dump.pl FB_public  2>warnings.$@  | uniq >$(SOURCE_FILES_DIR)$(DIR_NAME)Flybase.data.tmp
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $(SOURCE_FILES_DIR)$(DIR_NAME)Flybase.data.tmp $(SOURCE_FILES_DIR)$(DIR_NAME)Flybase.data
	@touch $@

tsv_flybase: map_flybase $(SOURCE_FILES_DIR)$(DIR_NAME)Flybase.data
	# Start in between stages socket
	@$(INBETWEENSTAGES) 7227 $(INBETWEENSTAGESPORT) &
	# Start anatomy id mapping socket
	@$(IDMAPPING) $(IDMAPPINGPORT) &
	@sleep 50
	# Prepare FlyBase tsv for insertion
	@perl FlyBase/prepare_data.pl  -bgee=$(BGEECMD) -Sport=$(INBETWEENSTAGESPORT) -Aport=$(IDMAPPINGPORT) -data=$(SOURCE_FILES_DIR)$(DIR_NAME)Flybase.data >$@.tmp 2>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
#TODO Data should be extracted from FlyMine in the future???
	@$(MV) $@.tmp $@


#TODO fix strain mapping - with Anne's file - on strain json and tsv files. Currently no C. elegans strain mapping in Anne's file
map_wormbase: get_annot
	# Get WormBase data
	#NOTE THINK to update WormBase dump file to last version
	@$(WGET) 'ftp://caltech.wormbase.org/pub/wormbase/expr_dump/$(WORMBASE_EXPR_DUMP).20201213' -O $(SOURCE_FILES_DIR)$(DIR_NAME)$(WORMBASE_EXPR_DUMP).tmp && $(MV) $(SOURCE_FILES_DIR)$(DIR_NAME)$(WORMBASE_EXPR_DUMP).tmp $(SOURCE_FILES_DIR)$(DIR_NAME)$(WORMBASE_EXPR_DUMP).ori 2>/dev/null  || rm -f $(SOURCE_FILES_DIR)$(DIR_NAME)$(WORMBASE_EXPR_DUMP).tmp
	# Get strain mapping from https://wormbase.org/species/all/strain/ natural isolates
	@$(WGET) 'https://wormbase.org/rest/widget/strain/all/natural_isolates?download=1&content-type=application%2Fjson' -O $(SOURCE_FILES_DIR)$(DIR_NAME)wormbase_natural_strains.json.tmp && $(MV) $(SOURCE_FILES_DIR)$(DIR_NAME)wormbase_natural_strains.json.tmp $(SOURCE_FILES_DIR)$(DIR_NAME)wormbase_natural_strains.json 2>/dev/null  || rm -f $(SOURCE_FILES_DIR)$(DIR_NAME)wormbase_natural_strains.json.tmp
	# Get complementary strain mapping from WormMine
	@perl WormBase/get_worm_wild_strains.pl >$(SOURCE_FILES_DIR)$(DIR_NAME)wormbase_wild_isolate.tsv
	# Map strain names
	@perl -e 'use lib ".."; use Utils; Utils::map_strain_names("$(SOURCE_FILES_DIR)$(DIR_NAME)$(WORMBASE_EXPR_DUMP).ori", "$(STRAIN_MAPPING_FILE)")' >$(SOURCE_FILES_DIR)$(DIR_NAME)$(WORMBASE_EXPR_DUMP)
	@$(RM) $(SOURCE_FILES_DIR)$(DIR_NAME)$(WORMBASE_EXPR_DUMP).ori
	@touch $@

#NOTE there are other worm species there! Keep only C. elegans!
tsv_wormbase: map_wormbase $(SOURCE_FILES_DIR)$(DIR_NAME)$(WORMBASE_EXPR_DUMP) $(SOURCE_FILES_DIR)$(DIR_NAME)wormbase_natural_strains.json $(SOURCE_FILES_DIR)$(DIR_NAME)wormbase_wild_isolate.tsv
	# Start in between stages socket
	@$(INBETWEENSTAGES) 6239 $(INBETWEENSTAGESPORT) &
	# Start anatomy id mapping socket
	@$(IDMAPPING) $(IDMAPPINGPORT) &
	@sleep 50
	# Prepare WormBase tsv for insertion
# See dump at ftp://caltech.wormbase.org/pub/wormbase/expr_dump/
	@perl WormBase/extract_insitu_data.pl  -bgee=$(BGEECMD) -Aport=$(IDMAPPINGPORT) -Sport=$(INBETWEENSTAGESPORT) -wormb_data=$(SOURCE_FILES_DIR)$(DIR_NAME)$(WORMBASE_EXPR_DUMP) -strain_map=$(SOURCE_FILES_DIR)$(DIR_NAME)wormbase_natural_strains.json -strain_map2=$(SOURCE_FILES_DIR)$(DIR_NAME)wormbase_wild_isolate.tsv >$@.tmp 2>warnings.$@
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
#TODO Data should be extracted from WormMine in the future!
	@$(MV) $@.tmp $@


insert_insitu: tsv_zfin  tsv_mgi  tsv_xenbase  tsv_BDGP  tsv_flybase  tsv_wormbase $(UBERON_SEX_INFO_FILE_PATH)
	# Insert all in situ data from well-formatted tsv
	@$(RM) $@.tmp warnings.$@
	@for tsv in `ls -1 tsv_*`; do \
		echo $$tsv; \
		>&2 echo $$tsv; \
		perl insert_in_situ.pl -bgee=$(BGEECMD) -tsv=$$tsv -sex_info=$(UBERON_SEX_INFO_FILE_PATH) >> $@.tmp 2>> warnings.$@; \
	done
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@


check_conditions: insert_insitu
	@echo "Conditions with anat. entity not existing in related species: " > $@.tmp
	@$(MYSQL) -e "select DISTINCT t1.speciesId, t1.conditionId, t1.exprMappedConditionId, t1.anatEntityId, t3.anatEntityName, t1.stageId, t4.stageName, t1.sex, t1.strain from cond as t1 left outer join anatEntityTaxonConstraint as t2 on t1.anatEntityId = t2.anatEntityId and (t2.speciesId is null or t1.speciesId = t2.speciesId) left outer join anatEntity AS t3 ON t3.anatEntityId = t1.anatEntityId LEFT OUTER JOIN stage AS t4 ON t1.stageId = t4.stageId left outer join inSituSpot as t10 on t1.conditionId = t10.conditionId where t2.anatEntityId is null AND t10.conditionId IS NOT NULL order by t1.speciesId;" >> $@.tmp 2> $@.warn
	@echo "Conditions with dev. stage not existing in related species: " >> $@.tmp
	@$(MYSQL) -e "select DISTINCT t1.speciesId, t1.conditionId, t1.exprMappedConditionId, t1.anatEntityId, t3.anatEntityName, t1.stageId, t4.stageName, t1.sex, t1.strain from cond as t1 left outer join stageTaxonConstraint as t2 on t1.stageId = t2.stageId and (t2.speciesId is null or t1.speciesId = t2.speciesId) left outer join anatEntity AS t3 ON t3.anatEntityId = t1.anatEntityId LEFT OUTER JOIN stage AS t4 ON t1.stageId = t4.stageId left outer join inSituSpot as t10 on t1.conditionId = t10.conditionId where t2.stageId is null AND t10.conditionId IS NOT NULL order by t1.speciesId;" >> $@.tmp 2>> $@.warn
	@$(MV) $@.tmp $@

#TODO replace inSituData by p-values
insert_expression: check_conditions insert_insitu
	# Insert the expression summary
	@perl insert_expression_in_situ.pl -bgee=$(BGEECMD) > $@.tmp 2>warnings.$@
	@echo "Searching for incorrect updates of inSituSpot/expression tables (there should be none)" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM inSituSpot WHERE reasonForExclusion = 'not excluded' AND expressionId IS NULL" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM inSituSpot WHERE reasonForExclusion != 'not excluded' AND expressionId IS NOT NULL" >> $@.tmp
	@echo >> $@.tmp
	@echo "Statistics for inSituSpot" >> $@.tmp
	@$(MYSQL) -e "SELECT speciesId, detectionFlag, inSituData, COUNT(*) AS inSituSpotCount FROM inSituSpot AS t1 INNER JOIN gene AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId GROUP BY speciesId, detectionFlag, inSituData ORDER BY speciesId, detectionFlag, inSituData" >> $@.tmp
	@echo "Statistics for expression table" >> $@.tmp
	@$(MYSQL) -e "SELECT speciesId, COUNT(*), COUNT(DISTINCT t1.bgeeGeneId), COUNT(DISTINCT t1.conditionId), COUNT(DISTINCT t3.anatEntityId, t3.stageId), COUNT(DISTINCT t3.anatEntityId, t3.stageId, t3.sex) FROM expression AS t1 INNER JOIN cond AS t3 ON t1.conditionId = t3.conditionId GROUP BY speciesId ORDER BY speciesId" >> $@.tmp
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@


$(VERIFICATIONFILE): tsv_zfin  tsv_mgi  tsv_xenbase  tsv_BDGP  tsv_flybase  tsv_wormbase  insert_insitu  check_conditions insert_expression
	# Check in situ data insertion and report statistics per dataSource
	@$(RM) $@.tmp
	@for src in `ls -1 tsv_*`; do \
		$(MYSQL) -e "SELECT t6.dataSourceName, t3.geneBioTypeName, IF(t1.expressionId IS NOT NULL, 1, 0) AS expressionState, COUNT(DISTINCT t1.inSituSpotId) AS spotCount, COUNT(DISTINCT t1.inSituSpotId)/COUNT(DISTINCT t1.bgeeGeneId) AS spotsByGene FROM inSituSpot AS t1 INNER JOIN gene AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId LEFT OUTER JOIN geneBioType AS t3 ON t2.geneBioTypeId = t3.geneBioTypeId LEFT OUTER JOIN inSituEvidence AS t4 ON t1.inSituEvidenceId=t4.inSituEvidenceId LEFT OUTER JOIN inSituExperiment AS t5 ON t4.inSituExperimentId=t5.inSituExperimentId LEFT OUTER JOIN dataSource AS t6 ON t6.dataSourceId=t5.dataSourceId WHERE t6.dataSourceName='$${src/tsv_/}' GROUP BY t6.dataSourceName, t3.geneBioTypeName, expressionState ORDER BY spotCount DESC" >> $@.tmp; \
		$(MYSQL) -e "SELECT COUNT(*) FROM inSituSpot AS t1 INNER JOIN inSituEvidence AS t2 ON t1.inSituEvidenceId=t2.inSituEvidenceId LEFT OUTER JOIN inSituExperiment AS t3 ON t2.inSituExperimentId=t3.inSituExperimentId LEFT OUTER JOIN dataSource AS t4 ON t3.dataSourceId=t4.dataSourceId WHERE t4.dataSourceName='$${src/tsv_/}'" >> $@.tmp; \
		echo >> $@.tmp; \
	done
	@$(MV) $@.tmp $@
# add the verification file to git
	-@$(GIT) add $@
	-@$(GIT) commit $@ -m "Commit verification file for In Situ"
	-@$(GIT) push

clean:
	-@$(RM) -R $(VERIFICATIONFILE) get_annot *_zfin  *_mgi  *_xenbase  *_BDGP  *_flybase  *.tmp  warnings.* *.out log.* check_conditions insert_expression
	-@$(GIT) rm $(VERIFICATIONFILE)
	-@$(GIT) commit $(VERIFICATIONFILE) -m "Removing verification file for In Situ"
	-@$(GIT) push

# special target to remove in situ from database. We make the name long on purpose, and does not
# incude it in the clean target, to not wipe the database accidentally
.PHONY := $(.PHONY) deleteInSitu
deleteInSitu:
	@$(MYSQL) -e "DELETE FROM inSituSpot; DELETE FROM inSituEvidence; DELETE FROM inSituExperiment;"

