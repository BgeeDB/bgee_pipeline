PIPELINEROOT := ../
DIR_NAME := ESTs/unigene/
include $(PIPELINEROOT)Makefile.common
EST_SP_TAXID := 7227, 7955, 8364, 9606, 10090

all: $(VERIFICATIONFILE)


get_annot:
	# Retrieve up-to-date annotation files from https://gitlab.sib.swiss/Bgee/expression-annotations submodule
	@$(CP) $(ANNOT_SUBMODULE_DIR)/Strains/StrainMapping.tsv                                $(STRAIN_MAPPING_FILE)
	@$(CP) $(ANNOT_SUBMODULE_DIR)/$(DIR_NAME)/species_EST_file.txt                         $(SOURCE_FILES_DIR)/$(DIR_NAME)/species_EST_file.txt
	@$(CP) $(ANNOT_SUBMODULE_DIR)/$(DIR_NAME)/annotation_libs_Drosophila_melanogaster.txt  $(SOURCE_FILES_DIR)/$(DIR_NAME)/annotation_libs_Drosophila_melanogaster.txt
	@$(CP) $(ANNOT_SUBMODULE_DIR)/$(DIR_NAME)/annotation_libs_Danio_rerio.txt              $(SOURCE_FILES_DIR)/$(DIR_NAME)/annotation_libs_Danio_rerio.txt
	@$(CP) $(ANNOT_SUBMODULE_DIR)/$(DIR_NAME)/annotation_libs_Homo_sapiens.txt             $(SOURCE_FILES_DIR)/$(DIR_NAME)/annotation_libs_Homo_sapiens.txt
	@$(CP) $(ANNOT_SUBMODULE_DIR)/$(DIR_NAME)/annotation_libs_Mus_musculus.txt             $(SOURCE_FILES_DIR)/$(DIR_NAME)/annotation_libs_Mus_musculus.txt
	@$(CP) $(ANNOT_SUBMODULE_DIR)/$(DIR_NAME)/annotation_libs_Xenopus_tropicalis.txt       $(SOURCE_FILES_DIR)/$(DIR_NAME)/annotation_libs_Xenopus_tropicalis.txt
	@$(CP) $(ANNOT_SUBMODULE_DIR)/ESTs/smirnadb/mapping_libs_smiRNAdb.csv                  $(SOURCE_FILES_DIR)/ESTs/smirnadb/mapping_libs_smiRNAdb.csv
	@$(CP) $(ANNOT_SUBMODULE_DIR)/ESTs/smirnadb/S.csv                                      $(SOURCE_FILES_DIR)/ESTs/smirnadb/S.csv
	@touch $@


unigene: get_annot $(EST_SPECIES_CURATION_FILEPATH)
	# Get species list from Bgee
	@$(MYSQL) -e "SELECT CONCAT(genus, '_', species) FROM species WHERE speciesId IN ($(EST_SP_TAXID))" | grep -v 'CONCAT' >$@.temp
	@if [[ -s $@.temp ]]; then echo -n; else echo "NO species"; exit 2; fi
	# Download UniGene species data set
	-@for sp in `cat $@.temp`; do \
		$(WGET) -nH -r --cut-dirs=3 -P $(SOURCE_FILES_DIR)/ESTs/  ftp://ftp.ncbi.nih.gov/repository/UniGene/$$sp/*.seq.uniq.gz; \
		$(WGET) -nH -r --cut-dirs=3 -P $(SOURCE_FILES_DIR)/ESTs/  ftp://ftp.ncbi.nih.gov/repository/UniGene/$$sp/*.data.gz; \
	done
	# Xenopus tropicalis UniGene file follow the former species name, so rename it for patterns used later
	@$(MV) $(SOURCE_FILES_DIR)/ESTs/Str.seq.uniq.gz $(SOURCE_FILES_DIR)/ESTs/Xt.seq.uniq.gz
	@ln -fs $(SOURCE_FILES_DIR)/ESTs/*.data.gz .
	@$(MV) $@.temp $@

library.report:
	# Download the latest version of library.report
	@$(WGET) 'ftp://ftp.ncbi.nih.gov/repository/UniLib/library.report' -O $(SOURCE_FILES_DIR)/ESTs/$@.tmp && $(MV) $(SOURCE_FILES_DIR)/ESTs/$@.tmp $(SOURCE_FILES_DIR)/ESTs/$@ 2>/dev/null  || rm -f $(SOURCE_FILES_DIR)/ESTs/$@.tmp
	@touch $@

ensembl_cdna:
	# Get species list from Bgee
	@$(MYSQL) -e "SELECT LOWER(CONCAT(genus, '_', species)) FROM species WHERE speciesId IN ($(EST_SP_TAXID))" | grep -v 'CONCAT' >$@.temp
	@if [[ -s $@.temp ]]; then echo -n; else echo "NO species"; exit 2; fi
	# Download all Ensembl cdna in fasta
	-@for sp in `cat $@.temp`; do \
		$(WGET) -nH -r --cut-dirs=5 -P $(SOURCE_FILES_DIR)/ESTs/  ftp://ftp.ensembl.org/pub/release-$(ENSRELEASE)/fasta/$$sp/cdna/*.cdna.all.fa.gz; \
	done
	# !!! NEED disk space to gunzip all those files !!! So gunzip only required ones later
	-@for cdna in $(SOURCE_FILES_DIR)/ESTs/*.cdna.all.fa.gz; do \
		base=$${cdna/**\//}; \
		tr=`echo $$base | sed -e 's@\..*@@; s@^\(.\).*_@\L\1@'`; \
		ln -fs $$cdna $${tr:0:4}.cdna.all.fa.gz; \
	done
	@$(MV) $@.temp $@

#TODO use external mapping for all? Without Ensembl or through UniProt ???
#TODO use better mapper than blastn! (blat, diamond, ???)
#${cdna,} and ${sp:0:2} are valid bash variable manipulations!
ensembl_mapping: ensembl_cdna unigene
	# Test for BLAST+
	@which makeblastdb
	@which blastn
	# Get species list from Bgee
	@$(MYSQL) -e "SELECT CONCAT(SUBSTR(genus, 1, 1), species) FROM species WHERE speciesId IN ($(EST_SP_TAXID))" | grep -v 'CONCAT' >$@.temp
	@if [[ -s $@.temp ]]; then echo -n; else echo "NO species"; exit 2; fi
	# Prepare for mapping
	-@for sp in `cat $@.temp`; do \
		$(GUNZIP) --to-stdout $(SOURCE_FILES_DIR)/ESTs/$${sp:0:2}.seq.uniq.gz >$${sp:0:2}.seq.uniq; \
		perl change_UniGene_headers.pl $${sp:0:2}.seq.uniq > $${sp:0:2}.seq.uniq_new_headers; \
		export cdna=$${sp:0:4}; \
		$(GUNZIP) --to-stdout $${cdna,}.cdna.all.fa.gz > $${cdna,}_raw; \
		perl change_Ensembl_headers.pl $${cdna,}_raw > $${cdna,}_cdna; \
		makeblastdb -dbtype nucl -input_type fasta -in $${cdna,}_cdna -out $${cdna,}_db -title $${cdna,}_db; \
		blastn -dust no -outfmt 6 -db $${cdna,}_db -query $${sp:0:2}.seq.uniq_new_headers -num_threads 4 -evalue 1e-10 -out $${cdna,}_cdna.results; \
		perl extract_results_for_all.pl $${cdna,}_cdna.results > $${cdna,}_mapping; \
		$(RM) $${sp:0:2}.seq.uniq $${sp:0:2}.seq.uniq_new_headers $${cdna,}.cdna.all.fa.gz $${cdna,}_raw $${cdna,}_cdna $${cdna,}_db* $${cdna,}_cdna.results; \
	done
	@$(MV) $@.temp $@


#NOTE This step is frozen because not updated this 2009
smiRNAdb: get_annot $(SMIRNADB_CURATION_MAPFILEPATH) $(SOURCE_FILES_DIR)/ESTs/smirnadb/S.csv
	# Check mapping file
	@if [[ -s $(SMIRNADB_CURATION_MAPFILEPATH) ]]; then echo -n; else echo "NO smiRNAdb mapping file"; exit 3; fi
	# Download smiRNAdb available report files
	@$(WGET) 'http://www.mirz.unibas.ch/Annotation/Report_Human.csv.bz2' -O $(SOURCE_FILES_DIR)/ESTs/smirnadb/Report_Human.csv.bz2.tmp && $(MV) $(SOURCE_FILES_DIR)/ESTs/smirnadb/Report_Human.csv.bz2.tmp $(SOURCE_FILES_DIR)/ESTs/smirnadb/Report_Human.csv.bz2 2>/dev/null  || rm -f $(SOURCE_FILES_DIR)/ESTs/smirnadb/Report_Human.csv.bz2.tmp
	@$(WGET) 'http://www.mirz.unibas.ch/Annotation/Report_Mouse.csv.bz2' -O $(SOURCE_FILES_DIR)/ESTs/smirnadb/Report_Mouse.csv.bz2.tmp && $(MV) $(SOURCE_FILES_DIR)/ESTs/smirnadb/Report_Mouse.csv.bz2.tmp $(SOURCE_FILES_DIR)/ESTs/smirnadb/Report_Mouse.csv.bz2 2>/dev/null  || rm -f $(SOURCE_FILES_DIR)/ESTs/smirnadb/Report_Mouse.csv.bz2.tmp
	@$(WGET) 'http://www.mirz.unibas.ch/Annotation/Report_Fish.csv.bz2'  -O $(SOURCE_FILES_DIR)/ESTs/smirnadb/Report_Fish.csv.bz2.tmp  && $(MV) $(SOURCE_FILES_DIR)/ESTs/smirnadb/Report_Fish.csv.bz2.tmp  $(SOURCE_FILES_DIR)/ESTs/smirnadb/Report_Fish.csv.bz2  2>/dev/null  || rm -f $(SOURCE_FILES_DIR)/ESTs/smirnadb/Report_Fish.csv.bz2.tmp
	@$(WGET) 'http://www.mirz.unibas.ch/Annotation/Report_Fly.csv.bz2'   -O $(SOURCE_FILES_DIR)/ESTs/smirnadb/Report_Fly.csv.bz2.tmp   && $(MV) $(SOURCE_FILES_DIR)/ESTs/smirnadb/Report_Fly.csv.bz2.tmp   $(SOURCE_FILES_DIR)/ESTs/smirnadb/Report_Fly.csv.bz2   2>/dev/null  || rm -f $(SOURCE_FILES_DIR)/ESTs/smirnadb/Report_Fly.csv.bz2.tmp
	# bunzip smiRNAdb report files
	@cd $(SOURCE_FILES_DIR)/ESTs/smirnadb/; bunzip2 -f Report_*.csv.bz2
	@touch $@

miRNA.dat: smiRNAdb
	# Get miRNA.dat file
	@$(WGET) 'ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.dat.gz' -O $(SOURCE_FILES_DIR)/ESTs/$@.gz.tmp && $(MV) $(SOURCE_FILES_DIR)/ESTs/$@.gz.tmp $(SOURCE_FILES_DIR)/ESTs/$@.gz 2>/dev/null  || rm -f $(SOURCE_FILES_DIR)/ESTs/$@.gz.tmp
	@$(GUNZIP) --to-stdout $(SOURCE_FILES_DIR)/ESTs/$@.gz > $@


#NOTE only speciesId from species with EST annotation files are listed here
EST: library.report unigene ensembl_mapping $(EST_SPECIES_CURATION_FILEPATH) $(UBERON_SEX_INFO_FILE_PATH)
	# Insert ESTs info in Bgee
	@$(MYSQL) -e "SELECT speciesId FROM species WHERE speciesId IN ($(EST_SP_TAXID))" | grep -v 'speciesId' >$@.temp
	@if [[ -s $@.temp ]]; then echo -n; else echo "NO species"; exit 3; fi
	@for sp in `cat $@.temp`; do \
		echo $$sp; \
		$(IDMAPPING)  $(IDMAPPINGPORT) &; \
		$(STGMAPPING) $(STGMAPPINGPORT) &; \
		sleep 50; \
		perl insert_est.pl -species=$$sp -bgee=$(BGEECMD) -Aport=$(IDMAPPINGPORT) -Sport=$(STGMAPPINGPORT) -ESTspecies=$(EST_SPECIES_CURATION_FILEPATH) -libReport=$(SOURCE_FILES_DIR)/ESTs/library.report -sex_info=$(UBERON_SEX_INFO_FILE_PATH) >> $@.temp 2>> warnings.$@; \
	done
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.temp $@


#TODO problems with for loop and ID MAPPINGs. Fixed with mappers killed ???
#NOTE only speciesId from species with miRNA annotation files are listed here
miRNA_est: miRNA.dat $(SMIRNADB_CURATION_MAPFILEPATH) $(UBERON_SEX_INFO_FILE_PATH)
	# Insert miRNA ESTs
	@$(MYSQL) -e "SELECT speciesId FROM species WHERE speciesId IN (7227, 7955, 9606, 10090)" | grep -v 'speciesId' >$@.temp
	@if [[ -s $@.temp ]]; then echo -n; else echo "NO species"; exit 4; fi
	@for sp in `cat $@.temp`; do \
		echo $$sp; \
		$(IDMAPPING)  $(IDMAPPINGPORT) &; \
		$(STGMAPPING) $(STGMAPPINGPORT) &; \
		sleep 50; \
		perl insert_miRNA_est.pl -species=$$sp -bgee=$(BGEECMD) -Aport=$(IDMAPPINGPORT) -Sport=$(STGMAPPINGPORT) -smiRNAdb=$(SMIRNADB_CURATION_MAPFILEPATH) -miRNA=miRNA.dat -sex_info=$(UBERON_SEX_INFO_FILE_PATH) >> $@.temp 2>> warnings.$@; \
	done
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.temp $@


#TODO redo this query without estLibraryCallQuality but with pValues!!!
TMP_TABLE_CLAUSE   := SELECT t10.expressionId, IF((SELECT COUNT(estLibraryId) FROM estLibraryExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.estLibraryCallQuality = 'high quality') >= 2, 'pstGold', IF((SELECT COUNT(estLibraryId) FROM estLibraryExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.estLibraryCallQuality = 'high quality') = 1, 'pstSilverFromHigh', IF((SELECT COUNT(estLibraryId) FROM estLibraryExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.estLibraryCallQuality = 'poor quality') >= 2, 'pstSilverFromLow', IF((SELECT COUNT(estLibraryId) FROM estLibraryExpression AS t11 WHERE t11.expressionId = t10.expressionId AND t11.estLibraryCallQuality = 'poor quality') = 1, 'pstBronze', 'IMPOSSIBLE')))) AS expressionCall FROM estLibraryExpression AS t10

check_conditions: EST miRNA_est
	@echo "Conditions with anat. entity not existing in related species: " > $@.tmp
	@$(MYSQL) -e "select DISTINCT t1.speciesId, t1.conditionId, t1.exprMappedConditionId, t1.anatEntityId, t3.anatEntityName, t1.stageId, t4.stageName, t1.sex, t1.strain from cond as t1 left outer join anatEntityTaxonConstraint as t2 on t1.anatEntityId = t2.anatEntityId and (t2.speciesId is null or t1.speciesId = t2.speciesId) left outer join anatEntity AS t3 ON t3.anatEntityId = t1.anatEntityId LEFT OUTER JOIN stage AS t4 ON t1.stageId = t4.stageId left outer join estLibrary as t10 on t1.conditionId = t10.conditionId where t2.anatEntityId is null AND t10.conditionId IS NOT NULL order by t1.speciesId;" >> $@.tmp 2> $@.warn
	@echo "Conditions with dev. stage not existing in related species: " >> $@.tmp
	@$(MYSQL) -e "select DISTINCT t1.speciesId, t1.conditionId, t1.exprMappedConditionId, t1.anatEntityId, t3.anatEntityName, t1.stageId, t4.stageName, t1.sex, t1.strain from cond as t1 left outer join stageTaxonConstraint as t2 on t1.stageId = t2.stageId and (t2.speciesId is null or t1.speciesId = t2.speciesId) left outer join anatEntity AS t3 ON t3.anatEntityId = t1.anatEntityId LEFT OUTER JOIN stage AS t4 ON t1.stageId = t4.stageId left outer join estLibrary as t10 on t1.conditionId = t10.conditionId where t2.stageId is null AND t10.conditionId IS NOT NULL order by t1.speciesId;" >> $@.tmp 2>> $@.warn
	@$(MV) $@.tmp $@

expression_est: check_conditions EST miRNA_est
	# Insert the expression summaries
	@perl insert_expression_est.pl -bgee=$(BGEECMD) > $@.tmp 2>warnings.$@
	@echo "Searching for incorrect updates of expressedSequenceTag/expression tables (there should be none)" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM expressedSequenceTag WHERE expressionId IS NULL" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM expressedSequenceTag AS t1 WHERE NOT EXISTS(SELECT 1 FROM estLibraryExpression AS t2 WHERE t2.expressionId = t1.expressionId)" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM estLibraryExpression AS t1 WHERE NOT EXISTS (SELECT 1 FROM expressedSequenceTag AS t2 WHERE t2.expressionId = t1.expressionId)" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM estLibraryExpression WHERE estCount = 0" >> $@.tmp
	@$(MYSQL) -e "SELECT * FROM estLibraryExpression WHERE estCount >= 7 AND estLibraryCallQuality != 'high quality' OR estCount < 7 AND estLibraryCallQuality != 'poor quality'" >> $@.tmp
	@$(MYSQL) -e "SELECT t1.* FROM estLibrary AS t1 WHERE NOT EXISTS (SELECT 1 FROM estLibraryExpression AS t2 WHERE t2.estLibraryId = t1.estLibraryId)" >> $@.tmp
	@echo >> $@.tmp
	@echo "Statistics for expression table" >> $@.tmp
	@$(MYSQL) -e "SELECT speciesId, expressionCall, COUNT(*), COUNT(DISTINCT t1.bgeeGeneId), COUNT(DISTINCT t1.conditionId), COUNT(DISTINCT t3.anatEntityId, t3.stageId), COUNT(DISTINCT t3.anatEntityId, t3.stageId, t3.sex) FROM expression AS t1 INNER JOIN cond AS t3 ON t1.conditionId = t3.conditionId INNER JOIN ($(TMP_TABLE_CLAUSE)) AS tmpTable ON t1.expressionId = tmpTable.expressionId GROUP BY speciesId, expressionCall ORDER BY speciesId, expressionCall" >> $@.tmp
	@if [[ ! -s warnings.$@ ]]; then $(RM) warnings.$@; fi
	@$(MV) $@.tmp $@


$(VERIFICATIONFILE): EST miRNA_est check_conditions expression_est
	# Check EST data insertion
	@$(MYSQL) -e "SELECT t2.speciesId, COUNT(t1.estId) FROM expressedSequenceTag AS t1 INNER JOIN gene AS t2 ON t1.geneId = t2.geneId GROUP BY t2.speciesId" > $@.tmp
	@$(MYSQL) -e "SELECT speciesId FROM species ORDER BY speciesId" | grep -v 'speciesId' >species
	@for sp in `cat species`; do \
		$(MYSQL) -e "SELECT t2.speciesId, t3.geneBioTypeName, IF(t1.expressionId IS NOT NULL, 1, 0) AS expressionState, COUNT(DISTINCT t1.estId) AS estCount, COUNT(DISTINCT t1.estId)/COUNT(DISTINCT t1.geneId) AS estsByGene FROM expressedSequenceTag AS t1 INNER JOIN gene AS t2 ON t1.geneId = t2.geneId LEFT OUTER JOIN geneBioType AS t3 ON t2.geneBioTypeId = t3.geneBioTypeId WHERE t2.speciesId = $$sp GROUP BY t2.speciesId, t3.geneBioTypeName, expressionState ORDER BY estCount DESC" >> $@.tmp; \
	done
	@$(RM) species
	@$(MV) $@.tmp $@
	# Add the verification file to GIT
	-@$(GIT) add $@
	-@$(GIT) commit $@ -m "Commit verification file for ESTs/miRNA"
	@$(GIT) push

clean:
	-@$(RM) -R $(VERIFICATIONFILE) *.data *.data.gz S.xls Report_*.csv.bz2 Dm.seq.uniq* Dm.seq.uniq_new_headers* my_db* *_cdna.results* \
		*_cdna ensembl_mapping EST miRNA.dat check_conditions expression_est miRNA_est species *.temp smiRNAdb *.tmp unigene warnings.*
	-@$(GIT) rm $(VERIFICATIONFILE)
	-@$(GIT) commit $(VERIFICATIONFILE) -m "Removing verification file for ESTs/miRNA"
	@$(GIT) push

# special target to remove ESTs from database. We make the name long on purpose, and does not
# incude it in the clean target, to not wipe the database accidentally
.PHONY := $(.PHONY) deleteEST
deleteEST:
	@$(MYSQL) -e "DELETE FROM estLibraryToKeyword; DELETE FROM expressedSequenceTag; DELETE FROM estLibrary;"

