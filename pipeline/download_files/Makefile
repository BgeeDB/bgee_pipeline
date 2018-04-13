PIPELINEROOT := ../
DIR_NAME := download_files/
include $(PIPELINEROOT)Makefile.common

all: $(VERIFICATIONFILE)

######################################################
# CREATE VARIABLES FOR INSERTING DOWNLOAD FILES INFO
#######################################################
# Creates the variables needed to insert information about generated download files, 
# using Makefile foreach function.
# The order of the file name prefixes is important, it defines the preferred 
# display order of the groups.
SPECIES_GROUP_VAR_PREFIXES := HUMAN MOUSE ZEBRAFISH FLY C_ELEGANS CHIMPANZEE BONOBO GORILLA MACACA RAT COW PIG HORSE RABBIT DOG CAT CAVY HEDGEHOG OPOSSUM PLATYPUS CHICKEN ANOLIS XENOPUS D_ANA D_MOJ D_PSE D_SIM D_VIR D_YAK

# We will need to replace spaces produce by the Makefile foreach method with $(ENTRY_SEP), 
# we need the following variables: 
empty:=
space:= $(empty) $(empty)

############
# Variables to create a Map species data group names -> species IDs 
###########
# Create a variable storing the parameter group names to species for single species data groups. 
# We replace spaces in display names with __, to be replaced afterwards with an escaped space 
# (because, see next line)
SINGLE_SPE_GROUP_NAME_TO_SPECIES := $(strip $(foreach prefix, $(SPECIES_GROUP_VAR_PREFIXES), $(subst $(space),__,$($(prefix)_DISPLAY_NAME))$(KEY_VAL_SEP)$($(prefix)_ID)))
# Now, we replace the spaces produce by the Makefile 'foreach' method with $(ENTRY_SEP)
SINGLE_SPE_GROUP_NAME_TO_SPECIES := $(subst $(space),$(ENTRY_SEP),$(SINGLE_SPE_GROUP_NAME_TO_SPECIES))
# We put back the spaces in display names
SINGLE_SPE_GROUP_NAME_TO_SPECIES := $(subst __,$(space),$(SINGLE_SPE_GROUP_NAME_TO_SPECIES))

# create the "real" variable MULTI_SPE_GROUP_NAME_TO_SPECIES, if we had data in all necessay species
# TODO: add many more groups, try all possible interesting groupings
# MULTI_SPE_GROUP_NAME_TO_SPECIES := $(HUMAN_DISPLAY_NAME)/$(MOUSE_DISPLAY_NAME)$(KEY_VAL_SEP)$(HUMAN_ID)$(VAL_SEP)$(MOUSE_ID)$(ENTRY_SEP)Murinae$(KEY_VAL_SEP)$(MOUSE_ID)$(VAL_SEP)$(RAT_ID)$(ENTRY_SEP)Primates$(KEY_VAL_SEP)$(PRIMATES_IDS)$(ENTRY_SEP)Theria$(KEY_VAL_SEP)$(THERIA_IDS)$(ENTRY_SEP)Mammalia$(KEY_VAL_SEP)$(MAMMALIA_IDS)$(ENTRY_SEP)Amniota$(KEY_VAL_SEP)$(AMNIOTA_IDS)$(ENTRY_SEP)Tetrapoda$(KEY_VAL_SEP)$(TETRAPODA_IDS)$(ENTRY_SEP)Euteleostomi$(KEY_VAL_SEP)$(EUTELEOSTOMI_IDS)

# create temporary variables because we don't have the multi-species present/absent expression files for now, allowing to include more species
CURRENT_THERIA_IDS := $(HUMAN_ID)$(VAL_SEP)$(CHIMPANZEE_ID)$(VAL_SEP)$(GORILLA_ID)$(VAL_SEP)$(MACACA_ID)$(VAL_SEP)$(RAT_ID)$(VAL_SEP)$(MOUSE_ID)$(VAL_SEP)$(COW_ID)$(VAL_SEP)$(OPOSSUM_ID)
CURRENT_MAMMALIA_IDS := $(CURRENT_THERIA_IDS)$(VAL_SEP)$(PLATYPUS_ID)
CURRENT_AMNIOTA_IDS := $(CURRENT_MAMMALIA_IDS)$(VAL_SEP)$(CHICKEN_ID)
CURRENT_TETRAPODA_IDS := $(CURRENT_AMNIOTA_IDS)$(VAL_SEP)$(XENOPUS_ID)

MULTI_SPE_GROUP_NAME_TO_SPECIES := $(HUMAN_DISPLAY_NAME)/$(MOUSE_DISPLAY_NAME)$(KEY_VAL_SEP)$(HUMAN_ID)$(VAL_SEP)$(MOUSE_ID)$(ENTRY_SEP)Murinae$(KEY_VAL_SEP)$(MOUSE_ID)$(VAL_SEP)$(RAT_ID)$(ENTRY_SEP)$(MACACA_DISPLAY_NAME)/$(CHIMPANZEE_DISPLAY_NAME)$(KEY_VAL_SEP)$(MACACA_ID)$(VAL_SEP)$(CHIMPANZEE_ID)$(ENTRY_SEP)Theria$(KEY_VAL_SEP)$(CURRENT_THERIA_IDS)$(ENTRY_SEP)Mammalia$(KEY_VAL_SEP)$(CURRENT_MAMMALIA_IDS)$(ENTRY_SEP)Amniota$(KEY_VAL_SEP)$(CURRENT_AMNIOTA_IDS)$(ENTRY_SEP)Tetrapoda$(KEY_VAL_SEP)$(CURRENT_TETRAPODA_IDS)


GROUP_NAME_TO_SPECIES := $(SINGLE_SPE_GROUP_NAME_TO_SPECIES)
# FIXME: enable when multi-species file generation will be possible
#$(ENTRY_SEP)$(MULTI_SPE_GROUP_NAME_TO_SPECIES)

############
# Variables to create a Map species data group names -> requested file categories 
###########
EXPR_SIMPLE_CAT := expr_simple
EXPR_COMPLETE_CAT := expr_advanced
DIFF_EXPR_ANAT_SIMPLE_CAT := diff_expr_anatomy_simple
DIFF_EXPR_ANAT_COMPLETE_CAT := diff_expr_anatomy_complete
DIFF_EXPR_DEV_SIMPLE_CAT := diff_expr_dev_simple
DIFF_EXPR_DEV_COMPLETE_CAT := diff_expr_dev_complete
ALL_CALLS_CATS := $(EXPR_SIMPLE_CAT)$(VAL_SEP)$(EXPR_COMPLETE_CAT)
# FIXME: enable when diff expression file generation will be possible
#$(VAL_SEP)$(DIFF_EXPR_ANAT_SIMPLE_CAT)$(VAL_SEP)$(DIFF_EXPR_ANAT_COMPLETE_CAT)$(VAL_SEP)$(DIFF_EXPR_DEV_SIMPLE_CAT)$(VAL_SEP)$(DIFF_EXPR_DEV_COMPLETE_CAT)

RNA_SEQ_ANNOT_CAT := rnaseq_annot
RNA_SEQ_DATA_CAT := rnaseq_data
AFFY_ANNOT_CAT := affy_annot
AFFY_DATA_CAT := affy_data
ALL_RAW_DATA_CATS := $(RNA_SEQ_ANNOT_CAT)$(VAL_SEP)$(RNA_SEQ_DATA_CAT)$(VAL_SEP)$(AFFY_ANNOT_CAT)$(VAL_SEP)$(AFFY_DATA_CAT)
ORTHOLOG_CAT := ortholog

# Create a variable storing the parameter group names to file categories for single species data groups. 
# We replace spaces in display names with __, to be replaced afterwards with an escaped space 
# (because, see next line)
SINGLE_SPE_GROUP_NAME_TO_CAT := $(strip $(foreach prefix, $(SPECIES_GROUP_VAR_PREFIXES), $(subst $(space),__,$($(prefix)_DISPLAY_NAME))$(KEY_VAL_SEP)$(ALL_CALLS_CATS)$(VAL_SEP)$(ALL_RAW_DATA_CATS)))
# Now, we replace the spaces produce by the Makefile 'foreach' method with $(ENTRY_SEP)
SINGLE_SPE_GROUP_NAME_TO_CAT := $(subst $(space),$(ENTRY_SEP),$(SINGLE_SPE_GROUP_NAME_TO_CAT))
# We put back the spaces in display names
SINGLE_SPE_GROUP_NAME_TO_CAT := $(subst __,$(space),$(SINGLE_SPE_GROUP_NAME_TO_CAT))

# Same for multi-species
# create the "real" variable MULTI_SPE_GROUP_NAME_TO_CAT, if we had data in all necessay species
# TODO: add many more groups, try all possible interesting groupings
# MULTI_SPE_GROUP_NAME_TO_CAT := $(HUMAN_DISPLAY_NAME)/$(MOUSE_DISPLAY_NAME)$(KEY_VAL_SEP)$(ALL_CALLS_CATS)$(VAL_SEP)$(ORTHOLOG_CAT)$(ENTRY_SEP)Murinae$(KEY_VAL_SEP)$(ALL_CALLS_CATS)$(VAL_SEP)$(ORTHOLOG_CAT)$(ENTRY_SEP)Primates$(KEY_VAL_SEP)$(ALL_CALLS_CATS)$(VAL_SEP)$(ORTHOLOG_CAT)$(ENTRY_SEP)Theria$(KEY_VAL_SEP)$(ALL_CALLS_CATS)$(VAL_SEP)$(ORTHOLOG_CAT)$(ENTRY_SEP)Mammalia$(KEY_VAL_SEP)$(ALL_CALLS_CATS)$(VAL_SEP)$(ORTHOLOG_CAT)$(ENTRY_SEP)Amniota$(KEY_VAL_SEP)$(ALL_CALLS_CATS)$(VAL_SEP)$(ORTHOLOG_CAT)$(ENTRY_SEP)Tetrapoda$(KEY_VAL_SEP)$(ALL_CALLS_CATS)$(VAL_SEP)$(ORTHOLOG_CAT)$(ENTRY_SEP)Euteleostomi$(KEY_VAL_SEP)$(ALL_CALLS_CATS)$(VAL_SEP)$(ORTHOLOG_CAT)

# create temporary variables because we don't have the multi-species present/absent expression files for now, allowing to include more species
MULTI_SPE_GROUP_NAME_TO_CAT := $(HUMAN_DISPLAY_NAME)/$(MOUSE_DISPLAY_NAME)$(KEY_VAL_SEP)$(DIFF_EXPR_ANAT_SIMPLE_CAT)$(VAL_SEP)$(DIFF_EXPR_ANAT_COMPLETE_CAT)$(VAL_SEP)$(ORTHOLOG_CAT)$(ENTRY_SEP)Murinae$(KEY_VAL_SEP)$(DIFF_EXPR_ANAT_SIMPLE_CAT)$(VAL_SEP)$(DIFF_EXPR_ANAT_COMPLETE_CAT)$(VAL_SEP)$(ORTHOLOG_CAT)$(ENTRY_SEP)$(MACACA_DISPLAY_NAME)/$(CHIMPANZEE_DISPLAY_NAME)$(KEY_VAL_SEP)$(DIFF_EXPR_ANAT_SIMPLE_CAT)$(VAL_SEP)$(DIFF_EXPR_ANAT_COMPLETE_CAT)$(VAL_SEP)$(ORTHOLOG_CAT)$(ENTRY_SEP)Theria$(KEY_VAL_SEP)$(DIFF_EXPR_ANAT_SIMPLE_CAT)$(VAL_SEP)$(DIFF_EXPR_ANAT_COMPLETE_CAT)$(VAL_SEP)$(ORTHOLOG_CAT)$(ENTRY_SEP)Mammalia$(KEY_VAL_SEP)$(DIFF_EXPR_ANAT_SIMPLE_CAT)$(VAL_SEP)$(DIFF_EXPR_ANAT_COMPLETE_CAT)$(VAL_SEP)$(ORTHOLOG_CAT)$(ENTRY_SEP)Amniota$(KEY_VAL_SEP)$(DIFF_EXPR_ANAT_SIMPLE_CAT)$(VAL_SEP)$(DIFF_EXPR_ANAT_COMPLETE_CAT)$(VAL_SEP)$(ORTHOLOG_CAT)$(ENTRY_SEP)Tetrapoda$(KEY_VAL_SEP)$(DIFF_EXPR_ANAT_SIMPLE_CAT)$(VAL_SEP)$(DIFF_EXPR_ANAT_COMPLETE_CAT)$(VAL_SEP)$(ORTHOLOG_CAT)


GROUP_NAME_TO_CAT := $(SINGLE_SPE_GROUP_NAME_TO_CAT)
# FIXME: enable when multi-species file generation will be possible
#$(ENTRY_SEP)$(MULTI_SPE_GROUP_NAME_TO_CAT)

############
# Variables to create a Map species data group names -> replacement string (to replace REPLACE_PATTERN in file names)
###########

# Create a variable storing the parameter group names to replacements for single species data groups. 
# We replace spaces in display names with __, to be replaced afterwards with an escaped space 
# (because, see next line)
SINGLE_SPE_GROUP_NAME_TO_REPLACEMENT := $(strip $(foreach prefix, $(SPECIES_GROUP_VAR_PREFIXES), $(subst $(space),__,$($(prefix)_DISPLAY_NAME))$(KEY_VAL_SEP)$(subst $(space),_,$($(prefix)_SCIENTIFIC_NAME))))
# Now, we replace the spaces produce by the Makefile 'foreach' method with $(ENTRY_SEP)
SINGLE_SPE_GROUP_NAME_TO_REPLACEMENT := $(subst $(space),$(ENTRY_SEP),$(SINGLE_SPE_GROUP_NAME_TO_REPLACEMENT))
# We put back the spaces in display names
SINGLE_SPE_GROUP_NAME_TO_REPLACEMENT := $(subst __,$(space),$(SINGLE_SPE_GROUP_NAME_TO_REPLACEMENT))

MULTI_SPE_GROUP_NAME_TO_REPLACEMENT := $(HUMAN_DISPLAY_NAME)/$(MOUSE_DISPLAY_NAME)$(KEY_VAL_SEP)$(subst $(space),_,$(HUMAN_DISPLAY_NAME))_$(subst $(space),_,$(MOUSE_DISPLAY_NAME))$(ENTRY_SEP)Murinae$(KEY_VAL_SEP)murinae$(ENTRY_SEP)$(MACACA_DISPLAY_NAME)/$(CHIMPANZEE_DISPLAY_NAME)$(KEY_VAL_SEP)$(subst $(space),_,$(MACACA_DISPLAY_NAME))_$(subst $(space),_,$(CHIMPANZEE_DISPLAY_NAME))$(ENTRY_SEP)Theria$(KEY_VAL_SEP)theria$(ENTRY_SEP)Mammalia$(KEY_VAL_SEP)mammalia$(ENTRY_SEP)Amniota$(KEY_VAL_SEP)amniota$(ENTRY_SEP)Tetrapoda$(KEY_VAL_SEP)tetrapoda

GROUP_NAME_TO_REPLACEMENT := $(SINGLE_SPE_GROUP_NAME_TO_REPLACEMENT)
# FIXME enable when multi-species file generation will be possible
#$(ENTRY_SEP)$(MULTI_SPE_GROUP_NAME_TO_REPLACEMENT)

############
# Variables to create a Map category -> file name pattern, for single species data groups, and a Map category -> file name pattern, for multi-species data groups
########### 
REPLACE_PATTERN := {REPLACE}

ZIP_EXT := .zip

SINGLE_SPE_CAT_TO_FILE_NAME_PATTERN := $(EXPR_SIMPLE_CAT)$(KEY_VAL_SEP)$(EXPR_CALLS_PATH)$(REPLACE_PATTERN)$(EXPR_SIMPLE_SUF)$(ZIP_EXT)$(ENTRY_SEP)$(EXPR_SIMPLE_CAT)$(KEY_VAL_SEP)$(EXPR_CALLS_PATH)$(REPLACE_PATTERN)$(EXPR_SIMPLE_DEV_SUF)$(ZIP_EXT)$(ENTRY_SEP)$(EXPR_SIMPLE_CAT)$(KEY_VAL_SEP)$(EXPR_CALLS_PATH)$(REPLACE_PATTERN)$(EXPR_COMPLETE_SUF)$(ZIP_EXT)$(ENTRY_SEP)$(EXPR_COMPLETE_CAT)$(KEY_VAL_SEP)$(EXPR_CALLS_PATH)$(REPLACE_PATTERN)$(EXPR_COMPLETE_DEV_SUF)$(ZIP_EXT)$(ENTRY_SEP)$(RNA_SEQ_ANNOT_CAT)$(KEY_VAL_SEP)$(RNA_SEQ_EXPR_VALUES_PATH)$(REPLACE_PATTERN)/$(REPLACE_PATTERN)$(RNA_SEQ_EXP_LIB_SUF_ZIP)$(ENTRY_SEP)$(RNA_SEQ_DATA_CAT)$(KEY_VAL_SEP)$(RNA_SEQ_EXPR_VALUES_PATH)$(REPLACE_PATTERN)/$(REPLACE_PATTERN)$(RNA_SEQ_COUNT_SUF_ZIP)$(ENTRY_SEP)$(AFFY_ANNOT_CAT)$(KEY_VAL_SEP)$(AFFYMETRIX_EXPR_VALUES_PATH)$(REPLACE_PATTERN)/$(REPLACE_PATTERN)$(AFFYMETRIX_EXP_CHIP_SUF_ZIP)$(ENTRY_SEP)$(AFFY_DATA_CAT)$(KEY_VAL_SEP)$(AFFYMETRIX_EXPR_VALUES_PATH)$(REPLACE_PATTERN)/$(REPLACE_PATTERN)$(AFFYMETRIX_PROBESET_SUF_ZIP)
# FIXME: enable when diff-expr file generation will be possible
#$(DIFF_EXPR_ANAT_SIMPLE_CAT)$(KEY_VAL_SEP)$(DIFF_EXPR_CALLS_PATH)$(REPLACE_PATTERN)$(DIFF_EXPR_ANATOMY_SIMPLE_SUF)$(ZIP_EXT)$(ENTRY_SEP)$(DIFF_EXPR_ANAT_COMPLETE_CAT)$(KEY_VAL_SEP)$(DIFF_EXPR_CALLS_PATH)$(REPLACE_PATTERN)$(DIFF_EXPR_ANATOMY_COMPLETE_SUF)$(ZIP_EXT)$(ENTRY_SEP)$(DIFF_EXPR_DEV_SIMPLE_CAT)$(KEY_VAL_SEP)$(DIFF_EXPR_CALLS_PATH)$(REPLACE_PATTERN)$(DIFF_EXPR_DEV_SIMPLE_SUF)$(ZIP_EXT)$(ENTRY_SEP)$(DIFF_EXPR_DEV_COMPLETE_CAT)$(KEY_VAL_SEP)$(DIFF_EXPR_CALLS_PATH)$(REPLACE_PATTERN)$(DIFF_EXPR_DEV_COMPLETE_SUF)$(ZIP_EXT)$(ENTRY_SEP)


# we don't have all categories for multi-species for now
# MULTI_SPE_CAT_TO_FILE_NAME_PATTERN := $(EXPR_SIMPLE_CAT)$(KEY_VAL_SEP)$(MULTI_SPE_EXPR_CALLS_PATH)$(REPLACE_PATTERN)$(MULTI_SPE_EXPR_SIMPLE_SUF)$(ZIP_EXT)$(ENTRY_SEP)$(EXPR_COMPLETE_CAT)$(KEY_VAL_SEP)$(MULTI_SPE_EXPR_CALLS_PATH)$(REPLACE_PATTERN)$(MULTI_SPE_EXPR_COMPLETE_SUF)$(ZIP_EXT)$(ENTRY_SEP)$(DIFF_EXPR_ANAT_SIMPLE_CAT)$(KEY_VAL_SEP)$(MULTI_SPE_DIFF_EXPR_CALLS_PATH)$(REPLACE_PATTERN)$(MULTI_SPE_DIFF_EXPR_ANATOMY_SIMPLE_SUF)$(ZIP_EXT)$(ENTRY_SEP)$(DIFF_EXPR_ANAT_COMPLETE_CAT)$(KEY_VAL_SEP)$(MULTI_SPE_DIFF_EXPR_CALLS_PATH)$(REPLACE_PATTERN)$(MULTI_SPE_DIFF_EXPR_ANATOMY_COMPLETE_SUF)$(ZIP_EXT)$(ENTRY_SEP)$(DIFF_EXPR_DEV_SIMPLE_CAT)$(KEY_VAL_SEP)$(MULTI_SPE_DIFF_EXPR_CALLS_PATH)$(REPLACE_PATTERN)$(MULTI_SPE_DIFF_EXPR_DEV_SIMPLE_SUF)$(ZIP_EXT)$(ENTRY_SEP)$(DIFF_EXPR_DEV_COMPLETE_CAT)$(KEY_VAL_SEP)$(MULTI_SPE_DIFF_EXPR_CALLS_PATH)$(REPLACE_PATTERN)$(MULTI_SPE_DIFF_EXPR_DEV_COMPLETE_SUF)$(ZIP_EXT)

MULTI_SPE_CAT_TO_FILE_NAME_PATTERN := $(ORTHOLOG_CAT)$(KEY_VAL_SEP)$(ORTHOLOGS_PATH)$(REPLACE_PATTERN)$(ORTHOLOGS_SUF)$(ZIP_EXT)$(ENTRY_SEP)$(DIFF_EXPR_ANAT_SIMPLE_CAT)$(KEY_VAL_SEP)$(MULTI_SPE_DIFF_EXPR_CALL_PATH)$(REPLACE_PATTERN)$(MULTI_SPE_DIFF_EXPR_ANATOMY_SIMPLE_SUF)$(ZIP_EXT)$(ENTRY_SEP)$(DIFF_EXPR_ANAT_COMPLETE_CAT)$(KEY_VAL_SEP)$(MULTI_SPE_DIFF_EXPR_CALL_PATH)$(REPLACE_PATTERN)$(MULTI_SPE_DIFF_EXPR_ANATOMY_COMPLETE_SUF)$(ZIP_EXT)

############
# Variables to create specific files
###########
ORGAN_PARAM := ANAT_ENTITY_ID
STAGE_PARAM := DEV_STAGE_ID

######################################################
# VARIABLES CREATED, INSERT
######################################################

# =========================
# Processed data download files
# =========================
$(OUTPUT_DIR)generate_processed_files:
	@perl generate_ref_expr_files.pl -bgee=$(BGEECMD) -speciesArg=- -affyDir=$(OUTPUT_DIR)affymetrix/ -estDir=- -inSituDir=- -rnaSeqDir=$(OUTPUT_DIR)rna_seq_processed/ -bgeeVersion=bgee_v14 > $@.tmp 2>$@.warnings
	@$(MV) $@.tmp $@

# =========================
# Gene expression call download files
# =========================
#TODO remove usage of each ID and use variable.
$(OUTPUT_DIR)generate_anatomy_call_files:
	@$(JAVA) GenerateBasicExprFile $(D_VIR_ID)$(ENTRY_SEP)$(D_ANA_ID)$(ENTRY_SEP)$(D_YAK_ID)$(ENTRY_SEP)$(D_SIM_ID)$(ENTRY_SEP)$(D_MOJ_ID)$(ENTRY_SEP)$(D_PSE_ID)$(ENTRY_SEP)$(CAVY_ID)$(ENTRY_SEP)$(ANOLIS_ID)$(ENTRY_SEP)$(HEDGEHOG_ID)$(ENTRY_SEP)$(RABBIT_ID)$(ENTRY_SEP)$(DOG_ID)$(ENTRY_SEP)$(CAT_ID)$(ENTRY_SEP)$(HORSE_ID)$(ENTRY_SEP)$(PLATYPUS_ID)$(ENTRY_SEP)$(PIG_ID)$(ENTRY_SEP)$(OPOSSUM_ID)$(ENTRY_SEP)$(CHICKEN_ID)$(ENTRY_SEP)$(GORILLA_ID)$(ENTRY_SEP)$(C_ELEGANS_ID)$(ENTRY_SEP)$(BONOBO_ID)$(ENTRY_SEP)$(CHIMPANZEE_ID)$(ENTRY_SEP)$(XENOPUS_ID)$(ENTRY_SEP)$(COW_ID)$(ENTRY_SEP)$(RAT_ID)$(ENTRY_SEP)$(FLY_ID)$(ENTRY_SEP)$(MACACA_ID)$(ENTRY_SEP)$(ZEBRAFISH_ID)$(ENTRY_SEP)$(MOUSE_ID)$(ENTRY_SEP)$(HUMAN_ID) $(EXPR_SIMPLE_CAT)$(ENTRY_SEP)$(EXPR_COMPLETE_CAT) $(DOWNLOAD_FILE_ROOT_PATH)$(EXPR_CALLS_PATH) $(ORGAN_PARAM) > $@.tmp 2>$@.warnings
	@$(MV) $@.tmp $@

$(OUTPUT_DIR)generate_anatomy_development_call_files:
	@$(JAVA) GenerateBasicExprFile $(D_VIR_ID)$(ENTRY_SEP)$(D_ANA_ID)$(ENTRY_SEP)$(D_YAK_ID)$(ENTRY_SEP)$(D_SIM_ID)$(ENTRY_SEP)$(D_MOJ_ID)$(ENTRY_SEP)$(D_PSE_ID)$(ENTRY_SEP)$(CAVY_ID)$(ENTRY_SEP)$(ANOLIS_ID)$(ENTRY_SEP)$(HEDGEHOG_ID)$(ENTRY_SEP)$(RABBIT_ID)$(ENTRY_SEP)$(DOG_ID)$(ENTRY_SEP)$(CAT_ID)$(ENTRY_SEP)$(HORSE_ID)$(ENTRY_SEP)$(PLATYPUS_ID)$(ENTRY_SEP)$(PIG_ID)$(ENTRY_SEP)$(OPOSSUM_ID)$(ENTRY_SEP)$(CHICKEN_ID)$(ENTRY_SEP)$(GORILLA_ID)$(ENTRY_SEP)$(C_ELEGANS_ID)$(ENTRY_SEP)$(BONOBO_ID)$(ENTRY_SEP)$(CHIMPANZEE_ID)$(ENTRY_SEP)$(XENOPUS_ID)$(ENTRY_SEP)$(COW_ID)$(ENTRY_SEP)$(RAT_ID)$(ENTRY_SEP)$(FLY_ID)$(ENTRY_SEP)$(MACACA_ID)$(ENTRY_SEP)$(ZEBRAFISH_ID)$(ENTRY_SEP)$(MOUSE_ID)$(ENTRY_SEP)$(HUMAN_ID) $(EXPR_SIMPLE_CAT)$(ENTRY_SEP)$(EXPR_COMPLETE_CAT) $(DOWNLOAD_FILE_ROOT_PATH)$(EXPR_CALLS_PATH) $(ORGAN_PARAM)$(ENTRY_SEP)$(STAGE_PARAM) > $@.tmp 2>$@.warnings
	@$(MV) $@.tmp $@

#TODO check it, not used because of usage of different servers
$(OUTPUT_DIR)zip_call_files:
	@find $(DOWNLOAD_FILE_ROOT_PATH)$(EXPR_CALLS_PATH) -type f -not -name "*.zip" -not -name "*.tmp" -execdir zip -9 '{}.zip' '{}' \;

# =========================
# Insert info about download files
# =========================
#TODO add $(OUTPUT_DIR)zip_call_files
#TODO instead of using insert_download_file_info.sql to insert data into downloadFile table, fix InsertSpeciesDataGroups
$(OUTPUT_DIR)insert_download_files_info: $(OUTPUT_DIR)generate_processed_files $(OUTPUT_DIR)generate_anatomy_call_files $(OUTPUT_DIR)generate_anatomy_development_call_files
	@$(MYSQL)  insert_download_file_info.sql > $@.tmp 2>$@.warnings
	#@$(JAVA) InsertSpeciesDataGroups "$(GROUP_NAME_TO_SPECIES)" "$(GROUP_NAME_TO_CAT)" "$(GROUP_NAME_TO_REPLACEMENT)" "$(SINGLE_SPE_CAT_TO_FILE_NAME_PATTERN)" "$(MULTI_SPE_CAT_TO_FILE_NAME_PATTERN)" $(DOWNLOAD_FILE_ROOT_PATH) >> $@.tmp 2>>$@.warnings
	@$(MYSQL) -e "SELECT (SELECT COUNT(*) FROM speciesDataGroup) as \`data group count\`, (SELECT COUNT(distinct speciesId) FROM speciesToDataGroup) as \`species in data group count\`, (SELECT COUNT(*) FROM downloadFile) as \`download file count\`" >> $@.tmp 2>>$@.warnings
	@$(MV) $@.tmp $@

######################################################
######################################################

# =========================
# XRefs from UniProt
# =========================
# IDs of data source Uniprot/SPTREMBL and Uniprot/SWISSPROT in our database 
# (see pipeline/db_creation/insert_data_sources.sql) 
# Very important that these are correct.
TREMBL_ID := 4
SWISSPROT_ID := 5

# XRef file: 
XREF_FILE := $(OUTPUT_DIR)XRefBgee.txt

# The SQL query run to generate the XRefFile, in readable format: 
# 
# SELECT DISTINCT CONCAT(t1.XRefId, '   DR   Bgee; ', uniqueMapping.uniqueGeneId, '; -.') AS xrefRow 
# FROM geneXRef AS t1 
# INNER JOIN (
# -- Some xRefIds are mapped to multiple genes, we retrieve a unique gene ID with expression data
#     SELECT t1.XRefId, (SELECT t4.geneId FROM geneXRef AS t2 
#                        INNER JOIN gene AS t4 ON t2.bgeeGeneId = t4.bgeeGeneId
#                        WHERE t2.XRefId = t1.XRefId 
#                        AND EXISTS (SELECT 1 FROM expression as t3 WHERE t3.bgeeGeneId = t2.bgeeGeneId)
#                        ORDER BY t4.geneId LIMIT 1) AS uniqueGeneId 
#     FROM geneXRef AS t1 
# ) AS uniqueMapping ON t1.XRefId = uniqueMapping.XRefId 
# WHERE t1.dataSourceId IN ($(TREMBL_ID), $(SWISSPROT_ID)) AND t1.XRefId IS NOT NULL AND t1.XRefId != '' AND t1.XRefId NOT LIKE '%\_%' 
# AND uniqueMapping.uniqueGeneId IS NOT NULL 
# ORDER BY uniqueMapping.uniqueGeneId, t1.XRefId;

XREF_SQL_QUERY := SELECT DISTINCT CONCAT(t1.XRefId, '   DR   Bgee; ', uniqueMapping.uniqueGeneId, '; -.') AS xRefRow FROM geneXRef AS t1 INNER JOIN (SELECT t1.XRefId, (SELECT t4.geneId FROM geneXRef AS t2 INNER JOIN gene AS t4 ON t2.bgeeGeneId = t4.bgeeGeneId WHERE t2.XRefId = t1.XRefId AND EXISTS (SELECT 1 FROM expression as t3 WHERE t3.bgeeGeneId = t2.bgeeGeneId) ORDER BY t4.geneId LIMIT 1) AS uniqueGeneId FROM geneXRef AS t1) AS uniqueMapping ON t1.XRefId = uniqueMapping.XRefId WHERE t1.dataSourceId IN ($(TREMBL_ID), $(SWISSPROT_ID)) AND t1.XRefId IS NOT NULL AND t1.XRefId != '' AND t1.XRefId NOT LIKE '%\_%' AND uniqueMapping.uniqueGeneId IS NOT NULL ORDER BY uniqueMapping.uniqueGeneId, t1.XRefId;

$(OUTPUT_DIR)xref_uniprot_info: 
	@$(MYSQL) -e "$(XREF_SQL_QUERY)" > $(XREF_FILE).tmp 2>$@.warnings
# Remove header line
	@sed 1d $(XREF_FILE).tmp > $(XREF_FILE).tmp.noheader && $(MV) $(XREF_FILE).tmp.noheader $(XREF_FILE)
	@$(RM) $(XREF_FILE).tmp
# create a verification file
# First, we display the targeted data sources, so that we can detect 
# incorrect data source IDs used (TREMBL_ID and SWISSPROT_ID)
	@echo "Data sources used for mapping \nIt is very important that it corresponds to TREMBL and SWISSPROT: " > $@.temp
	@$(MYSQL) -e "SELECT * FROM dataSource WHERE dataSourceId IN ($(TREMBL_ID), $(SWISSPROT_ID))" >> $@.temp  2>>$@.warnings
	@echo "\n\nNumber of expected mappings in XRef file: " >> $@.temp
	@$(MYSQL) -e "SELECT COUNT(DISTINCT t1.XRefId) FROM geneXRef AS t1 WHERE dataSourceId IN ($(TREMBL_ID), $(SWISSPROT_ID)) AND EXISTS (SELECT 1 FROM expression AS t2 WHERE t2.bgeeGeneId = t1.bgeeGeneId)" >> $@.temp 2>>$@.warnings
	@echo "\nNumber of lines in actual XRef file: " >> $@.temp
	@wc -l $(XREF_FILE) >> $@.temp
	@echo "\nNumber of mappings per species: " >> $@.temp
	@$(MYSQL) -e "SELECT t3.speciesId, t4.speciesCommonName, CONCAT(t4.genus, ' ', t4.species) AS species, COUNT(DISTINCT t1.XRefId) FROM geneXRef AS t1 INNER JOIN gene AS t3 ON t1.bgeeGeneId = t3.bgeeGeneId INNER JOIN species AS t4 ON t3.speciesId = t4.speciesId WHERE t1.dataSourceId IN ($(TREMBL_ID), $(SWISSPROT_ID)) AND EXISTS (SELECT 1 FROM expression AS t2 WHERE t2.bgeeGeneId = t1.bgeeGeneId) GROUP BY t3.speciesId" >> $@.temp 2>>$@.warnings
	@echo "\n\nNumber of ambiguous mappings per species (we forced the choice of one mapping in these cases): " >> $@.temp
	@$(MYSQL) -e "SELECT multipleGeneCount.speciesId, t4.speciesCommonName, CONCAT(t4.genus, ' ', t4.species) AS species, COUNT(*) FROM (SELECT t1.XRefId, t2.speciesId FROM geneXRef AS t1 INNER JOIN gene AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId WHERE t1.dataSourceId IN ($(TREMBL_ID), $(SWISSPROT_ID)) AND EXISTS (SELECT 1 FROM expression as t3 WHERE t3.bgeeGeneId = t1.bgeeGeneId) GROUP BY t1.XRefId HAVING COUNT(DISTINCT t1.bgeeGeneId) > 1) AS multipleGeneCount INNER JOIN species AS t4 ON multipleGeneCount.speciesId = t4.speciesId  GROUP BY multipleGeneCount.speciesId" >> $@.temp 2>>$@.warnings
	@echo "\nTotal number of ambiguous mappings: " >> $@.temp
	@$(MYSQL) -e "SELECT COUNT(*) FROM (SELECT t1.XRefId FROM geneXRef AS t1 WHERE t1.dataSourceId IN ($(TREMBL_ID), $(SWISSPROT_ID)) AND EXISTS (SELECT 1 FROM expression as t3 WHERE t3.bgeeGeneId = t1.bgeeGeneId) GROUP BY t1.XRefId HAVING COUNT(DISTINCT t1.bgeeGeneId) > 1) AS multipleGeneCount" >> $@.temp 2>>$@.warnings
	@$(MV) $@.temp $@

# =========================
# Download files with rank info
# =========================
$(OUTPUT_DIR)rank_files_info: 
# Create subdirectories if they don't exist
	@mkdir -p $(OUTPUT_DIR)$(RANKS_PER_ANAT); mkdir -p $(OUTPUT_DIR)$(RANKS_PER_CONDITION)
# Launch generation of files.
# First, files with data grouped per anat. entities
	@$(JAVA) GenerateRankFile generateRankFiles - $(UBERONFILEPATH) true $(OUTPUT_DIR)$(RANKS_PER_ANAT) > $@.temp 2>> $@.warnings
# Now, files with all conditions
	@$(JAVA) GenerateRankFile generateRankFiles - $(UBERONFILEPATH) false $(OUTPUT_DIR)$(RANKS_PER_CONDITION) >> $@.temp 2>> $@.warnings
# zip all files
	@find $(OUTPUT_DIR)$(RANKS_PER_ANAT) -type f -not -name "*.zip" -not -name "*.tmp" -execdir zip -9 '{}.zip' '{}' \;
	@find $(OUTPUT_DIR)$(RANKS_PER_CONDITION) -type f -not -name "*.zip" -not -name "*.tmp" -execdir zip -9 '{}.zip' '{}' \;
# generate a step verification file
	@echo "Content of directory $(OUTPUT_DIR)$(RANKS_PER_ANAT)" >> $@.temp
	@ls -lh $(OUTPUT_DIR)$(RANKS_PER_ANAT) >> $@.temp
	@echo -e "\nContent of directory $(OUTPUT_DIR)$(RANKS_PER_CONDITION)" >> $@.temp
	@ls -lh $(OUTPUT_DIR)$(RANKS_PER_CONDITION) >> $@.temp
	@$(MV) $@.temp $@
	
	

# ======================
# Verification file
# ======================
$(VERIFICATIONFILE): $(OUTPUT_DIR)insert_download_files_info $(OUTPUT_DIR)xref_uniprot_info $(OUTPUT_DIR)rank_files_info $(OUTPUT_DIR)generate_processed_files $(OUTPUT_DIR)generate_anatomy_call_files $(OUTPUT_DIR)generate_anatomy_development_call_files $(OUTPUT_DIR)zip_call_files

clean:
	-@$(RM) -R $(VERIFICATIONFILE) $(OUTPUT_DIR)xref_uniprot_info.tmp $(OUTPUT_DIR)xref_uniprot_info $(XREF_FILE).tmp.noheader $(XREF_FILE).tmp $(XREF_FILE) $(OUTPUT_DIR)insert_download_files_info $(OUTPUT_DIR)rank_files_info $(OUTPUT_DIR)generate_processed_files $(OUTPUT_DIR)generate_anatomy_call_files $(OUTPUT_DIR)generate_anatomy_development_call_files $(OUTPUT_DIR)zip_call_files
