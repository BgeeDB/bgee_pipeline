PIPELINEROOT := ../
include $(PIPELINEROOT)Makefile.common

all: $(VERIFICATIONFILE)

CONVERTER     := xlscat
CONVERTEROPT  := -s '\t' -S
EXTRACLEANING := perl -i -pe 's/\t +/\t/g; s/ +\t/\t/g; s/  +/ /g'

# Check that all required somewhere are available
global_check:
# clients/interpretors
	@which java           > $@.tmp
	@which perl          >> $@.tmp
	@which mysql         >> $@.tmp
	@which R             >> $@.tmp
	@which python3       >> $@.tmp
# Ensembl API
	@perl -MBio::EnsEMBL::Registry  -e 1 >> $@.tmp
# other Perl modules
	@perl -MArchive::Tar            -e 1 >> $@.tmp
	@perl -MArchive::Zip            -e 1 >> $@.tmp
	@perl -MArray::Utils            -e 1 >> $@.tmp
	@perl -MBio::SeqIO              -e 1 >> $@.tmp
	@perl -MCpanel::JSON::XS        -e 1 >> $@.tmp
	@perl -MData::Dumper            -e 1 >> $@.tmp
	@perl -MDBI                     -e 1 >> $@.tmp
	@perl -MDigest::SHA             -e 1 >> $@.tmp
	@perl -MExporter                -e 1 >> $@.tmp
	@perl -MFile::Basename          -e 1 >> $@.tmp
	@perl -MFile::Copy              -e 1 >> $@.tmp
	@perl -MFile::Find              -e 1 >> $@.tmp
	@perl -MFile::Path              -e 1 >> $@.tmp
	@perl -MFile::Slurp             -e 1 >> $@.tmp
	@perl -MFile::Spec              -e 1 >> $@.tmp
	@perl -MFindBin                 -e 1 >> $@.tmp
	@perl -MGetopt::Long            -e 1 >> $@.tmp
	@perl -MHTTP::Request           -e 1 >> $@.tmp
	@perl -MIO::Compress::Gzip      -e 1 >> $@.tmp
	@perl -MIO::Socket              -e 1 >> $@.tmp
	@perl -MIO::Uncompress::Gunzip  -e 1 >> $@.tmp
	@perl -MLWP::UserAgent          -e 1 >> $@.tmp
	@perl -MLWP::Simple             -e 1 >> $@.tmp
	@perl -MList::Compare           -e 1 >> $@.tmp
	@perl -MList::MoreUtils         -e 1 >> $@.tmp
	@perl -MList::Util              -e 1 >> $@.tmp
	@perl -MParallel::ForkManager   -e 1 >> $@.tmp
	@perl -MParallel::Loops         -e 1 >> $@.tmp
	@perl -MSort::Naturally         -e 1 >> $@.tmp
	@perl -MSpreadsheet::Read       -e 1 >> $@.tmp
	@perl -MSpreadsheet::ParseExcel -e 1 >> $@.tmp
	@perl -MText::CSV               -e 1 >> $@.tmp
	@perl -MSpreadsheet::ParseXLSX  -e 1 >> $@.tmp
	@perl -MSpreadsheet::WriteExcel -e 1 >> $@.tmp
	@perl -MTest::More              -e 1 >> $@.tmp
	@perl -MTime::HiRes             -e 1 >> $@.tmp
	@perl -MTime::localtime         -e 1 >> $@.tmp
	@perl -MWebservice::InterMine   -e 1 >> $@.tmp
	@perl -MXML::Fast               -e 1 >> $@.tmp
# R libraries
	@R -e 'library("affy")'              >> $@.tmp  2>/dev/null
	@R -e 'library("AnnotationDbi")'     >> $@.tmp  2>/dev/null
	@R -e 'library("BgeeCall")'          >> $@.tmp  2>/dev/null
	@R -e 'library("Biobase")'           >> $@.tmp  2>/dev/null
	@R -e 'library("biomaRt")'           >> $@.tmp  2>/dev/null
	@R -e 'library("Biostrings")'        >> $@.tmp  2>/dev/null
	@R -e 'library("BUSpaRse")'          >> $@.tmp  2>/dev/null
	@R -e 'library("data.table")'        >> $@.tmp  2>/dev/null
	@R -e 'library("devtools")'          >> $@.tmp  2>/dev/null
	@R -e 'library("dplyr")'             >> $@.tmp  2>/dev/null
	@R -e 'library("DropletUtils")'      >> $@.tmp  2>/dev/null
	@R -e 'library("edgeR")'             >> $@.tmp  2>/dev/null
	@R -e 'library("fdrtool")'           >> $@.tmp  2>/dev/null
	@R -e 'library("forcats")'           >> $@.tmp  2>/dev/null
	@R -e 'library("gcrma")'             >> $@.tmp  2>/dev/null
	@R -e 'library("genefilter")'        >> $@.tmp  2>/dev/null
	@R -e 'library("GenomicFeatures")'   >> $@.tmp  2>/dev/null
	@R -e 'library("ggExtra")'           >> $@.tmp  2>/dev/null
	@R -e 'library("gghighlight")'       >> $@.tmp  2>/dev/null
	@R -e 'library("ggplot2")'           >> $@.tmp  2>/dev/null
	@R -e 'library("gridExtra")'         >> $@.tmp  2>/dev/null
	@R -e 'library("HCAExplorer")'       >> $@.tmp  2>/dev/null
	@R -e 'library("HelpersMG")'         >> $@.tmp  2>/dev/null
	@R -e 'library("LaplacesDemon")'     >> $@.tmp  2>/dev/null
	@R -e 'library("lattice")'           >> $@.tmp  2>/dev/null
	@R -e 'library("limma")'             >> $@.tmp  2>/dev/null
	@R -e 'library("magrittr")'          >> $@.tmp  2>/dev/null
	@R -e 'library("Matrix")'            >> $@.tmp  2>/dev/null
	@R -e 'library("mclust")'            >> $@.tmp  2>/dev/null
	@R -e 'library("methods")'           >> $@.tmp  2>/dev/null
	@R -e 'library("plyr")'              >> $@.tmp  2>/dev/null
	@R -e 'library("R.utils")'           >> $@.tmp  2>/dev/null
	@R -e 'library("RCurl")'             >> $@.tmp  2>/dev/null
	@R -e 'library("readr")'             >> $@.tmp  2>/dev/null
	@R -e 'library("reshape2")'          >> $@.tmp  2>/dev/null
	@R -e 'library("rjson")'             >> $@.tmp  2>/dev/null
	@R -e 'library("Seurat")'            >> $@.tmp  2>/dev/null
	@R -e 'library("stringr")'           >> $@.tmp  2>/dev/null
	@R -e 'library("tibble")'            >> $@.tmp  2>/dev/null
	@R -e 'library("tools")'             >> $@.tmp  2>/dev/null
	@R -e 'library("tximport")'          >> $@.tmp  2>/dev/null
# R annotator script
	@R -e 'library("SRAdb")'             >> $@.tmp  2>/dev/null
	@R -e 'library("gsubfn")'            >> $@.tmp  2>/dev/null
# For topAnat
	@R -e 'library("Rgraphviz")'         >> $@.tmp  2>/dev/null
	@R -e 'library("rJava")'             >> $@.tmp  2>/dev/null
	@R -e 'library("Runiversal")'        >> $@.tmp  2>/dev/null
	@R -e 'library("topGO")'             >> $@.tmp  2>/dev/null
# Extra R libraries required during cel/MAS5 processing
	@R -e 'library("celeganscdf")'       >> $@.tmp  2>/dev/null
	@R -e 'library("celegansprobe")'     >> $@.tmp  2>/dev/null
	@R -e 'library("drosgenome1cdf")'    >> $@.tmp  2>/dev/null
	@R -e 'library("drosgenome1probe")'  >> $@.tmp  2>/dev/null
	@R -e 'library("drosophila2cdf")'    >> $@.tmp  2>/dev/null
	@R -e 'library("drosophila2probe")'  >> $@.tmp  2>/dev/null
	@R -e 'library("hgfocuscdf")'        >> $@.tmp  2>/dev/null
	@R -e 'library("hgfocusprobe")'      >> $@.tmp  2>/dev/null
	@R -e 'library("hgu95acdf")'         >> $@.tmp  2>/dev/null
	@R -e 'library("hgu95aprobe")'       >> $@.tmp  2>/dev/null
	@R -e 'library("hgu95bcdf")'         >> $@.tmp  2>/dev/null
	@R -e 'library("hgu95bprobe")'       >> $@.tmp  2>/dev/null
	@R -e 'library("hgu95ccdf")'         >> $@.tmp  2>/dev/null
	@R -e 'library("hgu95cprobe")'       >> $@.tmp  2>/dev/null
	@R -e 'library("hgu95dcdf")'         >> $@.tmp  2>/dev/null
	@R -e 'library("hgu95dprobe")'       >> $@.tmp  2>/dev/null
	@R -e 'library("hgu95ecdf")'         >> $@.tmp  2>/dev/null
	@R -e 'library("hgu95eprobe")'       >> $@.tmp  2>/dev/null
	@R -e 'library("hgu95av2cdf")'       >> $@.tmp  2>/dev/null
	@R -e 'library("hgu95av2probe")'     >> $@.tmp  2>/dev/null
	@R -e 'library("hgu133acdf")'        >> $@.tmp  2>/dev/null
	@R -e 'library("hgu133aprobe")'      >> $@.tmp  2>/dev/null
	@R -e 'library("hgu133a2cdf")'       >> $@.tmp  2>/dev/null
	@R -e 'library("hgu133a2probe")'     >> $@.tmp  2>/dev/null
	@R -e 'library("hgu133bcdf")'        >> $@.tmp  2>/dev/null
	@R -e 'library("hgu133bprobe")'      >> $@.tmp  2>/dev/null
	@R -e 'library("hgu133plus2cdf")'    >> $@.tmp  2>/dev/null
	@R -e 'library("hgu133plus2probe")'  >> $@.tmp  2>/dev/null
	@R -e 'library("htmg430acdf")'       >> $@.tmp  2>/dev/null
	@R -e 'library("htmg430aprobe")'     >> $@.tmp  2>/dev/null
	@R -e 'library("hu6800cdf")'         >> $@.tmp  2>/dev/null
	@R -e 'library("hu6800probe")'       >> $@.tmp  2>/dev/null
	@R -e 'library("mgu74acdf")'         >> $@.tmp  2>/dev/null
	@R -e 'library("mgu74aprobe")'       >> $@.tmp  2>/dev/null
	@R -e 'library("mgu74av2cdf")'       >> $@.tmp  2>/dev/null
	@R -e 'library("mgu74av2probe")'     >> $@.tmp  2>/dev/null
	@R -e 'library("mgu74bv2cdf")'       >> $@.tmp  2>/dev/null
	@R -e 'library("mgu74bv2probe")'     >> $@.tmp  2>/dev/null
	@R -e 'library("mgu74cv2cdf")'       >> $@.tmp  2>/dev/null
	@R -e 'library("mgu74cv2probe")'     >> $@.tmp  2>/dev/null
	@R -e 'library("moe430acdf")'        >> $@.tmp  2>/dev/null
	@R -e 'library("moe430aprobe")'      >> $@.tmp  2>/dev/null
	@R -e 'library("moe430bcdf")'        >> $@.tmp  2>/dev/null
	@R -e 'library("moe430bprobe")'      >> $@.tmp  2>/dev/null
	@R -e 'library("mouse4302cdf")'      >> $@.tmp  2>/dev/null
	@R -e 'library("mouse4302probe")'    >> $@.tmp  2>/dev/null
	@R -e 'library("mouse430a2cdf")'     >> $@.tmp  2>/dev/null
	@R -e 'library("mouse430a2probe")'   >> $@.tmp  2>/dev/null
	@R -e 'library("mu11ksubacdf")'      >> $@.tmp  2>/dev/null
	@R -e 'library("mu11ksubaprobe")'    >> $@.tmp  2>/dev/null
	@R -e 'library("mu11ksubbcdf")'      >> $@.tmp  2>/dev/null
	@R -e 'library("mu11ksubbprobe")'    >> $@.tmp  2>/dev/null
	@R -e 'library("rae230acdf")'        >> $@.tmp  2>/dev/null
	@R -e 'library("rae230aprobe")'      >> $@.tmp  2>/dev/null
	@R -e 'library("rgu34acdf")'         >> $@.tmp  2>/dev/null
	@R -e 'library("rgu34aprobe")'       >> $@.tmp  2>/dev/null
	@R -e 'library("rhesuscdf")'         >> $@.tmp  2>/dev/null
	@R -e 'library("rhesusprobe")'       >> $@.tmp  2>/dev/null
	@R -e 'library("zebrafishcdf")'      >> $@.tmp  2>/dev/null
	@R -e 'library("zebrafishprobe")'    >> $@.tmp  2>/dev/null
# Python libraries for Zebrafish intermine
	@python3 -c "import intermine"       >> $@.tmp  2>/dev/null
	@python3 -c "import simplejson"      >> $@.tmp  2>/dev/null
# Python libraries for download_files/analyze_all_files.ipynb
	@python3 -c "import IPython"         >> $@.tmp  2>/dev/null
	@python3 -c "import matplotlib"      >> $@.tmp  2>/dev/null
	@python3 -c "import numpy"           >> $@.tmp  2>/dev/null
	@python3 -c "import os"              >> $@.tmp  2>/dev/null
	@python3 -c "import pandas"          >> $@.tmp  2>/dev/null
	@python3 -c "import seaborn"         >> $@.tmp  2>/dev/null
# Python libraries for orthologs_paralogs_app
	@python3 -c "import SPARQLWrapper"   >> $@.tmp  2>/dev/null
	@python3 -c "import rdflib"          >> $@.tmp  2>/dev/null
	@python3 -c "import isodate"         >> $@.tmp  2>/dev/null
	@python3 -c "import pyparsing"       >> $@.tmp  2>/dev/null
	@python3 -c "import six"             >> $@.tmp  2>/dev/null
	@python3 -c "import configparser"    >> $@.tmp  2>/dev/null
	@python3 -c "import datetime"        >> $@.tmp  2>/dev/null
	@python3 -c "import distutils"       >> $@.tmp  2>/dev/null
	@python3 -c "import multiprocessing" >> $@.tmp  2>/dev/null
	@python3 -c "import pathlib"         >> $@.tmp  2>/dev/null
	@python3 -c "import typing"          >> $@.tmp  2>/dev/null
	@python3 -c "import ast"             >> $@.tmp  2>/dev/null
	@python3 -c "import errno"           >> $@.tmp  2>/dev/null
	@python3 -c "import fcntl"           >> $@.tmp  2>/dev/null
	@python3 -c "import fileinput"       >> $@.tmp  2>/dev/null
	@python3 -c "import io"              >> $@.tmp  2>/dev/null
	@python3 -c "import json"            >> $@.tmp  2>/dev/null
	@python3 -c "import logging"         >> $@.tmp  2>/dev/null
	@python3 -c "import os"              >> $@.tmp  2>/dev/null
	@python3 -c "import shutil"          >> $@.tmp  2>/dev/null
	@python3 -c "import sys"             >> $@.tmp  2>/dev/null
	@python3 -c "import getopt"          >> $@.tmp  2>/dev/null
	@python3 -c "import time"            >> $@.tmp  2>/dev/null
## Binaries
# network
	@which wget          >> $@.tmp
	@which rsync         >> $@.tmp
	@which git           >> $@.tmp
# others
	@which pipenv        >> $@.tmp
	@which mvn           >> $@.tmp
	@which gunzip        >> $@.tmp
	@which bunzip2       >> $@.tmp
	@which makeblastdb   >> $@.tmp
	@which blastn        >> $@.tmp
	@which $(CONVERTER)  >> $@.tmp
# tmpdir with enough free disk space available?
	@mkdir -p $(TMPDIR)
	@if [[ -d $(TMPDIR) ]]; then echo -n; else echo "Directory [$(TMPDIR)] does not exist"  >> $@.tmp; exit 6; fi
	@if [[ -w $(TMPDIR) ]]; then echo -n; else echo "Directory [$(TMPDIR)] is not writable" >> $@.tmp; exit 7; fi
	@echo "Think to sync the submodule(s) to get last (develop branch) annotations!"
	@$(MV) $@.tmp $@


# Build the Bgee Java application needed for some steps
$(PIPELINEROOT)$(BGEEPIPELINEJAR):
	@cd $(PIPELINEROOT)$(BGEEAPPDIR)m2_local_repository/; chmod +x mvn-install-file.sh; ./mvn-install-file.sh; $(RM) -r ~/.m2/repository/org/bbop
	@cd $(PIPELINEROOT)$(BGEEAPPDIR); mvn clean install --projects bgee-dao-api,bgee-dao-sql,bgee-pipeline; $(CP) bgee-pipeline/target/$(BGEEPIPELINEJARNAME) ../pipeline/pipeline/$(BGEEPIPELINEJAR)


## Install a converter to convert TSV files
#check_converter:
## Check if perl xls converter is there
#	@perl -MSpreadsheet::Read -e1 >$@.tmp
#	@which $(CONVERTER) >> $@.tmp
#	@$(MV) $@.tmp $@
#
#
## Export Affymetrix annotation files
#convert_annotAffy: check_converter $(PIPELINEROOT)$(AFFY_CURATION_FILEPATH)
#	# Get xls sheetnames
#	@$(CONVERTER) -i $(PIPELINEROOT)$(AFFY_CURATION_FILEPATH) >$@.tmp 2>&1
#	@echo >>$@.tmp
#	# Convert affymetrixChip       in TSV
#	@$(CONVERTER) $(CONVERTEROPT) `grep ' affymetrixChip' $@.tmp       | sed -e 's/^.* - //; s/: \[.*$$//'` $(PIPELINEROOT)$(AFFY_CURATION_FILEPATH) > $(PIPELINEROOT)$(AFFY_CHIP_FILEPATH)
#	@$(EXTRACLEANING) $(PIPELINEROOT)$(AFFY_CHIP_FILEPATH)
#	@head -5 $(PIPELINEROOT)$(AFFY_CHIP_FILEPATH)              >  $@.temp
#	@echo >>$@.temp
#	# Convert microarrayExperiment in TSV
#	@$(CONVERTER) $(CONVERTEROPT) `grep ' microarrayExperiment' $@.tmp | sed -e 's/^.* - //; s/: \[.*$$//'` $(PIPELINEROOT)$(AFFY_CURATION_FILEPATH) > $(PIPELINEROOT)$(MICROARRAY_EXPERIMENT_FILEPATH)
#	@$(EXTRACLEANING) $(PIPELINEROOT)$(MICROARRAY_EXPERIMENT_FILEPATH)
#	@head -5 $(PIPELINEROOT)$(MICROARRAY_EXPERIMENT_FILEPATH)  >> $@.temp
#	@echo >>$@.temp
#	# Convert chipType             in TSV
#	@$(CONVERTER) $(CONVERTEROPT) `grep ' chipType' $@.tmp             | sed -e 's/^.* - //; s/: \[.*$$//'` $(PIPELINEROOT)$(AFFY_CURATION_FILEPATH) > $(PIPELINEROOT)$(AFFY_CHIPTYPE_FILEPATH)
#	@$(EXTRACLEANING) $(PIPELINEROOT)$(AFFY_CHIPTYPE_FILEPATH)
#	@head -5 $(PIPELINEROOT)$(AFFY_CHIPTYPE_FILEPATH)          >> $@.temp
#	@$(CAT) $@.tmp $@.temp > $@
#	@$(RM) $@.tmp $@.temp
#
## Export RNA_seq annotation files
#convert_annotRNASeq: check_converter $(PIPELINEROOT)$(RNASEQ_CURATION_FILEPATH)
#	# Get xls sheetnames
#	@$(CONVERTER) -i $(PIPELINEROOT)$(RNASEQ_CURATION_FILEPATH) >$@.tmp 2>&1
#	@echo >>$@.tmp
#	# Convert RNAseqExperiment in TSV
#	@$(CONVERTER) $(CONVERTEROPT) `grep ' RNAseqExperiment' $@.tmp | sed -e 's/^.* - //; s/: \[.*$$//'` $(PIPELINEROOT)$(RNASEQ_CURATION_FILEPATH) > $(PIPELINEROOT)$(RNASEQ_EXPERIMENT_FILEPATH)
#	@$(EXTRACLEANING) $(PIPELINEROOT)$(RNASEQ_EXPERIMENT_FILEPATH)
#	@head -5 $(PIPELINEROOT)$(RNASEQ_EXPERIMENT_FILEPATH) >  $@.temp
#	@echo >>$@.temp
#	# Convert RNAseqLibrary    in TSV
#	@$(CONVERTER) $(CONVERTEROPT) `grep ' RNAseqLibrary'    $@.tmp | sed -e 's/^.* - //; s/: \[.*$$//'` $(PIPELINEROOT)$(RNASEQ_CURATION_FILEPATH) > $(PIPELINEROOT)$(RNASEQ_LIB_FILEPATH)
#	@$(EXTRACLEANING) $(PIPELINEROOT)$(RNASEQ_LIB_FILEPATH)
#	@head -5 $(PIPELINEROOT)$(RNASEQ_LIB_FILEPATH)       >> $@.temp
#	@$(CAT) $@.tmp $@.temp > $@
#	@$(RM) $@.tmp $@.temp

#FIXME used?
# Export homology annotation file
$(PIPELINEROOT)$(SIMILARITYFILEPATH): $(CONVERTER) $(PIPELINEROOT)$(SIMILARITY_CURATION_FILEPATH)
	@$(CONVERTER) $(CONVERTEROPT) 01 $(PIPELINEROOT)$(SIMILARITY_CURATION_FILEPATH) > $@
	@$(EXTRACLEANING) $@

build_container:
	@which singularity
	sudo singularity build ${PIPELINE_CONTAINER_FILE} ${PIPELINE_CONTAINER_NAME}.def
	@touch $@

test_container: build_container
# clients/interpretors
	@$(PIPELINE_CONTAINER_CMD) which java           > $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which perl          >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which mysql         >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which R             >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which python3       >> $@.tmp
# Ensembl API
	@$(PIPELINE_CONTAINER_CMD) perl -MBio::EnsEMBL::Registry  -e 1 >> $@.tmp
# other Perl modules
	@$(PIPELINE_CONTAINER_CMD) perl -MArchive::Tar            -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MArchive::Zip            -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MArray::Utils            -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MBio::SeqIO              -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MCpanel::JSON::XS        -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MData::Dumper            -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MDBI                     -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MDigest::SHA             -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MExporter                -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MFile::Basename          -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MFile::Copy              -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MFile::Find              -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MFile::Path              -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MFile::Slurp             -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MFile::Spec              -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MFindBin                 -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MGetopt::Long            -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MHTTP::Request           -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MIO::Compress::Gzip      -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MIO::Socket              -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MIO::Uncompress::Gunzip  -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MLWP::Simple             -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MLWP::UserAgent          -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MList::Compare           -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MList::MoreUtils         -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MList::Util              -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MParallel::ForkManager   -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MParallel::Loops         -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MSort::Naturally         -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MSpreadsheet::Read       -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MSpreadsheet::ParseExcel -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MText::CSV               -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MSpreadsheet::ParseXLSX  -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MSpreadsheet::WriteExcel -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MTest::More              -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MTime::HiRes             -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MTime::localtime         -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MWebservice::InterMine   -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MXML::Fast               -e 1 >> $@.tmp
# R libraries
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("affy")'              >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("AnnotationDbi")'     >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("BgeeCall")'          >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("Biobase")'           >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("biomaRt")'           >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("Biostrings")'        >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("BUSpaRse")'          >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("data.table")'        >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("devtools")'          >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("dplyr")'             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("DropletUtils")'      >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("edgeR")'             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("fdrtool")'           >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("forcats")'           >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("gcrma")'             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("genefilter")'        >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("GenomicFeatures")'   >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("ggExtra")'           >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("gghighlight")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("ggplot2")'           >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("gridExtra")'         >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("HCAExplorer")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("HelpersMG")'         >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("LaplacesDemon")'     >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("lattice")'           >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("limma")'             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("magrittr")'          >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("Matrix")'            >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("mclust")'            >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("methods")'           >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("plyr")'              >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("R.utils")'           >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("RCurl")'             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("readr")'             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("reshape2")'          >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("rjson")'             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("Seurat")'            >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("stringr")'           >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("tibble")'            >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("tools")'             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("tximport")'          >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("UniprotR")'          >> $@.tmp  2>/dev/null
# R annotator script
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("SRAdb")'             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("gsubfn")'            >> $@.tmp  2>/dev/null
# For topAnat
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("Rgraphviz")'         >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("rJava")'             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("Runiversal")'        >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("topGO")'             >> $@.tmp  2>/dev/null
# Extra R libraries required during cel/MAS5 processing
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("celeganscdf")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("celegansprobe")'     >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("drosgenome1cdf")'    >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("drosgenome1probe")'  >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("drosophila2cdf")'    >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("drosophila2probe")'  >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgfocuscdf")'        >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgfocusprobe")'      >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu95acdf")'         >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu95aprobe")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu95bcdf")'         >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu95bprobe")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu95ccdf")'         >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu95cprobe")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu95dcdf")'         >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu95dprobe")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu95ecdf")'         >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu95eprobe")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu95av2cdf")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu95av2probe")'     >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu133acdf")'        >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu133aprobe")'      >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu133a2cdf")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu133a2probe")'     >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu133bcdf")'        >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu133bprobe")'      >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu133plus2cdf")'    >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hgu133plus2probe")'  >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("htmg430acdf")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("htmg430aprobe")'     >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hu6800cdf")'         >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("hu6800probe")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("mgu74acdf")'         >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("mgu74aprobe")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("mgu74av2cdf")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("mgu74av2probe")'     >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("mgu74bv2cdf")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("mgu74bv2probe")'     >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("mgu74cv2cdf")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("mgu74cv2probe")'     >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("moe430acdf")'        >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("moe430aprobe")'      >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("moe430bcdf")'        >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("moe430bprobe")'      >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("mouse4302cdf")'      >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("mouse4302probe")'    >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("mouse430a2cdf")'     >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("mouse430a2probe")'   >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("mu11ksubacdf")'      >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("mu11ksubaprobe")'    >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("mu11ksubbcdf")'      >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("mu11ksubbprobe")'    >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("rae230acdf")'        >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("rae230aprobe")'      >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("rgu34acdf")'         >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("rgu34aprobe")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("rhesuscdf")'         >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("rhesusprobe")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("zebrafishcdf")'      >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -e 'library("zebrafishprobe")'    >> $@.tmp  2>/dev/null
# Python libraries for Zebrafish intermine
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import intermine"       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import simplejson"      >> $@.tmp  2>/dev/null
# Python libraries for download_files/analyze_all_files.ipynb
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import IPython"         >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import matplotlib"      >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import numpy"           >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import os"              >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import pandas"          >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import seaborn"         >> $@.tmp  2>/dev/null
# Python libraries for Sagane's annotation scripts
# TODO ...
# Python libraries for Sagane's H5AD scripts
# TODO ...
# Python libraries for orthologs_paralogs_app
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import SPARQLWrapper"   >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import rdflib"          >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import isodate"         >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import pyparsing"       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import six"             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import configparser"    >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import datetime"        >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import distutils"       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import multiprocessing" >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import pathlib"         >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import typing"          >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import ast"             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import errno"           >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import fcntl"           >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import fileinput"       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import io"              >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import json"            >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import logging"         >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import os"              >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import shutil"          >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import sys"             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import getopt"          >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import time"            >> $@.tmp  2>/dev/null
# others
	@$(PIPELINE_CONTAINER_CMD) which pipenv        >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which mvn           >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which gunzip        >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which bunzip2       >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which makeblastdb   >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which blastn        >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which kallisto      >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which bustools      >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which fastp         >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which fasterq-dump  >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which gtf_to_fasta  >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which $(CONVERTER)  >> $@.tmp
	@$(MV) $@.tmp $@

$(VERIFICATIONFILE): global_check $(PIPELINEROOT)$(SIMILARITYFILEPATH) $(PIPELINEROOT)$(BGEEPIPELINEJAR) convert_annotAffy convert_annotRNASeq
# create empty file or clean existing file
	@touch $@.tmp
	@> $@.tmp
# verifications for homology annotation file
	@echo ---------- $(SIMILARITYFILE)--------- >> $@.tmp
	@{ echo "Line count: \c" && wc -l < $(PIPELINEROOT)$(SIMILARITYFILEPATH); } >> $@.tmp
	@echo 5 first lines: >> $@.tmp
	@head -5 $(PIPELINEROOT)$(SIMILARITYFILEPATH) >> $@.tmp
# verifications for Affymetrix annotation file
	@echo ---------- AFFYMETRIX ANNOTATIONS --------- >> $@.tmp
	@{ echo "Line count in microarray experiments: \c" && wc -l < $(PIPELINEROOT)$(MICROARRAY_EXPERIMENT_FILEPATH); } >> $@.tmp
	@{ echo "Line count in affymetrix chips: \c" && wc -l < $(PIPELINEROOT)$(AFFY_CHIP_FILEPATH); } >> $@.tmp
	@{ echo "Line count in affymetrix chip types: \c" && wc -l < $(PIPELINEROOT)$(AFFY_CHIPTYPE_FILEPATH); } >> $@.tmp
	@echo 5 first lines of each file: >> $@.tmp
	@$(CAT) convert_annotAffy >> $@.tmp
	@$(MV) $@.tmp $@

clean:
	-@$(RM) global_check $(CONVERTER) $(PIPELINEROOT)$(SIMILARITYFILEPATH) $(PIPELINEROOT)$(BGEEPIPELINEJAR) $(VERIFICATIONFILE) convert_annotAffy convert_annotRNASeq check_converter *.tmp *.temp

