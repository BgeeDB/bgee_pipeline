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
# other Perl modules
	@perl -MArchive::Tar            -e 1 >> $@.tmp
	@perl -MArchive::Zip            -e 1 >> $@.tmp
	@perl -MArray::Utils            -e 1 >> $@.tmp
	@perl -MBio::SeqIO              -e 1 >> $@.tmp
	@perl -MCpanel::JSON::XS        -e 1 >> $@.tmp
	@perl -MData::Dumper            -e 1 >> $@.tmp
	@perl -MData::Dump              -e 1 >> $@.tmp
	@perl -MDBI                     -e 1 >> $@.tmp
	@perl -MDBD::mysql              -e 1 >> $@.tmp
	@perl -MDBD::Pg                 -e 1 >> $@.tmp
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
	@perl -MIO::Compress::Gzip      -e 1 >> $@.tmp
	@perl -MIO::Socket              -e 1 >> $@.tmp
	@perl -MIO::Uncompress::Gunzip  -e 1 >> $@.tmp
	@perl -MJSON                    -e 1 >> $@.tmp
	@perl -MJSON::XS                -e 1 >> $@.tmp
	@perl -MLWP::UserAgent          -e 1 >> $@.tmp
	@perl -MLWP::Simple             -e 1 >> $@.tmp
	@perl -MList::Compare           -e 1 >> $@.tmp
	@perl -MList::MoreUtils         -e 1 >> $@.tmp
	@perl -MList::Util              -e 1 >> $@.tmp
	@perl -MParallel::ForkManager   -e 1 >> $@.tmp
	@perl -MSort::Naturally         -e 1 >> $@.tmp
	@perl -MSpreadsheet::Read       -e 1 >> $@.tmp
	@perl -MSpreadsheet::ParseExcel -e 1 >> $@.tmp
	@perl -MText::CSV               -e 1 >> $@.tmp
	@perl -MSpreadsheet::ParseXLSX  -e 1 >> $@.tmp
	@perl -MSpreadsheet::WriteExcel -e 1 >> $@.tmp
	@perl -MTest::More              -e 1 >> $@.tmp
	@perl -MTie::IxHash             -e 1 >> $@.tmp
	@perl -MTime::HiRes             -e 1 >> $@.tmp
	@perl -MTime::localtime         -e 1 >> $@.tmp
	@perl -MTime::Piece             -e 1 >> $@.tmp
	@perl -MTry::Tiny               -e 1 >> $@.tmp
	@perl -MWebservice::InterMine   -e 1 >> $@.tmp
	@perl -MXML::Fast               -e 1 >> $@.tmp
# R libraries
	@R -e 'library("BgeeCall")'          >> $@.tmp  2>/dev/null
	@R -e 'library("Biostrings")'        >> $@.tmp  2>/dev/null
	@R -e 'library("BUSpaRse")'          >> $@.tmp  2>/dev/null
	@R -e 'library("data.table")'        >> $@.tmp  2>/dev/null
	@R -e 'library("dplyr")'             >> $@.tmp  2>/dev/null
	@R -e 'library("doParallel")'        >> $@.tmp  2>/dev/null
	@R -e 'library("DropletUtils")'      >> $@.tmp  2>/dev/null
	@R -e 'library("edgeR")'             >> $@.tmp  2>/dev/null
	@R -e 'library("foreach")'           >> $@.tmp  2>/dev/null
	@R -e 'library("GenomicFeatures")'   >> $@.tmp  2>/dev/null
	@R -e 'library("ggExtra")'           >> $@.tmp  2>/dev/null
	@R -e 'library("gghighlight")'       >> $@.tmp  2>/dev/null
	@R -e 'library("ggplot2")'           >> $@.tmp  2>/dev/null
	@R -e 'library("gridExtra")'         >> $@.tmp  2>/dev/null
	@R -e 'library("HelpersMG")'         >> $@.tmp  2>/dev/null
	@R -e 'library("LaplacesDemon")'     >> $@.tmp  2>/dev/null
	@R -e 'library("limma")'             >> $@.tmp  2>/dev/null
	@R -e 'library("Matrix")'            >> $@.tmp  2>/dev/null
	@R -e 'library("mclust")'            >> $@.tmp  2>/dev/null
	@R -e 'library("plyr")'              >> $@.tmp  2>/dev/null
	@R -e 'library("R.utils")'           >> $@.tmp  2>/dev/null
	@R -e 'library("RCurl")'             >> $@.tmp  2>/dev/null
	@R -e 'library("readr")'             >> $@.tmp  2>/dev/null
	@R -e 'library("reshape2")'          >> $@.tmp  2>/dev/null
	@R -e 'library("rjson")'             >> $@.tmp  2>/dev/null
	@R -e 'library("RMySQL")'            >> $@.tmp  2>/dev/null
	@R -e 'library("rtracklayer")'       >> $@.tmp  2>/dev/null
	@R -e 'library("Seurat")'            >> $@.tmp  2>/dev/null
	@R -e 'library("stringr")'           >> $@.tmp  2>/dev/null
	@R -e 'library("tools")'             >> $@.tmp  2>/dev/null
	@R -e 'library("tximport")'          >> $@.tmp  2>/dev/null
	@R -e 'library("UniProt.ws")'        >> $@.tmp  2>/dev/null
# R annotator script
	@R -e 'library("SRAdb")'             >> $@.tmp  2>/dev/null
	@R -e 'library("gsubfn")'            >> $@.tmp  2>/dev/null
# For topAnat
	@R -e 'library("Rgraphviz")'         >> $@.tmp  2>/dev/null
	@R -e 'library("rJava")'             >> $@.tmp  2>/dev/null
	@R -e 'library("Runiversal")'        >> $@.tmp  2>/dev/null
	@R -e 'library("topGO")'             >> $@.tmp  2>/dev/null
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
#	@python3 -c "import distutils"       >> $@.tmp  2>/dev/null #FIXME no more in Python 3.12+
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
# Python libraries for H5AD script
	@python3 -c "import pandas"          >> $@.tmp  2>/dev/null
	@python3 -c "import mysql.connector" >> $@.tmp  2>/dev/null
	@python3 -c "import numpy"           >> $@.tmp  2>/dev/null
	@python3 -c "import scipy"           >> $@.tmp  2>/dev/null
	@python3 -c "import anndata"         >> $@.tmp  2>/dev/null
	@python3 -c "import tqdm"            >> $@.tmp  2>/dev/null
	@python3 -c "import argparse"        >> $@.tmp  2>/dev/null
# Sagane's scripts
	@python3 -c "import collections"     >> $@.tmp  2>/dev/null
	@python3 -c "import concurrent"      >> $@.tmp  2>/dev/null
	@python3 -c "import csv"             >> $@.tmp  2>/dev/null
	@python3 -c "import difflib"         >> $@.tmp  2>/dev/null
	@python3 -c "import GEOparse"        >> $@.tmp  2>/dev/null
	@python3 -c "import h5py"            >> $@.tmp  2>/dev/null
	@python3 -c "import math"            >> $@.tmp  2>/dev/null
	@python3 -c "import networkx"        >> $@.tmp  2>/dev/null
	@python3 -c "import pandas"          >> $@.tmp  2>/dev/null
	@python3 -c "import pyvis"           >> $@.tmp  2>/dev/null
	@python3 -c "import random"          >> $@.tmp  2>/dev/null
	@python3 -c "import re"              >> $@.tmp  2>/dev/null
	@python3 -c "import requests"        >> $@.tmp  2>/dev/null
	@python3 -c "import scanpy"          >> $@.tmp  2>/dev/null
	@python3 -c "import statistics"      >> $@.tmp  2>/dev/null
	@python3 -c "import subprocess"      >> $@.tmp  2>/dev/null
	@python3 -c "import tempfile"        >> $@.tmp  2>/dev/null
	@python3 -c "import xml"             >> $@.tmp  2>/dev/null
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
	@which unxz          >> $@.tmp
	@which makeblastdb   >> $@.tmp
	@which blastn        >> $@.tmp
	@which kallisto      >> $@.tmp
	@which bustools      >> $@.tmp
	@which fastp         >> $@.tmp
	@which fasterq-dump  >> $@.tmp
	@which gtf_to_fasta  >> $@.tmp
	@which bamtofastq    >> $@.tmp
	@which esearch       >> $@.tmp
	@which xtract        >> $@.tmp
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

build_container: ${PIPELINE_CONTAINER_NAME}.def
	@which apptainer
	apptainer build ${PIPELINE_CONTAINER_FILE} ${PIPELINE_CONTAINER_NAME}.def 2>&1 > $@.tmp
	@$(MV) $@.tmp $@

test_container: build_container
# clients/interpretors
	@$(PIPELINE_CONTAINER_CMD) which java           > $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which perl          >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which mysql         >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which R             >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which python3       >> $@.tmp
# other Perl modules
	@$(PIPELINE_CONTAINER_CMD) perl -MArchive::Tar            -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MArchive::Zip            -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MArray::Utils            -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MBio::SeqIO              -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MCpanel::JSON::XS        -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MData::Dumper            -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MData::Dump              -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MDBI                     -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MDBD::mysql              -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MDBD::Pg                 -e 1 >> $@.tmp
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
	@$(PIPELINE_CONTAINER_CMD) perl -MIO::Compress::Gzip      -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MIO::Socket              -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MIO::Uncompress::Gunzip  -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MJSON                    -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MJSON::XS                -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MLWP::Simple             -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MLWP::UserAgent          -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MList::Compare           -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MList::MoreUtils         -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MList::Util              -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MParallel::ForkManager   -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MSort::Naturally         -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MSpreadsheet::Read       -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MSpreadsheet::ParseExcel -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MText::CSV               -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MSpreadsheet::ParseXLSX  -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MSpreadsheet::WriteExcel -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MTest::More              -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MTie::IxHash             -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MTime::HiRes             -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MTime::localtime         -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MTime::Piece             -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MTry::Tiny               -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MWebservice::InterMine   -e 1 >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) perl -MXML::Fast               -e 1 >> $@.tmp
# R libraries
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("BgeeCall")'          >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("Biostrings")'        >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("BUSpaRse")'          >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("data.table")'        >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("dplyr")'             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("doParallel")'        >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("DropletUtils")'      >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("edgeR")'             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("foreach")'           >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("GenomicFeatures")'   >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("ggExtra")'           >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("gghighlight")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("ggplot2")'           >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("gridExtra")'         >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("HelpersMG")'         >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("LaplacesDemon")'     >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("limma")'             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("Matrix")'            >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("mclust")'            >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("plyr")'              >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("R.utils")'           >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("RCurl")'             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("readr")'             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("reshape2")'          >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("rjson")'             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("RMySQL")'            >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("rtracklayer")'       >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("Seurat")'            >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("stringr")'           >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("tools")'             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("tximport")'          >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("UniProt.ws")'        >> $@.tmp  2>/dev/null
# R annotator script
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("SRAdb")'             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("gsubfn")'            >> $@.tmp  2>/dev/null
# For topAnat
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("Rgraphviz")'         >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("rJava")'             >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("Runiversal")'        >> $@.tmp  2>/dev/null
	@$(PIPELINE_CONTAINER_CMD) R -q -e 'library("topGO")'             >> $@.tmp  2>/dev/null
# Python libraries for Zebrafish intermine
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import intermine"       >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import simplejson"      >> $@.tmp
# Python libraries for download_files/analyze_all_files.ipynb
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import IPython"         >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import matplotlib"      >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import numpy"           >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import os"              >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import pandas"          >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import seaborn"         >> $@.tmp
# Python libraries for Sagane's annotation scripts
# TODO ...
# Python libraries for Sagane's H5AD scripts
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import pandas"          >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import mysql.connector" >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import numpy"           >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import scipy"           >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import anndata"         >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import tqdm"            >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import argparse"        >> $@.tmp
# Python libraries for orthologs_paralogs_app
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import SPARQLWrapper"   >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import rdflib"          >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import isodate"         >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import pyparsing"       >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import six"             >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import configparser"    >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import datetime"        >> $@.tmp
#@$(PIPELINE_CONTAINER_CMD) python3 -c "import distutils"       >> $@.tmp #FIXME no more in Python 3.12+
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import multiprocessing" >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import pathlib"         >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import typing"          >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import ast"             >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import errno"           >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import fcntl"           >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import fileinput"       >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import io"              >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import json"            >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import logging"         >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import os"              >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import shutil"          >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import sys"             >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import getopt"          >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import time"            >> $@.tmp
# Sagane's scripts
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import collections"     >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import concurrent"      >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import csv"             >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import difflib"         >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import GEOparse"        >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import h5py"            >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import math"            >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import networkx"        >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import pandas"          >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import pyvis"           >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import random"          >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import re"              >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import requests"        >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import scanpy"          >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import statistics"      >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import subprocess"      >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import tempfile"        >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) python3 -c "import xml"             >> $@.tmp
# others
	@$(PIPELINE_CONTAINER_CMD) which pipenv        >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which mvn           >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which gunzip        >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which bunzip2       >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which unxz          >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which makeblastdb   >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which blastn        >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which kallisto      >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which bustools      >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which fastp         >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which fasterq-dump  >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which gtf_to_fasta  >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which bamtofastq    >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which esearch       >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which xtract        >> $@.tmp
	@$(PIPELINE_CONTAINER_CMD) which $(CONVERTER)  >> $@.tmp
	@$(MV) $@.tmp $@

$(VERIFICATIONFILE): global_check $(PIPELINEROOT)$(SIMILARITYFILEPATH) $(PIPELINEROOT)$(BGEEPIPELINEJAR) convert_annotRNASeq
# create empty file or clean existing file
	@touch $@.tmp
	@> $@.tmp
# verifications for homology annotation file
	@echo ---------- $(SIMILARITYFILE)--------- >> $@.tmp
	@{ echo "Line count: \c" && wc -l < $(PIPELINEROOT)$(SIMILARITYFILEPATH); } >> $@.tmp
	@echo 5 first lines: >> $@.tmp
	@head -5 $(PIPELINEROOT)$(SIMILARITYFILEPATH) >> $@.tmp
	@$(MV) $@.tmp $@

clean:
	-@$(RM) global_check $(CONVERTER) $(PIPELINEROOT)$(SIMILARITYFILEPATH) $(PIPELINEROOT)$(BGEEPIPELINEJAR) $(VERIFICATIONFILE) convert_annotRNASeq check_converter *.tmp *.temp

