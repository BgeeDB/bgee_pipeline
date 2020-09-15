#!/usr/bin/env perl

# Perl embedded modules
use strict;
use warnings;
use diagnostics;

use File::Slurp;
use lib '.';
use Utils;

my $in_anat    = $ARGV[0]  or die "\n\t$0 anatId_list anat_port stageId_list stage_port\n\n";
my $anat_port  = $ARGV[1]  or die "\n\t$0 anatId_list anat_port stageId_list stage_port\n\n";
#TODO later
#my $in_stage   = $ARGV[2]  or die "\n\t$0 anatId_list anat_port stageId_list stage_port\n\n";
#my $stage_port = $ARGV[3]  or die "\n\t$0 anatId_list anat_port stageId_list stage_port\n\n";


my @anatId  = read_file("$in_anat",  chomp=>1);
#TODO later
#my @stageId = read_file("$in_stage", chomp=>1);

#print "[$anatId[100]]\t[$stageId[100]]\n";

# Launch the organ stage mapping tool (using $(CUSTOM_UBERON_FILE_PATH) and $(DEV_STAGE_ONT_FILE_PATH))
#@$(IDMAPPING)  $(IDMAPPINGPORT) &
#@$(STGMAPPING) $(STGMAPPINGPORT) &
#java -Xmx32g -Dbgee.dao.jdbc.username=$(DBUSER) -Dbgee.dao.jdbc.password=$(DBPASS) -Dbgee.dao.jdbc.driver.names=com.mysql.jdbc.Driver,net.sf.log4jdbc.sql.jdbcapi.DriverSpy -Dbgee.dao.jdbc.url='jdbc:log4jdbc:mysql://$(DBHOST):$(DBPORT)/$(DBNAME)?enableQueryTimeouts=false&sessionVariables=net_write_timeout=86400,net_read_timeout=86400,wait_timeout=86400' -jar ../java/bgee-pipeline-15-with-dependencies.jar UberonSocketTool idMapping custom_composite.owl 14555 &
#
#nohup java -Xmx32g -jar ../java/bgee-pipeline-15-with-dependencies.jar  UberonSocketTool idMapping  ../source_files/uberon/composite-metazoan.owl  14555 &
#nohup java -Xmx32g -jar ../java/bgee-pipeline-15-with-dependencies.jar  UberonSocketTool idMapping  ../source_files/uberon/dev_stage_ontology.obo  13222 &
#@sleep 50 # sleep because mappers need time to load Uberon


#get_anatomy_mapping
my $doneAnat = Utils::get_anatomy_mapping(\@anatId,  $anat_port,  0);
#TODO later
#my $doneStg  = Utils::get_anatomy_mapping(\@stageId, $stage_port, 0);


exit 0;

