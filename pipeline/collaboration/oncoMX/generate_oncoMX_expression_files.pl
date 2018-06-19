#!/usr/bin/env perl

# Julien Wollbrett, September 2017

# This script is reponsible for generating the oncoMX expression files
# containing processed expression data (not the calls,
# but the TPM values).
# See opt debug message below for information about parameters.
# This script is based on the generate_ref_expr_files.pl script
# written by Frederic Bastian

#############################################################

use strict;
use warnings;
use diagnostics;
use Data::Dumper;

use Archive::Zip qw( :ERROR_CODES :CONSTANTS );
use Getopt::Long;
use File::Spec;
use File::Path qw(make_path remove_tree);
use File::Copy qw(move);
use List::Util qw(sum);
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;
use IO::Compress::Gzip qw(gzip $GzipError) ;

# forces a flush after every write or print, so the output appears as soon as it's generated rather than being buffered.
$| = 1;

# Define arguments & their default value
my ($bgee_connector)                = ('');
my ($speciesArg)                    = ('');
my ($devStageArg)                   = ('');
my ($outputDir)                     = ('');
my ($bgeeVersion)                   = ('');
my ($debug)                         = (0);
my %opts = ('bgee=s'               => \$bgee_connector,     # Bgee connector string
            'speciesArg=s'         => \$speciesArg,
            'devStageArg=s'        => \$devStageArg,
            'outputDir=s'          => \$outputDir,
            'bgeeVersion=s'        => \$bgeeVersion,
            'debug'                => \$debug
           );
######################################         
# 		Check arguments
######################################         

my $emptyArg = '-';
my $outputFileNamePrefix = "oncoMX_expression_";
my $outputFileNameSuffix = ".tsv";
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $speciesArg eq '' || $devStageArg eq '' || $outputDir eq '' || $bgeeVersion eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEE_CMD) -speciesArg=\$(SPECIES_ARG) -affyDir=\$(AFFY_DIR) -estDir=\$(EST_DIR) -inSituDir=\$(IN_SITU_DIR) -rnaSeqDir=\$(RNA_SEQ_DIR)
\t-bgee                 Bgee    connector string
\t-speciesArg           A comma-separated list of species IDs to generate files for, or '-' to generate files for all species
\t-devStageArg          Parent dev. stage of all stages that will be taken into account to create the file, or '-' to generate files without dev stage restriction
\t-outputDir            Path to the directory where oncoMX files will be generated (one file per species)
\t-bgeeVersion          The Bgee release for which the files are being generated, e.g. 'bgee_v14'.
\t-debug                more verbose output
\n";
    exit 1;
}

######################################         
# Retrieve information needed for all species
######################################         



# -----------------------
# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);
# species list
my @speciesList = ();
if ($speciesArg ne $emptyArg) {
    @speciesList = split(',', $speciesArg);
}


my %species;

my $speciesSql = 'SELECT speciesId, genus, species FROM species';
if (@speciesList) {
    $speciesSql .= ' WHERE speciesId IN (';
    for my $i (0 .. $#speciesList) {
        if ($i > 0) {
            $speciesSql .= ', ';
        }
        $speciesSql .= $speciesList[$i];
    }
    $speciesSql .= ')';
}
my $speciesStmt = $dbh->prepare($speciesSql);
$speciesStmt->execute()  or die $speciesStmt->errstr;
while ( my @data = $speciesStmt->fetchrow_array ){
    $species{$data[0]}{'name'} = $data[1].' '.$data[2];
}


######################################         
#for each species, generate the files
######################################         

for my $speId ( keys %species ){
    print "Generating file for species $speId...\n";
    if ($outputDir ne $emptyArg) {
        my $absDir = $outputDir;
        if ( !File::Spec->file_name_is_absolute($outputDir) ){
            $absDir = File::Spec->rel2abs($outputDir) ;
        }
#         my $maxTPMValue = computeMaxMeanTPM($speId);
        generateFile($speId, $species{$speId}{'name'}, $absDir);
    }else{
		print "path to output directory can't be empty";
    }
    print "Done generating files for species $speId.\n";
}

$dbh->disconnect;
exit 0;

######################################         
# SUBS
######################################         

sub generateFile {
    my @args = @_;
    my $speciesId      = $args[0];
    my $speciesName    = $args[1];
    my $fileDir       = $args[2];

    my $speciesNameForFile = $speciesName;
    $speciesNameForFile =~ s/ /_/g;
    my $tmpDirName = $speciesNameForFile.'_tmp';
    my $expFileName = $speciesNameForFile.'_ontoMX.tsv';
    # recreate directory for experiment and library information
    remove_tree(File::Spec->catdir( $fileDir, $tmpDirName));
    make_path(File::Spec->catdir( $fileDir, $tmpDirName));

    # path to processed expression data on FTP
    my $exprFilePath = $bgeeVersion.'/download/ontoMXFiles/'
        .$speciesNameForFile.'/';

######################################         
#		Retrieve all genes used in RNASeq experiments and having one Uniprot/SwissProt 
######################################         

    # ID in Bgee
    # $genes{bgeeGeneId}                = bgee gene ID
    # $genes{bgeeGeneId}{'uniProtId '}  = uniprot ID
    my %genes = retrieveSwissprotGenes($speciesId);
    
######################################         
# 		Now we retrieve all children stage IDs of the dev. stage Id $devStageArg
######################################         
    
    my %stages = retrieveChildrenDevStages($devStageArg, $speciesId);
	
######################################         
# 		Retrieve expression information for each gene in all conditions
######################################         

    # $tpms{anatEntityId}                     			 = anatomical entity ID
    # $tpms{anatEntityId}{'maxOfMeanTPM'}           	 = mean TPM value calculate using mean o
    # $tpms{anatEntityId}{stageId}             		     = developmental stage IDs
    # $tpms{anatEntityId}{stageId}{sex}                  = sex (Male or Female)
    # @tpms{anatEntityId}{stageId}{sex}{presence}        = List of all TPMs for present/absent genes in this condition
    #
    #
    #
    # $expr{'uniprotId'}                                 = uniprot ID
	# $expr{'anatEntityId'}                              = anatomical entity ID
	# $expr{'stageId'}                                   = developmental stage ID
	# $expr{'sex'}                                       = sex (Male or Female)
	# $expr{'presence'}                                  = presence or absence of expression
	# $expr{'tpm'}                                       = avg expression value for this gene in these conditions
	#
	# $tau{uniprotId}                                    = tau value for one uniprot/swissprot ID
	#
	#
	
	# each entry of @expr will corresponds to one line in one output file
	my @expr = ();
	# hash with uniprot/swissprot ID as key and corresponding tau as value
	my %tau;
	#max mean log2 TPM for all genes
	my $maxLog2TPM = 0;
	
    print 'Number of genes : '.(scalar keys %genes)."\n";
	print "Start retrieving TPM values and calculation of Tau for each gene of species $speciesId\n";
	for my $bgeeGeneId ( keys %genes ){
		# hash allowing to easily filter by condition of expression
		my %tpms;
		# max log2 TPM values for all anatomical entities
		my $maxXiTau = 0;
				
        my $sqlTPMs = 'SELECT c.anatEntityId, c.stageId, c.sex, res.detectionFlag, res.tpm '. 
        'FROM expression AS e INNER JOIN rnaSeqResult AS res ON e.expressionId = res.expressionId '.
        'INNER JOIN cond AS c ON c.conditionId = e.conditionId '.
        'WHERE e.bgeeGeneId = ? order by anatEntityId, stageId, sex;';

    	my $stmtTPMs = $dbh->prepare($sqlTPMs);
    	$stmtTPMs->execute($bgeeGeneId)  or die $stmtTPMs->errstr;
    	while ( my @data = $stmtTPMs->fetchrow_array ){
			push @{ $tpms{$data[0]}{$data[1]}{$data[2]}{$data[3]} }, $data[4];
        }
         # parse all entries of this hash of conditions in order to populate @expr
        for my $anatEntityId (keys %tpms){
        	$tpms{$anatEntityId}{'maxOfMeanTPM'} = 0;
        	for my $stageId (keys % {$tpms{$anatEntityId}}){
        		#test presence of the current anat. entity ID in the hash of allowed anat entities %stages
        		if($devStageArg eq $emptyArg || exists($stages{$stageId}) || ( !%stages )) {
        		    for my $sex (keys %{$tpms{$anatEntityId}{$stageId}}){
        		    	my %expr;
        		    	my $log2MeanTPM = 0;
        		    	$expr{'uniprotId'} = $genes{$bgeeGeneId};
        		    	$expr{'anatEntityId'} = $anatEntityId;
        		    	$expr{'stageId'} = $stageId;
        		    	$expr{'sex'} = $sex;
        		    	my $presence ;
        			    if ( exists($tpms{$anatEntityId}{$stageId}{$sex}{'present'})){
        			    	$presence = 'present';
        			    }elsif (exists($tpms{$anatEntityId}{$stageId}{$sex}{'absent'})){
        			    	$presence = 'absent';
        			    }else{
        			    	exit("no Presence/Absence expression data");
        			    }
						my @tpmValues = @{$tpms{$anatEntityId}{$stageId}{$sex}{$presence}};
						my $log2Sum = 0;
						for my $tpmValue (@tpmValues){
							#log2 transformation create negatives values if initial value is < 1.
							#That's why we add 1 at all tpm values in order to avoid values between 0 and 1
        			    	$log2Sum += log($tpmValue+1)/log(2);
        			    }
        				$log2MeanTPM = $log2Sum/($#tpmValues+1);
        				$expr{'presence'} = $presence;
        				$expr{'tpm'} = $log2MeanTPM;
        				push @expr, \%expr;
        			    #test if this mean TPM value is bigger than max mean TPM value for current anatomical entity
        			    if($log2MeanTPM > $tpms{$anatEntityId}{'maxOfMeanTPM'}){
        			    	$tpms{$anatEntityId}{'maxOfMeanTPM'} = $log2MeanTPM;
        			    }
        			}
        		}
        	}
        	if($tpms{$anatEntityId}{'maxOfMeanTPM'} > $maxXiTau){
       		   	$maxXiTau = $tpms{$anatEntityId}{'maxOfMeanTPM'};
        	}
        }
        if($maxLog2TPM < $maxXiTau){
        	 $maxLog2TPM = $maxXiTau
        }
######################################         
#       Tau calculation
###################################### 
        my $tau = 0;
        my $SpeAnatEntity;
        my $maxTpmValueForThisGene = 0;
        for my $anatEntityId (keys %tpms){
        	if($tpms{$anatEntityId}{'maxOfMeanTPM'} > $maxTpmValueForThisGene){
        		$maxTpmValueForThisGene = $tpms{$anatEntityId}{'maxOfMeanTPM'};
        		$SpeAnatEntity = $anatEntityId;
        	}
        	$tau += 1 - ($tpms{$anatEntityId}{'maxOfMeanTPM'} / $maxXiTau);
        }
        $tau =$tau/((scalar keys %tpms) - 1);
        $tau{($genes{$bgeeGeneId})}{'value'} = $tau;
        $tau{($genes{$bgeeGeneId})}{'anatEntityId'} = $SpeAnatEntity;
    }
    print "Finish retrieving TPM values and calculation of Tau for each gene of species $speciesId\n";

    
######################################        
#       Create output File
###################################### 
	print "Max mean of log2(tpm) value for all genes : $maxLog2TPM\n";
	print "Start expression score normalization and write oncoMX expression file for species $speciesId\n";

    my $fileName = "$fileDir"."$outputFileNamePrefix"."$speciesId"."$outputFileNameSuffix";
    truncate $fileName, 0;
    open(my $fh, '>>', $fileName) or die "Can't open file '$fileName' $!";
    #write header
    say $fh "uniprotID\tanatEntityID\tscore\tdevStageID\tsex\tPres/Abs";
    close $fh;
    #Write expression data
 	foreach my $expr (@expr){
       	my $normelizedGeneExpression = ($expr->{'tpm'})*5/$maxLog2TPM;
       	open(my $fh, '>>', $fileName) or die "Can't open file '$fileName' $!";
       	if($tau{$expr->{'uniprotId'}}{'value'} >= 0.8 && $tau{$expr->{'uniprotId'}}{'anatEntityId'} eq $expr->{'anatEntityId'}){
       			my $score = $normelizedGeneExpression+5;
				say $fh "$expr->{'uniprotId'}\t$expr->{'anatEntityId'}\t$score\t$expr->{'stageId'}\t$expr->{'sex'}\t$expr->{'presence'}";
		}else{
			say $fh "$expr->{'uniprotId'}\t$expr->{'anatEntityId'}\t$normelizedGeneExpression\t$expr->{'stageId'}\t$expr->{'sex'}\t$expr->{'presence'}";
		}
		close $fh;
    }      
     
    print "Finish writing oncoMX expression file for species $speciesId\n";
    print "Start compression of the oncoMX expression file\n";
    gzip $fileName => "$fileName.gz"
        or die "gzip failed: $GzipError\n";
    print "Finish compression of the oncoMX expression file\n";
}


sub retrieveSwissprotGenes {
	my %genes;
	my @args = @_;
	my $speciesId = $args[0];
	print 'Start retrieving mapping between bgee gene IDs and uniprot/swissprot IDs for genes from species '.$speciesId."\n";
	my $sqlGenes =
    'select distinct g.bgeeGeneId,xref.XRefId '.
    'from expression e inner join gene g ON e.bgeeGeneId = g.bgeeGeneId '.
    'inner join geneXRef xref ON xref.bgeeGeneId = g.bgeeGeneId '. 
    'inner join rnaSeqResult rs ON e.expressionId = rs.expressionId '.
    'where g.speciesId = ? AND xref.dataSourceId = 5 and xref.XRefId NOT LIKE \'%#_%\' escape \'#\' order by xref.XRefId;';

    my $stmtGenes = $dbh->prepare($sqlGenes);
    $stmtGenes->execute($speciesId)  or die $stmtGenes->errstr;

    while ( my @data = $stmtGenes->fetchrow_array ){
        $genes{$data[0]}= $data[1];
    }
    my $i = 0;
    for my $testGenes (keys %genes){
    	$i = $i + 1;
    }
    return %genes;
}

sub retrieveChildrenDevStages {
	my %stages;
	my $devStageArg = $_[0];
	my $speciesId = $_[1];
	if($devStageArg ne $emptyArg){
	    print 'Start retrieving children dev. stage IDs of '.$devStageArg."\n";
        my $sqlStages = 'select stageId from stage '.
        'where stageName LIKE CONCAT(\'%\',(select speciesCommonName from species where speciesId = ?),\'%\') '.
        'AND stageLeftBound >= (select stageLeftBound from stage where stageId = ?) '.
        'AND stageRightBound <= (select stageRightBound from stage where stageId = ?);';

        my $stmtStages = $dbh->prepare($sqlStages);
        $stmtStages->execute($speciesId, $devStageArg, $devStageArg)  or die $stmtStages->errstr;

        while ( my @data = $stmtStages->fetchrow_array ){
            $stages{$data[0]} = 1;
        }
        #add parent dev stage term
        $stages{$devStageArg} = 1;
        print 'Finish retrieving children dev. stage IDs of '.$devStageArg."\n\t- Number of dev. Stages : ".(scalar keys %stages)."\n";
	}
	return %stages;
}
