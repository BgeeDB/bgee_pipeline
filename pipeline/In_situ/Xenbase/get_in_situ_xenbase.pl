#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Julien Roux, created 02/12/09
# USAGE: perl insert_in_situ_xenbase.pl
#######################################

use Getopt::Long;
use File::Slurp;
use LWP::Simple;

use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;

$| = 1;

# Define arguments & their default value
my ($bgee_connector)  = ('');
my ($Sport, $Aport)   = (0, 0);
my ($debug, $src_dir) = (0, '');
my %opts = ('debug'     => \$debug,            # more verbose
            'bgee=s'    => \$bgee_connector,   # Bgee connector string
            'Sport=i'   => \$Sport,            # Stage mapper socket port
            'Aport=i'   => \$Aport,            # Anatomy mapper socket port
            'src_dir=s' => \$src_dir,          # source_files/ path
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $Sport == 0 || $Aport == 0 ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD) -Sport=\$(INBETWEENSTAGESPORT) -Aport=\$(IDMAPPINGPORT) -src_dir=\$(SOURCE_FILES_DIR)\$(DIR_NAME)
\t-bgee      Bgee    connector string
\t-Sport     Stage   mapper socket port
\t-Aport     Anatomy mapper socket port
\t-src_dir   source_files/ path
\t-debug     More verbose
\n";
    exit 1;
}

# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);


# quality: not provided by Xenbase yet
# absence of expression: idem


###########################################
# Read and insert in situ data from Xenbase
###########################################
my $selSrc = $dbh->prepare('SELECT dataSourceId FROM dataSource WHERE dataSourceName = ?');
$selSrc->execute('Xenbase')  or die $selSrc->errstr;
my $data_source_id = $selSrc->fetchrow_array;
$selSrc->finish;
die "Data source Xenbase not found \n"  if ( !defined $data_source_id );

my %ensembl_gene;
my $retrieveEnsembl = $dbh->prepare('SELECT g.geneId, x.XRefId FROM geneXRef x, gene g WHERE x.dataSourceId = ? AND x.bgeeGeneId = g.bgeeGeneId AND g.speciesId=8364');
$retrieveEnsembl->execute($data_source_id)  or die $retrieveEnsembl->errstr;
while ( my @data = $retrieveEnsembl->fetchrow_array ){
    #xenbase ID -> Ensembl ID
    $ensembl_gene{$data[1]} = $data[0];
}
$retrieveEnsembl->finish;

print "Reading remote files...\n"  if ( $debug );
print "\tGeneExpression_tropicalis.txt (can be long)\n"  if ( $debug );
my $pattern_id = 1;
my %experiments;
my %patterns;

# all expression data for xenopus tropicalis
my $statusCode = mirror('ftp://ftp.xenbase.org/pub/GenePageReports/GeneExpression_tropicalis.txt', "$src_dir/GeneExpression_tropicalis.txt");
if ( $statusCode != 200 && $statusCode != 304 ){ # OKs & Not Modified
    die "Couldn't get file [GeneExpression_tropicalis.txt]!\n";
}

#NOTE No more XAO:\d+ in anatomy or stage terms !
system("perl -i -pe 's/XAO([0-9]{7,})/XAO:\$1/g; s/,([0-9]{7,})/,XAO:\$1/g' $src_dir/GeneExpression_tropicalis.txt");


# Get all start-end stages
my @Stages;
my @Anat;
for my $line ( read_file("$src_dir/GeneExpression_tropicalis.txt", chomp=>1) ){
    my @tmp = split("\t", $line);
    next  if ( lc($tmp[6]) ne 'in situ hybridization' );
    next  if ( lc($tmp[2]) ne 'wild type' );
    next  if ( !exists $ensembl_gene{$tmp[0]} );

    # Development stages mapping
    my $start_stage = '';
    if ( lc($tmp[4]) eq 'unspecified stage' ){
        $start_stage = 'XAO:1000081'; # unspecified stage XtroDv:0000086 #unknown
    }
    else {
        my @tmp2 = split(' ', $tmp[4]);
        $start_stage = $tmp2[0];
    }

    my $end_stage = '';
    if ( lc($tmp[5]) eq 'unspecified stage' ){
        $end_stage = 'XAO:1000081'; # unspecified stage XtroDv:0000086 #unknown
    }
    else {
        my @tmp2 = split(' ', $tmp[5]);
        $end_stage = $tmp2[0];
    }
    push @Stages, $start_stage, $end_stage;


    # Anatomy id mapping
    for my $couple ( split(',', $tmp[3]) ){
        next  if ( $couple =~ /^ / ); #Not a XAO id
        my @tmp3 = split(' ', $couple);
        #NOTE field syntax changed, XAO string is omitted
        #XAO:0000097 mandibular arch,0000098 hyoid arch,0000223 otic placode
        if ( $tmp3[0] !~ /^XAO:/ && $tmp3[0] =~ /^\d{7}$/ ){
            $tmp3[0] = 'XAO:'.$tmp3[0];
        }
        # XAO:0000059 is obsolete but may be used in ftp://ftp.xenbase.org/pub/GenePageReports/GeneExpression_tropicalis.txt  should be replaced by XAO:0002000
        $tmp3[0] =~ s{XAO:0000059}{XAO:0002000}  if ( $tmp3[0]     eq 'XAO:0000059' );
        $tmp3[0] = 'XAO:0003003'                 if ( lc($tmp3[0]) eq 'unspecified' );
        push @Anat, $tmp3[0];
    }
}
my $doneStages = Utils::get_in_between_stages(\@Stages, $Sport);
my $doneAnat   = Utils::get_anatomy_mapping(\@Anat,     $Aport);


SPOT:
for my $line ( read_file("$src_dir/GeneExpression_tropicalis.txt", chomp=>1) ){
    my @tmp = split("\t", $line);

    # Conditions:
    # if it is an in situ made on a wild type genotype
    # if the gene is mapped to Ensembl
    next  if ( lc($tmp[6]) ne 'in situ hybridization' );
    next  if ( lc($tmp[2]) ne 'wild type' );
    if ( !exists $ensembl_gene{$tmp[0]} ){
        warn "No gene mapping for [$tmp[0]]\n";
        next;
    }


    # split field 4 using comma
    my %organs;
    for my $couple ( split(',', $tmp[3]) ){
        my @tmp3 = split(' ', $couple);
        # XAO:0000059 is obsolete but may be used in ftp://ftp.xenbase.org/pub/GenePageReports/GeneExpression_tropicalis.txt  should be replaced by XAO:0002000
        $tmp3[0] =~ s{XAO:0000059}{XAO:0002000}  if ( $tmp3[0]     eq 'XAO:0000059' );
        $tmp3[0] = 'XAO:0003003'                 if ( lc($tmp3[0]) eq 'unspecified' );
        my $mapped_organ = $doneAnat->{$tmp3[0]} || '';
        warn "Unmapped organ: [$tmp3[0]]\n"  if ( $mapped_organ eq '' );
        $organs{ $mapped_organ } = ();
    }

    my $start_stage;
    if ( lc($tmp[4]) eq 'unspecified stage' ){
        $start_stage = 'XAO:1000081'; # unspecified stage XtroDv:0000086 #unknown
    }
    else {
        my @tmp2 = split(' ', $tmp[4]);
        $start_stage = $tmp2[0];
    }

    my $end_stage;
    if ( lc($tmp[5]) eq 'unspecified stage' ){
        $end_stage = 'XAO:1000081'; # unspecified stage XtroDv:0000086 #unknown
    }
    else {
        my @tmp2 = split(' ', $tmp[5]);
        $end_stage = $tmp2[0];
    }


    # Retrieve all stages between start and end stages
    my %stages;
    # stage mapping
    my $inbetweenstages = $start_stage.','.$end_stage;
    my $stages = $doneStages->{$inbetweenstages} || '';
    if ( $stages eq '' ){
        warn "Problem with in between stages for [$inbetweenstages] [$stages] in [$tmp[0]]\n";
        next SPOT;
    }

    STAGE_RANGE:
    for my $stage ( split(/\t/, $stages) ){
        # Skip this experiment if in between stage problem
        if ( $stage =~ /Exception/ ){
            chomp $stages;
            warn "Problem for [$inbetweenstages] [$stages] in [$tmp[0]]";
            next STAGE_RANGE  if ( $stage =~ /Could not find any OWLClass corresponding/ );
        }

        $stages{$stage} = {%organs};
    }


    my $experiment = $tmp[8];
    my $evidence   = $tmp[7];
    if ( defined $tmp[10] ){
        $experiments{$experiment} = $tmp[10] eq 'XB-ART-' ? '' : $tmp[10];
    }
    else {
        $experiments{$experiment} = '';
    }

    $patterns{$pattern_id}->{'gene'}       = $ensembl_gene{$tmp[0]};
    $patterns{$pattern_id}->{'experiment'} = $experiment;
    $patterns{$pattern_id}->{'evidence'}   = $evidence;
    $patterns{$pattern_id}->{'expression'} = {%stages};
    $patterns{$pattern_id}->{'strain'}     = $tmp[2];         #NOTE only "wild type" is kept and used in the file!
    $patterns{$pattern_id}->{'sex'}        = 'not annotated'; #NOTE not directly annotated
    # %organs are values of %stages
    $pattern_id++;
}


##################
# Fill Bgee tables
##################
print 'Filling inSituExperiment table... '  if ( $debug );
die "No experiment retrieved! Problem in the mapping of genes?\n"  if ( scalar keys %experiments eq 0 ); # can happen if Xenbase changed a file for example


print 'Filling inSituEvidence and inSituSpot tables... '  if ( $debug );
my %inserted_evidence;

my $spot_id = 0;
print join("\t", '#data_source', qw(inSituExperimentId  inSituEvidenceId  organId  stageId  geneId  detectionFlag  inSituData  linked  speciesId  strain  sex  ExprDesc)), "\n";
for my $pattern ( keys %patterns ){
    # insert evidences (if not already inserted before)
    if ( !exists $inserted_evidence{$patterns{$pattern}->{'evidence'}} ){
        $inserted_evidence{$patterns{$pattern}->{'evidence'}}++;
    }

    # insert spots
    for my $stage ( keys %{$patterns{$pattern}->{'expression'}} ){
        for my $organ ( keys %{$patterns{$pattern}->{'expression'}->{$stage}} ){
            next  if ( $organ eq '' );
            $spot_id++;
            print join("\t", $data_source_id,
                             $patterns{$pattern}->{'experiment'},
                             $patterns{$pattern}->{'evidence'},
                             $organ,
                             $stage,
                             $patterns{$pattern}->{'gene'},
                             'present',
                             'high quality',
                             $experiments{ $patterns{$pattern}->{'experiment'} },
                             8364,
                             $patterns{$pattern}->{'strain'} eq 'wild type' ? $Utils::WILD_TYPE_STRAIN : $patterns{$pattern}->{'strain'},
                             $patterns{$pattern}->{'sex'},
                             '',
                      ), "\n";
        }
    }
}
$dbh->disconnect;
print "Done\n"  if ( $debug );

exit 0;

