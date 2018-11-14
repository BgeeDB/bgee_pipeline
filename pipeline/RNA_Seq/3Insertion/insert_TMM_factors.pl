#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;
use File::Basename;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;
$| = 1;
require('rna_seq_utils.pl');

# Julien Roux, created October 2015; updated Nov 2016
# Insert TMM normalization factor for each RNA-Seq librarie.
# -debug: if provided, run in verbose mode (printing the update/insert SQL queries, not executing them)
#############################################################

# Define arguments & their default value
my ($bgee_connector) = ('');
my ($tmm_results)    = ('');
my ($debug)          = (0); # = former -v option
my %opts = ('bgee=s'                => \$bgee_connector,   # Bgee connector string
            'tmm_results=s'         => \$tmm_results,     # /var/bgee/extra/pipeline/rna_seq/processed_TMM_bgee_v14/
            'debug|v'               => \$debug,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $tmm_results eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD) -ensRelease=\$(ENSRELEASE)  -tmm_results=\$(RNASEQTMMPATH)
\t-bgee             Bgee connector string
\t-tmm_results      TMM factors results directory
\t-debug            more verbose output
\n";
    exit 1;
}

# Bgee db connection
my $bgee = Utils::connect_bgee_db($bgee_connector);

####################
# READ TMM FACTORS #
####################
print "Reading TMM factors results from files... \n";

my %toInsert;
for my $file ( glob($tmm_results.'/*.tsv') ){
    $file = basename($file);

    open(my $IN0, '<', $tmm_results.$file)  or die "Can't read file [$file]\n";
    my $line = <$IN0>; #header
    while ( defined ($line = <$IN0>) ){
        chomp $line;
        my @tmp = split(/\t/, $line);
        # rnaSeqLibraryId -> TMM factor
        $toInsert{$tmp[1]}->{'tmmFactor'} = $tmp[2];
    }
    close $IN0;
}

######################
# INSERT TMM FACTORS #
######################
print "Inserting TMM factors...\n";

my $insTMM = $bgee->prepare('UPDATE rnaSeqLibrary SET tmmFactor = ? WHERE rnaSeqLibraryId = ?;');
for my $library ( keys %toInsert ){
    print "\t$library\n";

    # insert the TMM factor rounded to 6 decimals
    if ( $debug ){
        print "\tUPDATE rnaSeqLibrary SET tmmFactor = ", int( $toInsert{$library}->{'tmmFactor'}*(10**6) +.5 ) / (10**6), ' WHERE rnaSeqLibraryId = ', $library, ";\n";
    }
    else {
        $insTMM->execute(int( $toInsert{$library}->{'tmmFactor'}*(10**6) +.5 ) / (10**6), $library)  or die $insTMM->errstr;
    }
}
$insTMM->finish();
print "Done\n";

$bgee->disconnect;
exit 0;

