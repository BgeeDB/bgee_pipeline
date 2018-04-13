#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Roelli Patrick, created 12.06.12
# USAGE perl complete_mapping_from_bdgp_to_bgee.pl <old_mapping_file.tsv> <new_mapping_file.tsv>

use Getopt::Long;
use List::MoreUtils qw{uniq};
use File::Slurp;
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;

$| = 1;

# Define arguments & their default value
my ($bgee_connector, $bdgp_connector) = ('', '');
my ($debug, $old, $new, $annot_info)  = (0, '', '', '');
my %opts = ('debug'         => \$debug,            # more verbose
            'bgee=s'        => \$bgee_connector,   # Bgee connector string
            'bdgp=s'        => \$bdgp_connector,   # BDGP connector string
            'old=s'         => \$old,              # old_mapping_file.tsv
            'new=s'         => \$new,              # new_mapping_file.tsv
            'annot_info=s'  => \$annot_info,       # stage info
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $bdgp_connector eq '' || $new eq '' || $old eq '' || $annot_info eq '' ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD) -bdgp=\$(BDGPCMD) -old=old_mapping_file.tsv -new=new_mapping_file.tsv
\t-bgee        Bgee       connector string
\t-bdgp        BDGP local connector string
\t-debug       More verbose
\t-old         Old mapping file to read results from
\t-new         New mapping file to put  results into
\t-annot_info  Annotation info
\n";
    exit 1;
}

# Bgee db connection
my $dbh  = Utils::connect_bgee_db($bgee_connector);
# BDGP db connection
my $BDGP = Utils::connect_bgee_db($bdgp_connector);


################################
print "running BDGP mapping\n"  if ( $debug );

# Get all terms from BGDP using the id as the main key, stage as secondary
my %BDGP_terms;
# The $query1 gets all terms from BDGP and their corresponding stage
my $selTerm = $BDGP->prepare('SELECT DISTINCT S.name AS stage_name, T.go_term AS BDGP_term, T.id AS BDGP_id, S.id AS stage_id FROM term AS T
                          INNER JOIN annot_term AS AT ON AT.term_id = T.id
                          INNER JOIN annot AS A ON AT.annot_id = A.id
                          INNER JOIN stage AS S ON S.id = A.stage ORDER BY BDGP_id');
$selTerm->execute()  or die $selTerm->errstr;
while ( my @data = $selTerm->fetchrow_array ){
    # $data[2]=ids, $data[0]=stages
    $BDGP_terms{$data[2]}{$data[0]}{'stage_name'}    = $data[0];
    $BDGP_terms{$data[2]}{$data[0]}{'stage_id'}      = $data[3];
    $BDGP_terms{$data[2]}{$data[0]}{'BDGP_term'}     = $data[1];
    $BDGP_terms{$data[2]}{$data[0]}{'FBbt_id'}       = 'none';
    $BDGP_terms{$data[2]}{$data[0]}{'FBbt_name'}     = 'none';
    $BDGP_terms{$data[2]}{$data[0]}{'matches_count'} = 0;
    $BDGP_terms{$data[2]}{$data[0]}{'match_source'}  = 'nomatch';
}
$selTerm->finish;
$BDGP->disconnect;
# Can't count only the ids since there are multiple entry by id with different stages (might be fixed one day in BDGP, one can only wish!: 22.06.12)
my $count_bdgp_terms = 0;
for my $id ( keys %BDGP_terms ){
    $count_bdgp_terms += scalar(keys %{$BDGP_terms{$id}});
}

my %tsv_annot = %{ Utils::read_spreadsheet("$annot_info", "\t", 'csv', '', 1) };
my $annot;
for my $line ( 0..$#{$tsv_annot{'count'}} ){
    $annot->{ $tsv_annot{'stage_id'}[$line] .'--'. $tsv_annot{'organ name'}[$line] } = $tsv_annot{'count'}[$line];
}


# Extract all lines from the manual mapping TSV file using the name of the ID as a key in the hash
my %organs_mapping; # The hash has an array of two elements linked to each key-> BDGP_id : ([0]=FBbt_id, [1]=FBbt_name)
my %tsv = %{ Utils::read_spreadsheet("$old", "\t", 'csv', '"', 1) };
# file structure: stage_name    BGDP_term    BDGP_id    FBbt_id    FBbt_name    matches_count    match_source    stage_id    occurrence
for my $line ( 0..$#{$tsv{'stage_name'}} ){
    # If the FBbt_id is not = 'none' and this FBbt_id is defined in bgee
#    if ( $tsv{'FBbt_id'}[$line] ne 'none' && exists $organs_fbbt_name{ $tsv{'FBbt_id'}[$line] } ){
#        @{ $organs_mapping{ $tsv{'BDGP_id'}[$line] }{ $tsv{'stage_name'}[$line] } } = ($tsv{'FBbt_id'}[$line], $organs_fbbt_name{ $tsv{'FBbt_id'}[$line] }); # Array assigned to a hash with key = $BDGP_id
#    }
#    elsif ( $tsv{'FBbt_id'}[$line] ne 'none' ){
    if ( $tsv{'FBbt_id'}[$line] ne 'none' ){
        # The "stage" table is no longer already filled at this step, so no stage name, only stage id
        @{ $organs_mapping{ $tsv{'BDGP_id'}[$line] }{ $tsv{'stage_name'}[$line] } } = ($tsv{'FBbt_id'}[$line], ''); # Array assigned to a hash with key = $BDGP_id
    }
}

# Hashes created:
# 1) BDGP_terms with BDGP id as main key, stage as secondary, main hash holding the states of the mappings
## 2) organs_bgee with organ name as a key, will be used to search for matches
# 3) organs_mapping with BDGP id as a key, representing mappings found in the TSV mapping file

# just a counter to display log messages
my $count_f = 1;

my $total         = $count_bdgp_terms;
my $match_mapping = 0;
my $match_bgee    = 0;
# Check all terms from BDGP that match the mapping
for my $id ( keys %BDGP_terms ){
    for my $stage ( keys %{$BDGP_terms{$id}} ){
        if ( $debug && $count_f%100==0 ){ print $count_f."\n"; } # To see that the program is running each 100 terms it will print the term number
        if ( exists $organs_mapping{$id} ){
            $BDGP_terms{$id}{$stage}{'FBbt_id'}        = $organs_mapping{$id}{$stage}[0];
            $BDGP_terms{$id}{$stage}{'FBbt_name'}      = $organs_mapping{$id}{$stage}[1];
            $BDGP_terms{$id}{$stage}{'matches_count'} += 1;
            $BDGP_terms{$id}{$stage}{'match_source'}   = 'mapping';
            $total--;
            $match_mapping++;
        }
        $count_f++;
    }
}


# Sorting based on the match_source
my @nomatch;
my @bgeematch;
my @mappingmatch;
for my $id ( keys %BDGP_terms ){
    for my $stage ( keys %{$BDGP_terms{$id}} ){
        # if there is an ambiguity in the match (more than one match)
        # do not use this mapping
        if ( $BDGP_terms{$id}{$stage}{'matches_count'} > 1 ){
            $BDGP_terms{$id}{$stage}{'FBbt_id'}      = 'none';
            $BDGP_terms{$id}{$stage}{'FBbt_name'}    = 'none';
            $BDGP_terms{$id}{$stage}{'match_source'} = 'ambiguous';
        }

        if ( $BDGP_terms{$id}{$stage}{'match_source'} eq 'mapping' ){
            push(@mappingmatch, $id);
        }
        elsif ( $BDGP_terms{$id}{$stage}{'match_source'} eq 'bgee' ){
            push (@bgeematch, $id);
        }
        else{
            push (@nomatch, $id);
        }
    }
}


my @finalarray = (uniq sort @nomatch, uniq sort @bgeematch, uniq sort @mappingmatch);
# Create a new Excel workbook for the mapping (, add a worksheet and define a format)
my $output = join("\t", '#stage_name', qw(BGDP_term BDGP_id FBbt_id FBbt_name matches_count match_source stage_id occurence))."\n";
for my $id ( sort { $a <=> $b } @finalarray ){
    for my $stage ( sort keys %{$BDGP_terms{$id}} ){
        my $occurence = $annot->{ $BDGP_terms{$id}{$stage}{'stage_id'} .'--'. $BDGP_terms{$id}{$stage}{'BDGP_term'} } // '';
        $output .= join("\t", $BDGP_terms{$id}{$stage}{'stage_name'},
                              $BDGP_terms{$id}{$stage}{'BDGP_term'},
                              $id,
                              $BDGP_terms{$id}{$stage}{'FBbt_id'},
                              $BDGP_terms{$id}{$stage}{'FBbt_name'},
                              $BDGP_terms{$id}{$stage}{'matches_count'},
                              $BDGP_terms{$id}{$stage}{'match_source'},
                              $BDGP_terms{$id}{$stage}{'stage_id'},
                              $occurence,
                       )."\n";
    }
}
write_file($new, $output);

# summary
print "Number of terms retrieved from BDGP: $count_bdgp_terms\n"  if ( $debug );
# print "The number of elements from the mapping is: $count_mapping\n";
print "$match_mapping matched the mapping, $match_bgee matched bgee, $total unmatched\n"  if ( $debug );

exit 0;

