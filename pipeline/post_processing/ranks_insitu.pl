#!/usr/bin/env perl
use strict;
use warnings;
use diagnostics;

# Calculate rawRanks for inSitu expression data per gene and condition
# Based on the scoring and ranking logic from ranks_global.pl

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
use Getopt::Long;

$|=1;

# Define arguments and their default value
my ($bgee_connector) = ('');
my ($output_file)    = ('');
my ($cond_id)        = ('');
my %opts = ('bgee=s'   => \$bgee_connector, # Bgee connector string
            'output=s' => \$output_file,    # Output file path (optional)
            'cond=s'   => \$cond_id,        # Specific conditionId (optional)
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD)
\t-bgee    Bgee connector string (required)
\t-output  Output file path (optional, default: STDOUT)
\t-cond    Specific conditionId (optional, default: all conditions)
\n";
    exit 1;
}

# Connect to database
my $dbh = Utils::connect_bgee_db($bgee_connector);

print "Starting inSitu rank calculation...\n";

# Prepare output filehandle
my $out_fh;
if ($output_file) {
    open($out_fh, '>', $output_file) or die "Cannot open output file $output_file: $!\n";
} else {
    $out_fh = *STDOUT;
}

# Print header
print $out_fh "conditionId\tbgeeGeneId\tscoreSum\trawRank\n";

# Get list of conditions to process
my @cond_list;
if ($cond_id) {
    @cond_list = ($cond_id);
} else {
    # Get all conditions that have inSitu data
    my $cond_query = $dbh->prepare(
        "SELECT DISTINCT cond.conditionId 
         FROM cond
         INNER JOIN inSituSpot ON cond.conditionId = inSituSpot.conditionId
         WHERE inSituSpot.expressionId IS NOT NULL
         ORDER BY cond.conditionId"
    );
    $cond_query->execute() or die $cond_query->errstr;
    while (my @row = $cond_query->fetchrow_array()) {
        push @cond_list, $row[0];
    }
    $cond_query->finish();
}

print "Processing " . scalar(@cond_list) . " condition(s)...\n";

# Process each condition
my $cond_count = 0;
foreach my $cond_id (@cond_list) {
    $cond_count++;
    
    # Calculate score for each gene in this condition
    # Score calculation:
    # present - high quality = 1
    # present - low quality = 0.5
    # absent - low quality = -0.5
    # absent - high quality = -1
    my $score_query = $dbh->prepare(
        "SELECT STRAIGHT_JOIN inSituSpot.bgeeGeneId,
            SUM(
                IF(inSituSpot.detectionFlag = ? AND inSituSpot.inSituData = ?, 1,
                    IF(inSituSpot.detectionFlag = ?, 0.5,
                        IF(inSituSpot.detectionFlag = ? AND inSituSpot.inSituData = ?, -1,
                            IF(inSituSpot.detectionFlag = ?, -0.5, 0))))
            ) AS scoreSum
         FROM inSituSpot
         WHERE inSituSpot.conditionId = ?
         AND inSituSpot.expressionId IS NOT NULL
         GROUP BY inSituSpot.bgeeGeneId
         ORDER BY scoreSum DESC"
    );
    
    $score_query->execute(
        $Utils::PRESENT_CALL, $Utils::HIGH_QUAL,
        $Utils::PRESENT_CALL,
        $Utils::ABSENT_CALL, $Utils::HIGH_QUAL,
        $Utils::ABSENT_CALL,
        $cond_id
    ) or die $score_query->errstr;
    
    # Fetch all results for this condition
    my @results;
    while (my @row = $score_query->fetchrow_array()) {
        push @results, {
            gene_id => $row[0],
            score => $row[1]
        };
    }
    $score_query->finish();
    
    # Apply dense ranking
    # Dense ranking means genes with the same score get the same rank,
    # and the next rank is the previous rank + 1 (not skipping any ranks)
    my $current_rank = 0;
    my $previous_score = undef;
    
    foreach my $result (@results) {
        if (!defined $previous_score || $result->{score} != $previous_score) {
            $current_rank++;
            $previous_score = $result->{score};
        }
        
        # Output: conditionId, bgeeGeneId, scoreSum, rawRank
        print $out_fh join("\t", 
            $cond_id, 
            $result->{gene_id}, 
            $result->{score}, 
            $current_rank
        ) . "\n";
    }
    
    if ($cond_count % 100 == 0) {
        print "Processed $cond_count / " . scalar(@cond_list) . " conditions...\n";
    }
}

print "Completed! Processed $cond_count condition(s).\n";

# Close output file if needed
if ($output_file) {
    close($out_fh);
    print "Results written to: $output_file\n";
}

$dbh->disconnect();

exit 0;
