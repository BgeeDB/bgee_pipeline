package Utils;
#File Utils.pm

# Perl core modules
use strict;
use warnings;
use diagnostics;

use Array::Utils qw(:all);
use Bio::EnsEMBL::Registry; # Require Ensembl API
use DBI;
use File::Basename;
use File::Slurp;
use IO::Socket;
use List::MoreUtils qw(uniq);
use Spreadsheet::Read qw{ReadData row};

# flush after every write
$| = 1;

# Define some variables to be used through the whole pipeline to consistently capture
# expression calls and qualities.
# call types
our $PRESENT_CALL   = 'present';
our $ABSENT_CALL    = 'absent';
our $UNDEFINED_CALL = 'undefined';
our $MARGINAL_CALL  = 'marginal';
# qualities
our $LOW_QUAL  = 'poor quality'; # 'poor' instead of 'low' for historical reasons
our $HIGH_QUAL = 'high quality';
# Reasons for exclusion
our $CALL_NOT_EXCLUDED             = 'not excluded';
our $EXCLUDED_FOR_UNDEFINED        = 'undefined';
our $EXCLUDED_FOR_PRE_FILTERED     = 'pre-filtering';
our $EXCLUDED_FOR_NO_EXPR_CONFLICT = 'noExpression conflict';


# Define some variables to be used through the whole pipeline to consistently capture sex info,
# and some strain info.
our $FEMALE_SEX          = 'female';
our $MALE_SEX            = 'male';
our $HERMAPHRODITE_SEX   = 'hermaphrodite';
our $MIXED_SEXES         = 'mixed';
our $NOT_ANNOTATED_SEX   = 'not annotated';
our $NA_SEX              = 'NA';
our @ACCEPTABLE_SEX_INFO = ($MALE_SEX, $FEMALE_SEX, $HERMAPHRODITE_SEX, $MIXED_SEXES, $NOT_ANNOTATED_SEX, $NA_SEX);
our @NO_ANNOT_SEX_INFO   = ($NOT_ANNOTATED_SEX, $NA_SEX);

our $WILD_TYPE_STRAIN         = 'wild-type';
our $NA_STRAIN                = 'NA';
our $NOT_ANNOTATED_STRAIN     = 'not annotated';
our $CRD_STRAIN               = 'confidential_restricted_data';
our @STANDARDIZED_STRAIN_INFO = ($WILD_TYPE_STRAIN, $NA_STRAIN, $NOT_ANNOTATED_STRAIN, $CRD_STRAIN);
our @NO_ANNOT_STRAIN_INFO     = ($WILD_TYPE_STRAIN, $NOT_ANNOTATED_STRAIN, $NA_STRAIN, $CRD_STRAIN);


our $EST_DATA_TYPE     = 'est';
our $AFFY_DATA_TYPE    = 'affymetrix';
our $IN_SITU_DATA_TYPE = 'inSitu';
our $RNA_SEQ_DATA_TYPE = 'rnaSeq';
our @DATA_TYPES        = ($EST_DATA_TYPE, $AFFY_DATA_TYPE, $IN_SITU_DATA_TYPE, $RNA_SEQ_DATA_TYPE);


# Variables used for computing ranks
our $ANAT_ENTITY_PARAM    = 'anatEntityId';
our $STAGE_PARAM          = 'stageId';
our $SEX_PARAM            = 'sex';
our $STRAIN_PARAM         = 'strain';
our @CONDITION_PARAMETERS = ($ANAT_ENTITY_PARAM, $STAGE_PARAM, $SEX_PARAM, $STRAIN_PARAM);


# A sub to retrieve condition parameter combinations with no rank yet computed
# for the requested data type.
# This sub takes as argument an array where the first element is a connection to the Bgee database
# (see sub 'connect_bgee_db'), and the second element is the requested data type
# (see variable @DATA_TYPES for an array of valid data types).
# It returns a reference to an array of array references, where each inner array is a condition parameter combination,
# containing the desired condition parameters for this combination (see variable @CONDITION_PARAMETERS
# for an array of valid condition parameters).
sub get_cond_param_combinations {
    my ($dbh, $dataType)= @_;

    if ( !grep( /^$dataType$/, @DATA_TYPES ) ) {
        die "Unrecognized data type: $dataType\n";
    }

    my @condParamCombinations = ();
    my $sql = "
    SELECT DISTINCT
        IF(t1.anatEntityId IS NULL, 0, 1) AS anatEntityParam,
        IF(t1.stageId IS NULL, 0, 1) AS stageParam,
        IF(t1.sex IS NULL, 0, 1) AS sexParam,
        IF(t1.strain IS NULL, 0, 1) AS strainParam
    FROM globalCond AS t1 "
#    .
#    # identify condition parameter combinations with no rank already computed
#    #for the requested data type
#    "WHERE NOT EXISTS (
#        SELECT 1 FROM globalCond AS t2
#        WHERE ";
#    if ($dataType eq $EST_DATA_TYPE) {
#        $sql .= "estMaxRank IS NOT NULL";
#    } elsif ($dataType eq $AFFY_DATA_TYPE) {
#        $sql .= "affymetrixMaxRank IS NOT NULL";
#    } elsif ($dataType eq $IN_SITU_DATA_TYPE) {
#        $sql .= "inSituMaxRank IS NOT NULL";
#    } elsif ($dataType eq $RNA_SEQ_DATA_TYPE) {
#        $sql .= "rnaSeqMaxRank IS NOT NULL";
#    } else {
#        die "Unsupported data type: $dataType\n";
#    }
#    $sql .= " AND (t1.anatEntityId IS NULL AND t2.anatEntityId IS NULL OR
#                 t1.anatEntityId IS NOT NULL AND t2.anatEntityId IS NOT NULL) AND
#             (t1.stageId IS NULL AND t2.stageId IS NULL OR
#                 t1.stageId IS NOT NULL AND t2.stageId IS NOT NULL) AND
#             (t1.sex IS NULL AND t2.sex IS NULL OR
#                 t1.sex IS NOT NULL AND t2.sex IS NOT NULL) AND
#             (t1.strain IS NULL AND t2.strain IS NULL OR
#                 t1.strain IS NOT NULL AND t2.strain IS NOT NULL)
#    )"
    ;

    my $query = $dbh->prepare($sql);
    $query->execute()  or die $query->errstr;
    while ( my $dataRef = $query->fetchrow_hashref ){
        my @localCondParamComb = ();
        for my $column ( keys %{$dataRef} ){
            if ( $column eq 'anatEntityParam' ){
                if ( $dataRef->{$column} == 1 ) {
                    push @localCondParamComb, $ANAT_ENTITY_PARAM;
                }
            } elsif ( $column eq 'stageParam' ){
                if ( $dataRef->{$column} == 1 ) {
                    push @localCondParamComb, $STAGE_PARAM;
                }
            } elsif ( $column eq 'sexParam' ){
                if ( $dataRef->{$column} == 1 ) {
                    push @localCondParamComb, $SEX_PARAM;
                }
            } elsif ( $column eq 'strainParam' ){
                if ( $dataRef->{$column} == 1 ) {
                    push @localCondParamComb, $STRAIN_PARAM;
                }
            } else {
                die "Unsupported column: $column\n";
            }
        }
        if (!@localCondParamComb) {
            die "No condition parameter combination retrieved for row: "
                ."$_ $dataRef->{$_} - " for (keys %{$dataRef})."\n";
        }
        push @condParamCombinations, \@localCondParamComb;
    }

    return \@condParamCombinations;
}

# This sub returns the SQL "AND" conditions to select a specific condition parameter combination
# from the globalCond table. It takes as first argument a reference to an array containing
# the condition parameters to consider (see variable @CONDITION_PARAMETERS for an array of
# valid condition parameters), and as second argument the name of the table 'globalCond'
# in the SQL query. It returns a String containing 'AND' conditions to be used in the SQL query.
sub get_cond_param_comb_sql_clause {
    my ($condParamCombRef, $tableName)= @_;

    my @condParamComb = @{$condParamCombRef};
    my $sql = "";
    my $i = 0;
    my $paramFound = 0;
    for my $condParam (@CONDITION_PARAMETERS) {
        if ($i > 0) {
            $sql .= " AND ";
        }
        $sql .= $tableName.".".$condParam;
        if ( grep( /^$condParam$/, @condParamComb ) ) {
            $sql .= " IS NOT NULL";
            $paramFound = 1;
        } else {
            $sql .= " IS NULL";
        }
        $i++;
    }
    if (!$paramFound) {
        die "No valid condition parameter provided in the combination.\n";
    }

    return $sql;
}

sub get_fields_for_sql_select {
    my ($fields_ref) = @_;
    my @fields = @{ $fields_ref };

    my $sql = '';
    my $i = 0;
    for my $field ( @fields ){
        if ($i > 0) {
            $sql .= ', ';
        }
        $sql .= $field;
        $i++;
    }
    return $sql;
}

sub get_fields_for_sql_join {
    my ($fields_ref, $tableName1, $tableName2) = @_;
    my @fields = @{ $fields_ref };

    my $sql = '';
    my $i = 0;
    for my $field ( @fields ){
        if ($i > 0) {
            $sql .= ' AND ';
        }
        $sql .= $tableName1.'.'.$field.' = '.$tableName2.'.'.$field;
        $i++;
    }
    return $sql;
}

## Connectors, from parsed command line
sub connect_ensembl_registry {
    my ($ensembl_connector, $verbose) = @_;
    $verbose = $verbose || 0;

    # Get ensembl db parameters from command line
    my %ensembl = map { /^(.+?)=(.*?)$/; $1 => $2 }   # Split based on =  and assigned as key=>value
                  split('__', $ensembl_connector);    # Split based on __

    # Connection to Ensembl via Ensembl API/Registry
    my $reg = 'Bio::EnsEMBL::Registry';
    $reg->load_registry_from_db( -user    => $ensembl{'user'},
                                 -pass    => $ensembl{'pass'},
                                 -host    => $ensembl{'host'},
                                 -port    => $ensembl{'port'},
                                 -verbose => $verbose, # Verbose output
                               );

    return $reg;
}

sub connect_bgee_db {
    my ($bgee_connector)= @_;

    # Get bgee db parameters from command line
    my %bgee = map { /^(.+?)=(.+?)$/; $1 => $2 }   # Split based on =  and assigned as key=>value
               split('__', $bgee_connector);       # Split based on __

    # Bgee db connection
    my $dbh = DBI->connect("dbi:mysql:database=$bgee{'name'};host=$bgee{'host'};port=$bgee{'port'}",
                $bgee{'user'}, $bgee{'pass'})  or die $DBI::errstr;

    return $dbh;
}

sub connect_mgi_db {
    my ($mgi_connector)= @_;

    # Get MGI db parameters from command line
    my %mgi = map { /^(.+?)=(.+?)$/; $1 => $2 }   # Split based on =  and assigned as key=>value
              split('__', $mgi_connector);       # Split based on __

    # MGI db connection
    my $dbh = DBI->connect("dbi:Pg:database=$mgi{'name'};host=$mgi{'host'};port=$mgi{'port'}",
                $mgi{'user'}, $mgi{'pass'})  or die $DBI::errstr;

    return $dbh;
}


## spreadsheet reader (xls, xlsx, ods, xsc, tsv, csv)
# and put columns in hash: header == key, rows == list
## TO DO: this function is really messy. Needs to be cleaned up for clarity. It is also probably avoidable to write a temporary file
sub read_spreadsheet {
    my ($file, $separator, $parser, $quote, $sheet) = @_;

    my @infile  = read_file("$file");
    my $tmpfile = '/tmp/'.basename($file).$$.'.'.$parser;
    # FIXME Trick to avoid parsing commented lines
    # Should be done in a better way but the current parsing is by column, not by row
    # And the header line also starts by # but we want to keep it!
    write_file("$tmpfile", join('', $infile[0], grep { !/^$quote?#/ } @infile[1..$#infile]));

    #NOTE sometimes csv parser has to be forced because not detected as it (file extension ???)
    # quote => $quote  for no quote around fields to avoid conflicts with ' and " already there
    # strip => 3       to remove trailing- and leading-whitespace from every field
    my $book   = ReadData ("$tmpfile", sep => $separator, parser => $parser, quote => $quote, strip => 3);
    unlink "$tmpfile";
    #NOTE rows & columns indexes start at 1 to be Excel cell name compliant
    my @header = row($book->[$sheet], 1);         # Get header lines, from sheet 1
    my @Col    = @{$book->[$sheet]->{'cell'}};    # Get column content
    my %tsv;
    for my $col ( 0..$#header ){
        my $colu = $col;
        #NOTE To be able to assign header key name to column, header line must start with a #
        # or no easy way to guess if there are headers or not
        if ( $header[0] =~ /^#/ ){ # Real headers
            $colu = $header[$col];
            $colu =~ s{^#}{};
        }
        shift @{ $Col[$col+1] }; # Remove first col value because is always undef to start spreadsheet cell number at 1 as in real spreadsheet
        shift @{ $Col[$col+1] }  if ( $Col[$col+1][0] =~ /^#?$header[$col]$/ ); # Remove header, if any
        # Replace double or triple ' ' by a single ' '
        my @column = map { $_ =~ s/  +/ /g  if defined $_; $_ }
                     @{ $Col[$col+1] };
        $tsv{$colu} = \@column;
    }
    return \%tsv; # Return hash pointer/reference
}

sub trim {
    my ($string) = @_;

    if ( defined $string ){
        $string =~ s{^\s+|\s+$}{}g;
        $string =~ s{\s\s+}{ }g;
    }

    return $string;
}


# A sub to reconcile expressed/not expressed calls
# defineCallAndQuality was too harsh as a scoring scheme, since, e.g., only a single probeset
# as pst/low was sufficient to decrease the quality to low.
# This new scoring system should be more relaxed
# WARNING, METHOD NOT TO USE ANYMORE, THE REASONS FOR EXCLUSION HAVE CHANGED AS OF BGEE 14
sub defineCallAndQualityNew {
    # A ref to a hash of hashes with two keys:
    # 'call' and 'quality'
    # The elements can come for instance from a list of in situ spots,
    # or Affymetrix probesets, or Affymetrix summarized calls per chip, or RNA-Seq results
    my ($callsRef) = @_;

    # New scoring scheme is based on two scores:
    # - Sum presence evidences to give a presence score:
    #   pst high / pst low
    #   +1      /  +0.5
    # - Sum absence evidences to give an absence score
    #   abs high / abs low
    #   +1         +0.5

    my $pst_score = 0;
    my $abs_score = 0;
    my $pstHigh   = 0; # a boolean to know if there is at least on pst/high call

    for my $key ( keys %$callsRef ){
        my $flag    = $callsRef->{$key}->{'call'};
        my $quality = $callsRef->{$key}->{'quality'};

        if (    $flag eq 'present' && $quality eq 'high quality' ){
          $pst_score += 1;
          $pstHigh    = 1;
        }
        elsif ( $flag eq 'present' && $quality eq 'poor quality' ){
          $pst_score += 0.5;
        }
        elsif ( $flag eq 'absent'  && $quality eq 'high quality' ){
          $abs_score += 1;
        }
        elsif ( $flag eq 'absent'  && $quality eq 'poor quality' ){
          $abs_score += 0.5;
        }
    }
    my $summary;

    # New scoring scheme:
    # Scores can be 0, 0.5 and >=1.
    #
    # Note on pst/bronze quality: until Bgee 13 included, bronze quality used to be
    # no pst high + some pst low vs. some abs low and/or high => bronze qual. => not inserted
    # Now it is slightly changed: if no pst high + some pst,
    # BUT pst score > abs score -> pst low (previoulsy, this case was not "rescued").
    # Also, as a result, a sum of pst/low can only give a pst/high summary if abs score == 0
    # (previously, a sum of pst/low was never giving a pst/high summary even if abs score == 0,
    # and was always discarded as bronze quality if abs score > 0)
    #
    # Ratio rule:
    # if pst score >= 1 && abs score == 0                      -> pst high
    # if some pst/high               && pst score >  abs score -> pst high
    # if some pst/high               && pst score <= abs score -> pst low
    # if no pst/high && some pst/low && pst score >  abs score -> pst low (no abs + pst score = 0.5, or bronze "rescued")
    # if no pst/high && some pst/low && pst score <= abs score -> pst bronze qual (not inserted)
    # if pst score == 0              && abs score >= 1         -> abs high
    # if pst score == 0              && abs score == 0.5       -> abs low (not inserted)

    ## ratio rules
    if ( $pstHigh && $pst_score > $abs_score ||
         $pst_score >= 1 && $abs_score == 0 ) {
      $summary->{'call'}      = 'present';
      $summary->{'quality'}   = 'high quality';
      $summary->{'exclusion'} = $CALL_NOT_EXCLUDED;
    }
    elsif ( $pstHigh && $pst_score <= $abs_score ||
            $pst_score > 0 && $pst_score > $abs_score ) {
      $summary->{'call'}      = 'present';
      $summary->{'quality'}   = 'poor quality';
      $summary->{'exclusion'} = $CALL_NOT_EXCLUDED;
    }
    elsif ( $pst_score > 0 && $pst_score <= $abs_score ) {
      $summary->{'call'}      = undef;
      $summary->{'quality'}   = undef;
      $summary->{'exclusion'} = 'bronze quality'; #Note: not discarded anymore in Bgee 14
    }
    elsif ( $pst_score == 0 && $abs_score >= 1 ) {
      $summary->{'call'}      = 'absent';
      $summary->{'quality'}   = 'high quality';
      $summary->{'exclusion'} = $CALL_NOT_EXCLUDED;
    }
    elsif ( $pst_score == 0 && $abs_score > 0 && $abs_score < 1) {
      $summary->{'call'}      = undef;
      $summary->{'quality'}   = undef;
      $summary->{'exclusion'} = 'absent low quality'; #Note: not discarded anymore in Bgee 14
    }
    else {
      # only 'undefined' calls have been seen
      $summary->{'call'}      = undef;
      $summary->{'quality'}   = undef;
      $summary->{'exclusion'} = $EXCLUDED_FOR_UNDEFINED;
    }
    return $summary;
}

# A sub to extract relevant information, from a list of calls, for computing expression calls
# and qualities. This sub does not directly compute a call and a quality, but retrieve the information
# allowing to compute them, to be inserted into the database, for later on-the-fly computations.
# The calls should all be about the same gene and condition.
#
# Argument of the sub: reference to a hash of hashes of hashes
# experimentId -> evidenceId -> 'call'/'quality'
# * values of 'call' must be one of $Utils::PRESENT_CALL and  $Utils::ABSENT_CALL,
# other values not taken into account.
# * values of 'quality' must be one of $Utils::HIGH_QUAL and  $Utils::LOW_QUAL.
#
# Returned value: reference to an array with two values.
# First value: a string providing a reason for exclusion (e.g., $Utils::EXCLUDED_FOR_UNDEFINED).
# Second value: reference to a hash of hashes
# experimentId -> 'pstHighEvidenceCount'/'pstLowEvidenceCount'/'absHighEvidenceCount'/'absLowEvidenceCount'
#                 /'expCall'/'expCallQuality'.
# 'expCall' is one of $Utils::PRESENT_CALL and  $Utils::ABSENT_CALL, and provides
# the global expression call produced by this experiment for this gene-condition.
# 'expCallQuality' is one of $Utils::HIGH_QUAL and  $Utils::LOW_QUAL, and provides
# the global call quality produced by this experiment for this gene-condition.
sub summarizeExperimentCallAndQuality {
    my ($callsRef) = @_;

    my %summary = ();
    EXP: for my $expId ( keys %{$callsRef} ) {
        my $pstHighEvidenceCout = 0;
        my $pstLowEvidenceCout  = 0;
        my $absHighEvidenceCout = 0;
        my $absLowEvidenceCout  = 0;

        for my $evidenceId ( keys %{$callsRef->{$expId}} ){
            my $flag    = $callsRef->{$expId}->{$evidenceId}->{'call'};
            my $quality = $callsRef->{$expId}->{$evidenceId}->{'quality'};

            if (    $flag eq $PRESENT_CALL && $quality eq $HIGH_QUAL ){
                $pstHighEvidenceCout += 1;
            }
            elsif ( $flag eq $PRESENT_CALL && $quality eq $LOW_QUAL ){
                $pstLowEvidenceCout += 1;
            }
            elsif ( $flag eq $ABSENT_CALL  && $quality eq $HIGH_QUAL ){
                $absHighEvidenceCout += 1;
            }
            elsif ( $flag eq $ABSENT_CALL  && $quality eq $LOW_QUAL ){
                $absLowEvidenceCout += 1;
            }
        }
        if ($pstHighEvidenceCout == 0 && $pstLowEvidenceCout == 0 &&
               $absHighEvidenceCout == 0 && $absLowEvidenceCout == 0) {
            next EXP;
        }

        my $expCall     = undef;
        my $expCallQual = undef;

        if ($pstHighEvidenceCout > 0) {
            $expCall     = $PRESENT_CALL;
            $expCallQual = $HIGH_QUAL;
        }
        elsif ($pstLowEvidenceCout > 0) {
            $expCall     = $PRESENT_CALL;
            $expCallQual = $LOW_QUAL;
        }
        elsif ($absHighEvidenceCout > 0) {
            $expCall     = $ABSENT_CALL;
            $expCallQual = $HIGH_QUAL;
        }
        elsif ($absLowEvidenceCout > 0) {
            $expCall     = $ABSENT_CALL;
            $expCallQual = $LOW_QUAL;
        }
        if (!defined $expCall || !defined $expCallQual) {
            die "Impossible call summary for experiment $expId\n";
        }

        $summary{$expId}->{'pstHighEvidenceCount'} = $pstHighEvidenceCout;
        $summary{$expId}->{'pstLowEvidenceCount'}  = $pstLowEvidenceCout;
        $summary{$expId}->{'absHighEvidenceCount'} = $absHighEvidenceCout;
        $summary{$expId}->{'absLowEvidenceCount'}  = $absLowEvidenceCout;
        $summary{$expId}->{'expCall'}              = $expCall;
        $summary{$expId}->{'expCallQuality'}       = $expCallQual;
    }
    my $reasonForExclusion  = $CALL_NOT_EXCLUDED;
    if ( !%summary ) {
        $reasonForExclusion = $EXCLUDED_FOR_UNDEFINED;
    }

    return ($reasonForExclusion, \%summary);
}


# Sub-query to map the actually annotated stage to a not too granular stage
# (e.g., in human, we don't want to compare organs at stage "83-year old",
# but to consider all info at stage "80 year-old and over human stage").
sub get_stage_equivalences_query {
    return 'SELECT DISTINCT t1.stageId,
             (SELECT t2.stageId
                FROM stage AS t2 '.
                # Use taxon constraints to make sure to get a parent stage valid in the related species
                'INNER JOIN stageTaxonConstraint AS t2bis ON t2.stageId = t2bis.stageId '.
                # Get the stage itself or its parent (left bound - right bound)
                'WHERE t2.stageLeftBound <= t1.stageLeftBound AND t2.stageRightBound >= t1.stageRightBound '.
                # in the proper species
                'AND (t2bis.speciesId IS NULL OR t2bis.speciesId = t1bis.speciesId) '.
                # that is not too granular (tooGranular = 0),
                # and that is the closest to the annotated stage (left bound desc order limit 1)
                'AND t2.tooGranular= 0 ORDER BY t2.stageLeftBound DESC LIMIT 1
              ) AS stageIdToUse
              FROM stage AS t1 INNER JOIN stageTaxonConstraint AS t1bis on t1.stageId = t1bis.stageId
               ';
}

#argument of the sub: the DBI Connection object connected to the Bgee database
sub get_stage_equivalences {
    my ($dbh)= @_;
    my %stage_equivalences;
    # Bgee db connection
    my $query = $dbh->prepare(
        get_stage_equivalences_query()
    );

    $query->execute()  or die $query->errstr;
    while ( my @data = $query->fetchrow_array ){
        $stage_equivalences{$data[0]} = $data[1];
    }

    return \%stage_equivalences;
}

sub get_in_between_stages {
    my ($stages, $port, $debug) = @_;
    $debug   = $debug || 0;
    my $host = '127.0.0.1';
    my $done;

    my $server = IO::Socket::INET->new(PeerAddr => $host,
                                       PeerPort => $port,
                                      )  or die"\tUnable to connect to server $host with port $port\n";


    #TODO catch 'Exception' if problems from java socket/daemon
    for (my $i=0; $i<=$#$stages; $i+=2){
        my $start  = @$stages[$i];
        my $end    = @$stages[$i+1];
        my $stages = $start.','.$end;

        if ( !exists $done->{$stages} ){
            if ( $server->connected() ){
                $server->autoflush(1); # Send immediately

                # Send request
#                die "Problem sending [$stages]\n"  if ( send($server, "$stages\n", length("$stages\n"), MSG_NOSIGNAL) != length("$stages\n") );
                print $server "$stages\n";
                $server->flush();

                # Wait for answer
                my $result = <$server>;
                chomp $result;
                $result = ''  if ( $result =~ /Could not find/ || $result =~ /Exception/ || $result =~ /Incorrect number/ );
                $done->{$stages} = $result;
                print "NEW  $stages -> ", $done->{$stages}, "\n"  if ( $debug );
            }
            else {
                die "Socket not available\n";
            }
        }
        else { # exists $done->{$stages}
            print "DONE $stages -> ", $done->{$stages}, "\n"  if ( $debug );
        }
    }
    $server->close();

    return $done;
}

sub get_anatomy_mapping {
    my ($anatIds, $port, $debug) = @_;
    $debug   = $debug || 0;
    my $host = '127.0.0.1';
    my $done;

    my $server = IO::Socket::INET->new(PeerAddr => $host,
                                       PeerPort => $port,
                                      )  or die"\tUnable to connect to server $host with port $port\n";

    for my $anatId ( @$anatIds ){
        next  if ( !defined $anatId || $anatId eq '' );
        if ( !exists $done->{$anatId} ){
            if ( $server->connected() ){
                $server->autoflush(1); # Send immediately

                # Send request
#                die "Problem sending [$anatId]\n"  if ( send($server, "$anatId\n", length("$anatId\n"), MSG_NOSIGNAL) != length("$anatId\n") );
                print $server "$anatId\n";
                $server->flush();

                # Wait for answer
                my $result = <$server>;
                chomp $result;
                $result = ''  if ( $result =~ /Could not find/ || $result =~ /Exception/ || $result =~ /Incorrect number/ );
                $anatId =~ s{^ }{}; # Trick for WormBase anatomy ids that look too short ?!?
                $done->{$anatId} = $result;
                print "NEW  $anatId -> ", $done->{$anatId}, "\n"  if ( $debug );
            }
            else {
                die "Socket not available\n";
            }
        }
        else { #exists $done->{$anatId}
            print "DONE $anatId -> ", $done->{$anatId}, "\n"  if ( $debug );
        }
    }
    $server->close();

    return $done;
}


# Returns all conditions already inserted into the database.
# Returned as a hash where keys are created from anatEntityId/stageId/speciesId/sex/strain information,
# (see sub generate_condition_key), the associated value being a hash with two keys:
# conditionId and exprMappedConditionId.
# $conditions->{anatEntityId/stageId/speciesId/sex/strain}->{'conditionId'}           = conditionId and
# $conditions->{anatEntityId/stageId/speciesId/sex/strain}->{'exprMappedConditionId'} = exprMappedConditionId)
# conditionId: real ID of the condition
# exprMappedConditionId: ID of the corresponding not-too-granular condition,
# to use for insertion into the expression table, if the condition is too granular.
# So, if the condition is not too granular, conditionId = exprMappedConditionId.
# It is thus safe to always use exprMappedConditionId for insertion into the expression/noExpression tables.
# Argument of the method: connector to the database.
sub query_conditions {
    my ($dbh) = @_;

    my $cond = $dbh->prepare('SELECT conditionId, exprMappedConditionId,
            anatEntityId, stageId, speciesId, sex, sexInferred, strain FROM cond');
    $cond->execute()  or die $cond->errstr;
    my $cond_ref = $cond->fetchall_arrayref;

    my $conditions;
    map {
        my $condKey = generate_condition_key($_->[2], $_->[3], $_->[4], $_->[5], $_->[6], $_->[7]);
        $conditions->{ $condKey }->{'conditionId'}           = $_->[0];
        $conditions->{ $condKey }->{'exprMappedConditionId'} = $_->[1];
        $conditions->{ $condKey }->{'strain'}                = $_->[7];
        $conditions->{ $condKey }->{'speciesId'}             = $_->[4];
        $conditions->{ $condKey }->{ 'sexInference' }        = $_->[6];
    }  @{ $cond_ref };

    return $conditions;
}

# Returns all conditions already inserted into the database, and their corresponding not-too-granular exprmappedConditionId.
# Returns a hash with as keys: all conditions already inserted into the database.
# And as values: exprMappedConditionId, the ID of the corresponding not-too-granular condition,
# to use for insertion into the expression table, if the condition is too granular.
# Argument of the method: connector to the database.
## FIXME no, this method should not be used, chip should be queried and grouped per *mapped* condition
sub get_condition_equivalences {
    my ($dbh) = @_;

    my $cond = $dbh->prepare('SELECT conditionId, exprMappedConditionId FROM cond');
    $cond->execute()  or die $cond->errstr;
    my $cond_ref = $cond->fetchall_arrayref;

    my $conditions;
    map {
      $conditions->{ $_->[0] } = $_->[1];
    }  @{ $cond_ref };

    return $conditions;
}


# Load from a file the putative sex related to organs, e.g., an oocyte can be assigned
# either a female or hermaphrodite sex information.
# The file is a TSV file with 5 columns: Uberon ID, Uberon name, female (T/F values),
# male (T/F values), hermaphrodite (T/F values).
# Argument of the sub: path to the file.
# Example returned hash: $anatSexInfo->{'UBERON:00000XXX'} = ('female', 'male')
sub get_anat_sex_info {
    my ($sex_info_file) = @_;

    my $anatSexInfo;
    # TODO use column names rather than column order
    for my $line ( grep { !/^#/ && !/^Uberon ID/} read_file("$sex_info_file", chomp => 1) ){
        #Uberon ID   Uberon name female  male    hermaphrodite
        #AEO:0000013 single-cell tissue  F   F   F
        my @tmp = split(/\t/, $line);

        # basic check on column values
        if ( $tmp[2] ne 'T' && $tmp[2] ne 'F' || $tmp[3] ne 'T' && $tmp[3] ne 'F' ||
                     $tmp[4] ne 'T' && $tmp[4] ne 'F' ) {
            die "Unrecognized file format at line: $line\n";
        }

        $anatSexInfo->{ $tmp[0] } = ();
        # If all columns are false, it means the term is not part of any sex-related branch.
        # Then we assume it exists in all sex types.
        my $allSexFalse = ( $tmp[2] eq 'F' && $tmp[3] eq 'F' && $tmp[4] eq 'F' );
        push @{ $anatSexInfo->{ $tmp[0] } }, $FEMALE_SEX         if ( $tmp[2] eq 'T' || $allSexFalse );
        push @{ $anatSexInfo->{ $tmp[0] } }, $MALE_SEX           if ( $tmp[3] eq 'T' || $allSexFalse );
        push @{ $anatSexInfo->{ $tmp[0] } }, $HERMAPHRODITE_SEX  if ( $tmp[4] eq 'T' || $allSexFalse );
    }

    return $anatSexInfo;
}

# Returns the acceptable sex values for all species in Bgee.
# This allows automatic sex inference, to determine, e.g., that if an annotation is about an oocyte,
# and the species is mouse, then the sex has to be 'female', it can't be 'hermaphrodite'.
# argument of the method: connector to database.
# Example returned hash: $speciesSexInfo->{9606} = ('female', 'male')
sub get_species_sex_info {
    my ($dbh) = @_;

    my $speSex = $dbh->prepare('SELECT speciesId, sex FROM speciesToSex');
    $speSex->execute()  or die $speSex->errstr;
    my $speSex_ref = $speSex->fetchall_arrayref;

    my $speciesSexInfo;
    map {
        if (!exists $speciesSexInfo->{ $_->[0] }) {
            # create new array storing all acceptable sex values for this species
            $speciesSexInfo->{ $_->[0] } = [];
        }
        if ($_->[1] ne $MALE_SEX && $_->[1] ne $FEMALE_SEX && $_->[1] ne $HERMAPHRODITE_SEX) {
            die "Unrecognized sex info: $_->[1] for species: $_->[0]\n";
        }
        push @{ $speciesSexInfo->{ $_->[0] }}, $_->[1];
    }  @{ $speSex_ref };

    return $speciesSexInfo;
}

# Sub to infer sex for a given organ and species.
# Arguments:
# - hash containing acceptable sexes per organ (see sub get_annot_sex_info)
# - hash containing acceptable sexes per species (see sub get_species_sex_info)
# - the ID of the organ to infer sex for
# - the ID of the species to infer sex for
# Returned value: either nothing if no sex could be inferred, or a value corresponding to one of
# either $MALE_SEX, $FEMALE_SEX, or $HERMAPHRODITE_SEX
sub infer_sex {
    my ($anatSexInfo, $speciesSexInfo, $anatEntityId, $speciesId) = @_;

    # If we have no sex info for this organ, we return nothing
    if (!exists $anatSexInfo->{ $anatEntityId }) {
        return;
    }

    # we keep the intersection of sexes acceptable for the considered organ and considered species
    my @sexes = intersect( @{ $anatSexInfo->{ $anatEntityId }}, @{ $speciesSexInfo->{ $speciesId }} );

    # if the intersection contains only one element, then we have a sex inference!
    if (@sexes && scalar(@sexes) == 1) {
        my $inferredSex = $sexes[0];
        # simply check that the value is correct
        if ($inferredSex ne $MALE_SEX && $inferredSex ne $FEMALE_SEX && $inferredSex ne $HERMAPHRODITE_SEX) {
            die "Unrecognized inferred sex: $inferredSex\n";
            return;
        }
        return $inferredSex;
    }
    # Otherwise we could not choose which sex was correct, we return nothing
    return;
}

# Insert the requested condition into database and $conditions hash if not already existing,
# and return the condition in any case. If the condition is too granular for the expression table,
# the corresponding not-too-granular condition will also be created and inserted if not existing.
# This methods needs to be provided with the hashes returned by the subs 'get_stage_equivalences',
# 'get_anat_sex_info', and 'get_species_sex_info'.
sub insert_get_condition {
    my ($dbh, $conditions, $stage_equivalences, $anatEntityId, $stageId, $speciesId, $sex, $strain,
        $anatSexInfo, $speciesSexInfo, $expId, $geneId) = @_;

    $expId  = $expId  || ''; # In case $expId is not provided
    $geneId = $geneId || ''; # In case $expId is not provided

    # ====================================
    # SANITY CHECKS
    # ====================================
    if (!$anatEntityId) {
        die "Missing anatEntityId\n";
    }
    if (!$stageId) {
        die "Missing stageId\n";
    }
    if (!$speciesId) {
        die "Missing speciesId\n";
    }
    if (!$sex) {
        die "Missing sex\n";
    }
    if ( !grep( /^$sex$/, @ACCEPTABLE_SEX_INFO ) ) {
        die "Unrecognized sex info: $sex\n";
    }
    # Check strain vs. "standardized" names, to avoid, e.g., "Wild-Type" vs. "wild-type".
    # We compare strains with only alphanumerical chars, to detect cases such as
    # 'wild-type' vs. 'wild type'.
    my $strainToCheck = $strain;
    $strainToCheck =~ s/[^a-zA-Z0-9]//g;
    my $quotedStrainToCheck = quotemeta($strainToCheck);
    my $quotedStrain = quotemeta($strain);
    if ( grep( /$quotedStrainToCheck/i, map {s/[^a-zA-Z0-9]//g; } @STANDARDIZED_STRAIN_INFO ) &&
        !grep( /^$quotedStrain$/, @STANDARDIZED_STRAIN_INFO ) ) {
        warn "Incorrect standardized strain name detected: '$strain' for $speciesId\n";
    }

    # ====================================
    # SEX INFERENCE
    # ====================================
    # now try to infer sex info related to type of organ, e.g. "ovary" -> female, "testis" -> male
    # Some terms can be potentially assigned different sexes, e.g.,
    # CL:0000023 "oocyte" can belong to both female and hermaphroditic organisms, so we also need
    # some information about the species considered.
    my $sexNotInferred = 0;
    my $sexToUse       = $sex;
    my $sexInference   = $sexNotInferred; # use 0/1 for key generation and database insertion

    my $inferredSex    = infer_sex($anatSexInfo, $speciesSexInfo, $anatEntityId, $speciesId);
    # Check whether we have a discripency between annotated and inferred sex information
    # annotation will have preceedence over inference anyway
    if ( ($sex eq $MALE_SEX || $sex eq $FEMALE_SEX || $sex eq $HERMAPHRODITE_SEX || $sex eq  $MIXED_SEXES) &&
             $inferredSex && $sex ne $inferredSex) {
        warn "Inconsistent sex inference as compared to annotation, annotation will have preceedence. "
             , "Annotation: $anatEntityId - $stageId - $speciesId - $sex - $strain for experiment [$expId] and gene [$geneId]. Inferred sex: $inferredSex\n";
    }
    # XXX: not sure if we should consider 'NA'.
    if ( grep( /^$sex$/, @NO_ANNOT_SEX_INFO ) && $inferredSex) {
        $sexToUse     = $inferredSex;
        $sexInference = 1;
    }


    # ====================================
    # RETRIEVE CONDITION IF EXISTING
    # ====================================
    # Generate a unique condition key for the provided parameters
    my $condKey = generate_condition_key($anatEntityId, $stageId, $speciesId, $sexToUse, $sexInference, $strain);

    # If this condition is already available, nothing to do
    if ( defined $conditions->{$condKey} ){
        #TODO is it a good idea to pass back the modified hash since it is passed as a reference?
        #=> I now return the condition, Sebastien must be informed!
        return ($conditions->{ $condKey }, $conditions);
    }


    # ====================================
    # CREATE NEW CONDITION IF NOT EXISTING
    # ====================================
    # First, if sex was inferred and it is the first time we see this condition, we log a message
    if ( $sexInference ){
        print "Using inferred sex info '$inferredSex' for organ '$anatEntityId' and species '$speciesId' for experiment [$expId] and gene [$geneId]\n";
    }
    # OK, now, retrieve the max conditionId: To avoid to first insert the condition,
    # and then update the table to provide the exprMappedConditionId, we generate our own condition IDs,
    # rather than using the AUTO_INCREMENT of the table.
    # At the same time, we also make a check on the strains aready stored in existing conditions,
    # in the same function, to avoid iterating the hash several times.
    # XXX: WARNING, THIS IS NOT THREAD SAFE!
    # (well, worse case scenario, a query will fail if the ID returned was used meanwhile)
    # XXX: maybe there is something better to do than to iterate all conditions each time?
    my $condId = get_max_condition_id_and_check_strain($conditions, $strain) + 1;

    # REMAPPING TO NOT-TOO-GRANULAR CONDITION:
    #
    # * Now, we check whether the requested condition is too granular to be used in the expression table;
    # in that case, we need to retrieve or create the corresponding not-too-granular condition.
    #
    # * Also, we don't want to distinguish sex inference types for insertion into the expression table,
    # so we always use a value of '0' for $sexInference for insertion into expression table:
    # conditions with sex inferred are considered too granular, kind of.
    # Similarly, we do not want to distinguish between sex 'not annotated' and 'NA',
    # so we always map to 'not annotated'.
    #
    # * And also, we don't want to distinguish some strain types for insertion into expression table:
    # 'NA', 'not annotated', 'wild-type', 'confidential_restricted_data' are all mapped to
    # 'wild-type' in expression table (since everything is wild-type anyway in Bgee so far).

    my $mappedSexToUse = $sexToUse;
    if ( $sexToUse eq $NA_SEX ){
        $mappedSexToUse = $NOT_ANNOTATED_SEX;
    }
    my $mappedStrainToUse = $strain;
    if ( $strain eq $NA_STRAIN || $strain eq $NOT_ANNOTATED_STRAIN || $strain eq $CRD_STRAIN ){
        $mappedStrainToUse = $WILD_TYPE_STRAIN;
    }

    my $exprMappedCondKey = generate_condition_key($anatEntityId, $stage_equivalences->{ $stageId },
            $speciesId, $mappedSexToUse, $sexNotInferred, $mappedStrainToUse);
    my $exprMappedCondId = $condId;

    # condition too granular or with sex inferred or NA
    if ( $condKey ne $exprMappedCondKey ){
        # Does the condition, not-too-granular and with no sex inference and no NA, already exist?
        if ( defined $conditions->{$exprMappedCondKey} ){
            $exprMappedCondId = $conditions->{$exprMappedCondKey}->{'conditionId'};
        } else {
            # Does not exist, update the condition IDs and insert into database
            $exprMappedCondId = $condId;
            $condId = $exprMappedCondId + 1;
            # Not-too-granular conditions are mapped to themselves
            insert_condition($dbh, $exprMappedCondId, $exprMappedCondId,
                    $anatEntityId, $stage_equivalences->{ $stageId }, $speciesId,
                    $mappedSexToUse, $sexNotInferred, $mappedStrainToUse);

            # And update the $condition hash
            $conditions->{ $exprMappedCondKey }->{ 'conditionId' }           = $exprMappedCondId;
            $conditions->{ $exprMappedCondKey }->{ 'exprMappedConditionId' } = $exprMappedCondId;
            $conditions->{ $exprMappedCondKey }->{ 'strain' }                = $mappedStrainToUse;
            $conditions->{ $exprMappedCondKey }->{ 'speciesId' }             = $speciesId;
            $conditions->{ $exprMappedCondKey }->{ 'sexInference' }          = $sexNotInferred;
        }
    }
    # Assertion test: at this point, if sex was inferred or was equal to NA, we should always have a "mapped" condition
    if ( ($sexInference || $sexToUse eq $NA_SEX) && $condId == $exprMappedCondId ){
        die "Assertion error, inferred or NA sex conditions should be seen as too granular\n";
    }
    # Same for strains mapped to 'wild-type'
    if ( ($strain eq $NA_STRAIN || $strain eq $NOT_ANNOTATED_STRAIN || $strain eq $CRD_STRAIN) &&
          $condId == $exprMappedCondId ){
        die "Assertion error, strain info '$strain' should be seen as too granular\n";
    }

    # Now, we can insert the condition itself
    insert_condition($dbh, $condId, $exprMappedCondId, $anatEntityId, $stageId, $speciesId, $sexToUse, $sexInference, $strain);
    # And update the $condition hash
    $conditions->{ $condKey }->{ 'conditionId' }           = $condId;
    $conditions->{ $condKey }->{ 'exprMappedConditionId' } = $exprMappedCondId;
    $conditions->{ $condKey }->{ 'strain' }                = $strain;
    $conditions->{ $condKey }->{ 'speciesId' }             = $speciesId;
    $conditions->{ $condKey }->{ 'sexInference' }          = $sexInference;

    # ====================================
    # RETURN CONDITION
    # ====================================
    # returns "real" condition Id (may not be used in expression table) and modified hash of conditions
    #TODO is it a good idea to pass back the modified hash since it is passed as a reference?
    #=> I now return the condition, Sebastien must be informed!
    return ($conditions->{ $condKey }, $conditions);
}

# Basic sub to insert condition into database, with already generated conditionId
# (no use of AUTO_INCREMENT)
# If $sex or $sexInferred do not correspond to allowed values, die.
sub insert_condition {
    my ($dbh, $conditionId, $exprMappedConditionId, $anatEntityId, $stageId, $speciesId, $sex, $sexInferred, $strain) = @_;
    if ( !grep( /^$sex$/, @ACCEPTABLE_SEX_INFO ) ) {
        die "Incorrect sex value: $sex\n";
    }
    if ($sexInferred != 0 && $sexInferred != 1) {
        die "Incorrect sexInferred value: $sexInferred\n";
    }

    my $ins = $dbh->prepare('INSERT INTO cond (conditionId, exprMappedConditionId,
            anatEntityId, stageId, speciesId, sex, sexInferred, strain) VALUES (?, ?, ?, ?, ?, ?, ?, ?)');
    $ins->execute($conditionId, $exprMappedConditionId, $anatEntityId, $stageId, $speciesId, $sex, $sexInferred, $strain)
            or die $ins->errstr;
}

# Generate a key from anatEntityId/stageId/speciesId/sex/sexInferenceType/strain information.
# If $sex or $sexInferred do not correspond to allowed values, die.
sub generate_condition_key {
    my ($anatEntityId, $stageId, $speciesId, $sex, $sexInferred, $strain) = @_;
    if ( !grep( /^$sex$/, @ACCEPTABLE_SEX_INFO ) ) {
        die "Incorrect sex value: $sex\n";
    }
    if ($sexInferred != 0 && $sexInferred != 1) {
        die "Incorrect sexInferred value: $sexInferred\n";
    }

    return $anatEntityId.'--'.$stageId.'--'.$speciesId.'--'.$sex.'--'.$sexInferred.'--'.$strain;
}

# Retrieve the max condition ID used from the hash of conditions.
# First max ID if $conditions is empty is 0.
# At the same time, we take the opportunity to verify the strain used
# as compared to already stored strains. Yeah, it's not very related, I know...
# TODO I guess there is a better way to write it
sub get_max_condition_id_and_check_strain {
    my ($conditions, $strain) = @_;
    my $maxId = 0;
    my @warnDone = ();
    my $strainToCheck = $strain;
    $strainToCheck =~ s/[^a-zA-Z0-9]//g;

    for my $key ( keys %{ $conditions } ){
        if( $conditions->{ $key }->{ 'conditionId' } > $maxId ){
            $maxId = $conditions->{ $key }->{ 'conditionId' };
        }

        my $quotedCondStrain  = quotemeta($conditions->{ $key }->{ 'strain' });
        my $condStrainToCheck = $conditions->{ $key }->{ 'strain' };
        $condStrainToCheck    =~ s/[^a-zA-Z0-9]//g;
        if( !grep( /^$quotedCondStrain$/, @warnDone) &&
             lc $condStrainToCheck eq lc $strainToCheck &&
             $conditions->{ $key }->{ 'strain' } ne $strain ) {

            warn "Different strain names equal in lower case for alphanumerical chars: ".
                 "'$conditions->{ $key }->{ 'strain' }' and '$strain'. ".
                 "Maybe you should merge these conditions before inserting data into the expression tables.\n";
            push @warnDone, $conditions->{ $key }->{ 'strain' };
        }
    }
    return $maxId
}

# Get all bgeeGeneId & geneId
sub query_bgeeGene {
    my ($dbh, $speciesId) = @_;

    my $gene = $dbh->prepare('SELECT bgeeGeneId, geneId FROM gene WHERE speciesId = ?');
    if ( !$speciesId || $speciesId eq '' ){
        $gene = $dbh->prepare('SELECT bgeeGeneId, geneId FROM gene');
        $gene->execute()  or die $gene->errstr;
    }
    else {
        $gene->execute($speciesId)  or die $gene->errstr;
    }
    my $gene_ref = $gene->fetchall_arrayref;

    my $genes;
    map { $genes->{ $_->[1] } = $_->[0]; }  @{ $gene_ref };

    return $genes;
}


# An implementation of fractional ranking,
# see https://en.wikipedia.org/wiki/Ranking#Fractional_ranking_.28.221_2.5_2.5_4.22_ranking.29
# The data is assumed to be passed in an array of hashref {"id" => SOME_ID, "val" => VALUE}
# the key being the ID (of whatever is ranked) and the value is the value used for ranking.
#
# /!\ WARNING: INPUT DATA MUST BE SORTED /!\
# As the implementation assign the same rank ex-aequos that have the same value,
# but does not check that the ordering is consistent. Therefore, it works whether higher values have higher ranks or
# the inverse.
sub fractionnal_ranking  {
    my @data = @_;
    my $rank = 1;
    my $last = -1;
    my $first = 1;
    my %result;
    my @exaquo;
    my $count = 0;
    my $excount = 0;

    foreach (@data) {
        my $id = $_->{"id"};
        my $value = $_->{"val"};
        if ($value != $last && $first == 0) {
            $rank = ($count + 1) + 0.5 * (scalar( @exaquo ) - 1);
            foreach my $lid (@exaquo) {
                $result{$lid} = $rank;
            }
            $count += scalar( @exaquo );

            @exaquo = ();
            $excount = 0;
        }
        $last = $value;
        $exaquo[$excount++] = $id;
        $first = 0;
    }
    $rank = ($count + 1) + 0.5 * (scalar( @exaquo ) - 1);
    foreach my $lid (@exaquo) {
        $result{$lid} = $rank;
    }
    $count += scalar( @exaquo );
    return %result;
}

# An implementation of dense ranking,
# see https://en.wikipedia.org/wiki/Ranking#Dense_ranking_.28.221223.22_ranking.29
# The data is assumed to be passed in an array of hashref {"id" => SOME_ID, "val" => VALUE}
# the key being the ID (of whatever is ranked) and the value is the value used for ranking.
#
# /!\ WARNING: INPUT DATA MUST BE SORTED /!\
# As the implementation assign the same rank ex-aequos that have the same value,
# but does not check that the ordering is consistent. Therefore, it works whether higher values have higher ranks or
# the inverse.
sub dense_ranking  {
    my @data = @_;
    my $rank = 0;
    my $last = -1;
    my $first = 1;
    my %result;
    my @exaquo;
    my $count = 0;
    my $excount = 0;

    foreach (@data) {
        my $id = $_->{"id"};
        my $value = $_->{"val"};
        if ($value != $last || $first == 1) {
            $rank++;
        }
        $result{$id} = $rank;
        $last = $value;
        $first = 0;
    }
    return %result;
}

# Reverse a hash to a hash of arrays.  The values of the result hash are
# references to arrays containing all keys of the original hash having the same
# value (which becomes the corresponding key in the new hash).
# -> Original hash
# <- Resulting hash of arrays
# From http://www.volkerschatz.com/perl/snippets/revhash.html
sub revhash(\%)
{
    my ($orighash)= @_;
    my %result;

    for my $key (keys %$orighash) {
        push @{$result{$$orighash{$key}}}, $key;
    }
    return ( %result );
}


# Return hash of stageId contained in Bgee database
# with leftBound and rightBound stage associated
sub getBgeedbStages {
    my ($dbh) = @_;

    my %stages;
    my $selStage = $dbh->prepare('SELECT stageId, stageLeftBound, stageRightBound FROM stage');
    $selStage->execute()  or die $selStage->errstr;
    while ( my @data = $selStage->fetchrow_array ) {
        $stages{$data[0]}->{'leftBound'}  = $data[1];
        $stages{$data[0]}->{'rightBound'} = $data[2];
    }
    $selStage->finish;

    return \%stages;
}

# Return hash of anatEntityId contained in Bgee database
# with start and end stages associated
sub getBgeedbOrgans {
    my ($dbh) = @_;

    my %organs;
    my $selAnat = $dbh->prepare('SELECT anatEntityId, startStageId, endStageId FROM anatEntity');
    $selAnat->execute()  or die $selAnat->errstr;
    while ( my @data = $selAnat->fetchrow_array ) {
        $organs{$data[0]}->{'startStageId'} = $data[1];
        $organs{$data[0]}->{'endStageId'}   = $data[2];
    }
    $selAnat->finish;

    return \%organs;
}

1;

