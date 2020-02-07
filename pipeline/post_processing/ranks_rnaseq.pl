#!/usr/bin/env perl
use strict;
use warnings;
use Time::HiRes qw( time );
use diagnostics;

# Updates the ranks in rnaSeqResult table and RNA-Seq mean rank in the expression table.
# Philippe Moret, created Oct 2015.
# Frederic Bastian, updated June 2016.
# Frederic Bastian, updated Feb. 2017: adapt to new conditions and new schema in Bgee 14
# Frederic Bastian, updated Jan. 2020: parallelize rank computations

use Parallel::ForkManager;
my $parallel_jobs  = 15;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
use Getopt::Long;

$|=1;

# Define arguments and their default value
my ($bgee_connector) = ('');
my ($ranks_computed) = (0);
my %opts = ('bgee=s'         => \$bgee_connector, # Bgee connector string
            'ranks_computed' => \$ranks_computed,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0 -bgee=\$(BGEECMD)
\t-bgee             Bgee connector string
\t-ranks_computed   Skip generation of raw ranks per library
\n";
    exit 1;
}

my $dbh = Utils::connect_bgee_db($bgee_connector);

#Set to 0 in order to disable autocommit to optimize speed
my $auto = 0;

if ( $auto == 0 ){
    $dbh->{'AutoCommit'} = 0;
}

# Reasonning of the computations:
# 1) identify the valid set of genes that should be considered for ranking in all libraries:
# the set of all genes that received at least one read over all libraries in Bgee.
# 2) compute gene fractional ranks in table rnaSeqResult, for each RNA-Seq library, based on TPM values
# 3) update expression table: compute weighted mean of ranks per gene and condition,
# weighted by the number of distinct ranks in each library: we assume that libraries with higher
# number of distinct ranks have a higher power for ranking genes.
# Note that we do not "normalize" ranks between samples before computing the mean, as for Affymetrix data:
# all libraries are used to produce ranking over always the same set of genes in a given species,
# so the genomic coverage is always the same, and no "normalization" is required. The higher power
# at ranking genes of a library (for instance, thanks to a higher number of mapped reads)
# is taken into account by weighting the mean by the number of distinct ranks in the library,
# not by "normalizing" away libraries with lower max ranks; this would penalize conditions
# with a lower number of expressed genes, and thus with more ex-aequo ranked genes, corresponding
# to genes receiving 0 read.
# 4) also, insert max ranks per mapped condition in condition table
# (will allow to normalize ranks between conditions and data types), and sum of distinct rank counts
# per gene and condition, in expression table (used to compute weigthed mean over all data types in a condition)


##############################################
# IDENTIFY VALID GENES                       #
##############################################
# We rank all genes that have received at least one read in any condition.
# So we always rank the same set of gene in a given species over all libraries.
# We don't use a temp table to be able to close/open the connection for each condition parameter combination.
# So we have to drop the table at the end.
my $dropValidGenesStmt = $dbh->prepare('DROP TABLE IF EXISTS rnaSeqValidGenes');
my $validGenesStmt     = $dbh->prepare('CREATE TABLE rnaSeqValidGenes
                                            (PRIMARY KEY(bgeeGeneId))
                                            SELECT DISTINCT t1.bgeeGeneId
                                            FROM rnaSeqResult AS t1
                                            WHERE t1.reasonForExclusion = "'.$Utils::CALL_NOT_EXCLUDED.'"');

printf('Identifying set of valid genes for ranking: ');
$dropValidGenesStmt->execute()  or die $dropValidGenesStmt->errstr;
$validGenesStmt->execute()  or die $validGenesStmt->errstr;
printf("Done\n");


##############################################
# COMPUTE RANKS PER LIBRARY                  #
##############################################
sub compute_rank_lib_batch {
    my ($libBatchRef) = @_;
    my $batchLength = scalar @{ $libBatchRef };

    # Connection to database in the parent process must have been closed
    # before calling this sub, otherwise it will generate errors
    # (ForkManager and DBI don't go well together, see https://www.perlmonks.org/?node_id=752289)
    my $pm = new Parallel::ForkManager($parallel_jobs);
    # Returning ranks computed from each thread to the main thread for update
    # (reading and updating the table at the same time leads to deadlock issues).
    # We will put everything in memory
    my @batchLibRankInfo = ();  # for collecting responses
    $pm -> run_on_finish ( # called BEFORE the first call to start()
        sub {
            my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;

            # retrieve data structure from child
            # $data_structure_reference is a reference to a hash.
            # Key "libraryId" associated to a value being the ID of a RNA-Seq library
            # Key "ranks" associated to a value being a reference to a hash,
            # with ranks as keys, and reference to an array of gene IDs with that rank as value
            if (defined($data_structure_reference)) {
                push @batchLibRankInfo, $data_structure_reference;
            } else {  # problems occurring during storage or retrieval will throw a warning
                warn "No message received from child process $pid!\n";
            }
        }
    );

    # Compute ranks
    for my $k ( 0..$batchLength-1 ) {

        # Forks and returns the pid for the child
        my $pid = $pm->start and next;
        # Get a new database connection for each thread
        my $dbh_thread = Utils::connect_bgee_db($bgee_connector);
        #Set to 0 in order to disable autocommit to optimize speed
        if ( $auto == 0 ){
            $dbh_thread->{'AutoCommit'} = 0;
        }
        my $rnaSeqLibraryId = ${$libBatchRef}[$k];
        my $rnaSeqResultsStmt = $dbh_thread->prepare(
            'SELECT DISTINCT t1.bgeeGeneId, t1.tpm
             FROM rnaSeqResult AS t1 '.
             # join to table rnaSeqValidGenes to force
             # the selection of valid genes
             'INNER JOIN rnaSeqValidGenes AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId
             WHERE t1.rnaSeqLibraryId = ?
             AND t1.reasonForExclusion = "'.$Utils::CALL_NOT_EXCLUDED.'"
             ORDER BY t1.tpm DESC');

        $rnaSeqResultsStmt->execute($rnaSeqLibraryId)  or die $rnaSeqResultsStmt->errstr;

        my @results = map { {'id' => $_->[0], 'val' => $_->[1]} } @{$rnaSeqResultsStmt->fetchall_arrayref};
        $dbh_thread->disconnect();
        my %sorted = Utils::fractionnal_ranking(@results);
        # we get ranks as keys, with reference to an array of gene IDs with that rank as value
        my %reverseHash = Utils::revhash(%sorted);
        # We add the info of library ID in the returned hash to be used by the main thread
        my %returnedInfo;
        $returnedInfo{'libraryId'} = $rnaSeqLibraryId;
        $returnedInfo{'ranks'} = \%reverseHash;
        #print status
        printf("Lib: %s - %d/%d", $rnaSeqLibraryId, $k+1, $batchLength);
        # Terminates the child process and send the results to the main thread
        $pm->finish(0, \%returnedInfo);
    }
    $pm->wait_all_children;

    return \@batchLibRankInfo;
}

# We are going to update the ranks in batches to not overload memory,
# this sub will be called every x child processes terminate
sub update_ranks {
    my ($batchLibRankInfoRef, $done) = @_;

    # Reopen the connection in parent process that we disconnected before launching
    # the child processes
    $dbh_sub = Utils::connect_bgee_db($bgee_connector);
    if ( $auto == 0 ){
        $dbh_sub->{'AutoCommit'} = 0;
    }
    # if several genes at a same rank, we'll update them at once
    # with a 'bgeeGeneId IN (?,?, ...)' clause.
    # If only one gene at a given rank, updated with the prepared statement below.
    my $rankUpdateStart = 'UPDATE rnaSeqResult SET rank = ? WHERE rnaSeqLibraryId = ? and bgeeGeneId ';
    my $rnaSeqResultUpdateStmt = $dbh_sub->prepare($rankUpdateStart.'= ?');

    for my $libRankInfoRef (@{ $batchLibRankInfoRef }) {
        my $rnaSeqLibraryId = $libRankInfoRef->{'libraryId'};
        my $libRankRef      = $libRankInfoRef->{'ranks'};
        for my $rank ( keys %{ $libRankRef } ){
            my $geneIds_arrRef = $libRankRef->{$rank};
            my @geneIds_arr = @$geneIds_arrRef;
            my $geneCount = scalar @geneIds_arr;
            if ( $geneCount == 1 ){
                my $geneId = $geneIds_arr[0];
                $rnaSeqResultUpdateStmt->execute($rank, $rnaSeqLibraryId, $geneId)
                    or die $rnaSeqResultUpdateStmt->errstr;
            } else {
                my $query = $rankUpdateStart.'IN (';
                for ( my $i = 0; $i < $geneCount; $i++ ){
                    if ( $i > 0 ){
                        $query .= ', ';
                    }
                    $query .= '?';
                }
                $query .= ')';
                my $rnaSeqRankMultiUpdateStmt = $dbh_sub->prepare($query);
                $rnaSeqRankMultiUpdateStmt->execute($rank, $rnaSeqLibraryId, @geneIds_arr)
                    or die $rnaSeqRankMultiUpdateStmt->errstr;
            }
        }

        printf("Lib updated: %s\tIdx:%d", $rnaSeqLibraryId, $done);
        $done += 1;
    }
    # commit here if auto-commit is disabled
    if ( $auto == 0 ){
        $dbh_sub->commit()  or die('Failed commit');
    }
    # Important to disconnect the DBI connection open in parent process,
    # otherwise it will generate errors
    # (ForkManager and DBI don't go well together, see https://www.perlmonks.org/?node_id=752289)
    $dbh_sub->disconnect();
}



if ( !$ranks_computed ) {

    # Clean potentially already computed ranks
#    my $cleanRNASeq = $dbh->prepare("UPDATE rnaSeqResult SET rank = NULL");
#    my $cleanLib    = $dbh->prepare("UPDATE rnaSeqLibrary SET libraryMaxRank = NULL,
#                                                              libraryDistinctRankCount = NULL");
#    printf("Cleaning existing data: ");
#    $cleanRNASeq->execute() or die $cleanRNASeq->errstr;
#    $cleanLib->execute() or die $cleanLib->errstr;
#    printf("Done\n");

    # Queries to compute gene ranks per library.
    # We rank all genes that have received at least one read in any condition.
    # So we always rank the same set of gene in a given species over all libraries.
    # We assume that each library maps to only one species through its contained genes.
    my $rnaSeqLibStmt = $dbh->prepare('SELECT t1.rnaSeqLibraryId FROM rnaSeqLibrary AS t1
                                       WHERE EXISTS (SELECT 1 FROM rnaSeqResult AS t2
                                       WHERE t1.rnaSeqLibraryId = t2.rnaSeqLibraryId
                                       AND t2.readsCount > 0)');

    my $t0 = time();
    # Get the list of all rna-seq libraries
    $rnaSeqLibStmt->execute()  or die $rnaSeqLibStmt->errstr;

    my @libs = map { $_->[0] } @{$rnaSeqLibStmt->fetchall_arrayref};
    my $l = @libs;
    printf("Found %d libraries\n", $l);


    print("Rank computation per library...\n");
    # Disconnect the DBI connection open in parent process, otherwise it will generate errors
    # (ForkManager and DBI don't go well together, see https://www.perlmonks.org/?node_id=752289)
    $dbh->disconnect();
    # We are going to compute/store ranks using parallelization,
    # by batch of libraries not to overload memory
    my $libCountBeforeUpdate = 1000;
    my $count = 0;
    while ( my @next_libs = splice(@libs, 0, $libCountBeforeUpdate) ) {
        print($count*$libCountBeforeUpdate.' done/'.$l.' libraries, batch of '.scalar(@next_libs)
            ." libraries to analyze\n")
        my $rankBatchRef = compute_rank_lib_batch(\@next_libs);
        update_ranks($rankBatchRef, $count*$libCountBeforeUpdate);
        $count += 1;
    }
    print("Rank computation per library done\n");


    # ##############
    # Store max rank and number of distinct ranks per library.
    # Reopen connection that was closed
    $dbh = Utils::connect_bgee_db($bgee_connector);
    if ( $auto == 0 ){
        $dbh->{'AutoCommit'} = 0;
    }
    my $sql =
    "UPDATE rnaSeqLibrary AS t0
     INNER JOIN (
         SELECT t1.rnaSeqLibraryId, MAX(t1.rank) AS maxRank, COUNT(DISTINCT t1.rank) AS distinctRankCount
         FROM rnaSeqResult AS t1
         WHERE t1.reasonForExclusion NOT IN ('$Utils::EXCLUDED_FOR_PRE_FILTERED', '$Utils::EXCLUDED_FOR_UNDEFINED')
         GROUP BY t1.rnaSeqLibraryId
     ) AS ranks ON t0.rnaSeqLibraryId = ranks.rnaSeqLibraryId
     SET t0.libraryMaxRank = ranks.maxRank, t0.libraryDistinctRankCount = ranks.distinctRankCount
     WHERE EXISTS (
         SELECT 1 FROM rnaSeqResult AS t2
         WHERE t2.expressionId IS NOT NULL AND t2.rnaSeqLibraryId = t0.rnaSeqLibraryId
     )";

    $t0 = time();
    printf('Inserting max ranks and distinct rank counts in rnaSeqLibrary table...');
    my $maxRankLibStmt = $dbh->prepare($sql);
    $maxRankLibStmt->execute()  or die $maxRankLibStmt->errstr;
    printf("Done in %.2fs\n", (time() - $t0));


    if ( $auto == 0 ){
        $dbh->commit()  or die('Failed commit');
    }
}


####################################################################################
# COMPUTE WEIGHTED MEAN RANKS PER GENE-CONDITION #
####################################################################################

# Store max rank in each species, for later normalization between conditions, data types and species
my $dropMaxRankSpeciesStmt = $dbh->prepare('DROP TABLE IF EXISTS rnaSeqMaxSpecies');
my $maxRankSpeciesStmt  = $dbh->prepare('
CREATE TABLE rnaSeqMaxSpecies (PRIMARY KEY(speciesId))
    SELECT t2.speciesId, MAX(t1.libraryMaxRank) AS maxRank
    FROM rnaSeqLibrary AS t1
    INNER JOIN cond AS t2 ON t1.conditionId = t2.conditionId
    GROUP BY t2.speciesId');

my $t0 = time();
printf('Create temporary table with max rank per species: ');
$dropMaxRankSpeciesStmt->execute()  or die $dropMaxRankSpeciesStmt->errstr;
$maxRankSpeciesStmt->execute()      or die $maxRankSpeciesStmt->errstr;
printf("OK in %.2fs\n", (time() - $t0));
if ( $auto == 0 ){
    $dbh->commit()  or die('Failed commit');
}


# Retrieve codition parameter combinations to compute RNA-Seq ranks for.
my $condParamCombinationsArrRef = Utils::get_cond_param_combinations($dbh, $Utils::RNA_SEQ_DATA_TYPE);
print "All condition parameter combinations to compute for RNA-Seq data:\n";
print "@$_\n"  for @{$condParamCombinationsArrRef};

# connection will be opened/closed at each combination iteration to delete all temp tables
$dbh->disconnect();


for my $condParamCombArrRef ( @{$condParamCombinationsArrRef} ){
    my @condParamComb = @{$condParamCombArrRef};
    print "***** combination @condParamComb *****\n";

    $dbh = Utils::connect_bgee_db($bgee_connector);
    if ( $auto == 0 ){
        $dbh->{'AutoCommit'} = 0;
    }
    # Queries to first clean data
#    my $cleanExpr = $dbh->prepare("UPDATE globalExpression AS t1
#                                   INNER JOIN globalCond AS t2 ON t1.globalConditionId = t2.globalConditionId
#                                   SET rnaSeqMeanRank              = null,
#                                       rnaSeqMeanRankNorm          = null,
#                                       rnaSeqDistinctRankSum       = null,
#                                       rnaSeqGlobalMeanRank        = null,
#                                       rnaSeqGlobalMeanRankNorm    = null,
#                                       rnaSeqGlobalDistinctRankSum = null
#                                   WHERE ".Utils::get_cond_param_comb_sql_clause($condParamCombArrRef, "t2"));
#    my $cleanCond = $dbh->prepare("UPDATE globalCond SET rnaSeqMaxRank = null,
#                                                     rnaSeqGlobalMaxRank = null
#                                   WHERE ".Utils::get_cond_param_comb_sql_clause($condParamCombArrRef, "globalCond"));
#
#    printf("Cleaning existing data: ");
#    $cleanExpr->execute() or die $cleanExpr->errstr;
#    $cleanCond->execute() or die $cleanCond->errstr;
#    if ($auto == 0) {
#        $dbh->commit() or die("Failed commit");
#    }
#    printf("Done\n");


    # connection will be opened/closed at each combination iteration to delete all temp tables
    $dbh->disconnect();

    # we will compute two different rank information for each combination: one taking into account
    # all rank info mapped to a given condition ('self'), or all rank info mapped to a given condition
    # plus all its sub-conditions ('global').
    my $selfRanks   = 0;
    my $globalRanks = 0;

    while ( !$selfRanks || !$globalRanks ){
        if ( !$selfRanks ){
            print "** Computation of self ranks\n"
        } else {
            print "** Computation of global ranks\n"
        }
        $dbh = Utils::connect_bgee_db($bgee_connector);
        if ( $auto == 0 ){
            $dbh->{'AutoCommit'} = 0;
        }

        # Store an association between each globalCondition and the libraries considered in it
        my $sql =
        'CREATE TEMPORARY TABLE globalCondToLib (
             PRIMARY KEY(rnaSeqLibraryId, globalConditionId), INDEX(globalConditionId))
             SELECT DISTINCT t1.globalConditionId, t4.rnaSeqLibraryId '.
             # Retrieve the valid raw conditions mapped to each globalCondition
             'FROM globalCond AS t1
              INNER JOIN globalCondToCond AS t2 ON t1.globalConditionId = t2.globalConditionId ';
        if ( !$selfRanks ){
            $sql .= "AND t2.conditionRelationOrigin = 'self' ";
        } else {
            $sql .= "AND t2.conditionRelationOrigin IN ('self', 'descendant') ";
        }
        $sql .= 'INNER JOIN cond AS t3 ON t3.exprMappedConditionId = t2.conditionId '.
                # retrieve the libraries present in this globalCondition
                'INNER JOIN rnaSeqLibrary AS t4 ON t3.conditionId = t4.conditionId'.
                # Use only globalConditions for the requested condition parameter combination
                ' WHERE '.Utils::get_cond_param_comb_sql_clause($condParamCombArrRef, 't1');

        $t0 = time();
        printf('Creating temp table mapping RNA-Seq libraries to globalConditions: ');
        my $libToGlobalCondStmt = $dbh->prepare($sql);
        $libToGlobalCondStmt->execute()  or die $libToGlobalCondStmt->errstr;
        printf("Done in %.2fs\n", (time() - $t0));


        # No within-datatype normalization, because we consider that all libraries in each species
        # could access a same putative max rank, and because we don't normalize between species at this point,
        # we're only supposed to normalize samples in a same condition and species (for Affymetrix data).


        # ###################
        # Run computations per globalCondition
        # ###################

#       Prepare queries
        my $tableName = 'weightedMeanRank';
        # compute weighted mean normalized ranks, and sum of numbers of distinct ranks
        $sql = 'CREATE TEMPORARY TABLE '.$tableName.'
                SELECT STRAIGHT_JOIN
                      rnaSeqResult.bgeeGeneId,
                      SUM(rnaSeqResult.rank * rnaSeqLibrary.libraryDistinctRankCount)
                          /SUM(rnaSeqLibrary.libraryDistinctRankCount) AS meanRank,
                      SUM(rnaSeqLibrary.libraryDistinctRankCount) AS distinctRankCountSum
                FROM globalCondToLib
                INNER JOIN rnaSeqLibrary ON rnaSeqLibrary.rnaSeqLibraryId = globalCondToLib.rnaSeqLibraryId
                INNER JOIN rnaSeqResult ON rnaSeqResult.rnaSeqLibraryId = rnaSeqLibrary.rnaSeqLibraryId
                INNER JOIN rnaSeqValidGenes ON rnaSeqResult.bgeeGeneId = rnaSeqValidGenes.bgeeGeneId
                WHERE rnaSeqResult.expressionId IS NOT NULL AND globalCondToLib.globalConditionId = ?
                GROUP BY rnaSeqResult.bgeeGeneId';
        my $rnaSeqWeightedMeanStmt  = $dbh->prepare($sql);
        my $dropLibWeightedMeanStmt = $dbh->prepare('DROP TABLE '.$tableName);

        #update the expression table
        $sql = 'UPDATE '.$tableName.'
                STRAIGHT_JOIN globalExpression '.
                # we build the query this way in order to benefit from the clustered index
                # on (bgeeGeneId, globalConditionId) of the globalExpression table
               'ON globalExpression.bgeeGeneId = '.$tableName.'.bgeeGeneId
                    AND globalExpression.globalConditionId = ? ';
        if ( !$selfRanks ){
            $sql .= 'SET globalExpression.rnaSeqMeanRank = '.$tableName.'.meanRank,
                globalExpression.rnaSeqDistinctRankSum = '.$tableName.'.distinctRankCountSum';
        } else {
            $sql .= 'SET globalExpression.rnaSeqGlobalMeanRank = '.$tableName.'.meanRank,
                globalExpression.rnaSeqGlobalDistinctRankSum = '.$tableName.'.distinctRankCountSum';
        }
        my $expressionUpdateMeanRank = $dbh->prepare($sql);

        # Retrieve all global conditions
        my $queryConditions = $dbh->prepare('SELECT DISTINCT globalConditionId FROM globalCondToLib');
        $queryConditions->execute()  or die $queryConditions->errstr;
        my @conditions = ();
        while ( my @data = $queryConditions->fetchrow_array ){
            push(@conditions, $data[0]);
        }
        my $conditionCount = scalar(@conditions);
        print "Computing ranks per global condition, $conditionCount conditions retrieved.\n";

        my $i = 0;
        for my $globalConditionId ( @conditions ){
            $i++;
            $t0 = time();

#            printf("Creating temp table for weighted mean ranks: ");
            $rnaSeqWeightedMeanStmt->execute($globalConditionId)  or die $rnaSeqWeightedMeanStmt->errstr;
#            printf("OK in %.2fs\n", (time() - $t0));

#            $t0 = time();
#            printf("Updating expression table with normalized mean ranks and sum of numbers of distinct ranks: ");
            $expressionUpdateMeanRank->execute($globalConditionId)  or die $expressionUpdateMeanRank->errstr;
#            printf("OK in %.2fs\n", (time() - $t0));

            $dropLibWeightedMeanStmt->execute()  or die $dropLibWeightedMeanStmt->errstr;

            # commit here if auto-commit is disabled
            if ( $auto == 0 ){
                $dbh->commit()  or die('Failed commit');
            }
            if ( ($i / 100 - int($i / 100)) == 0 ){
                printf("$i conditions done.\n");
            }
        }

        # update the globalCond table
        # We apply the same max rank to all expression calls of a same species, in all gene-condition
        $sql = 'UPDATE globalCond '.
                # Retrieve the conditions with RNA-Seq data
                'INNER JOIN globalCondToLib ON globalCond.globalConditionId = globalCondToLib.globalConditionId '.
                # get the max of max ranks per species
                'INNER JOIN rnaSeqMaxSpecies ON globalCond.speciesId = rnaSeqMaxSpecies.speciesId ';
        if ( !$selfRanks ){
            $sql .= 'SET globalCond.rnaSeqMaxRank = rnaSeqMaxSpecies.maxRank';
        } else {
            $sql .= 'SET globalCond.rnaSeqGlobalMaxRank = rnaSeqMaxSpecies.maxRank';
        }
        $t0 = time();
        printf('Updating condition table with max ranks: ');
        my $updateExpressionMaxRank = $dbh->prepare($sql);
        $updateExpressionMaxRank->execute()  or die $updateExpressionMaxRank->errstr;
        printf("OK in %.2fs\n", (time() - $t0));

        if ( $auto == 0 ){
            $dbh->commit()  or die('Failed commit');
        }

        # closing the connection will destroy the temporary tables
        $dbh->disconnect();

        if ( !$selfRanks ){
            $selfRanks = 1;
        } else {
            $globalRanks = 1;
        }
    }
}


# reopen the connection that was closed at the last parameter combination iteration,
# to drop non-temporary tables
$dbh = Utils::connect_bgee_db($bgee_connector);
if ( $auto == 0 ){
    $dbh->{'AutoCommit'} = 0;
}
$dropValidGenesStmt     = $dbh->prepare('DROP TABLE rnaSeqValidGenes');
$dropValidGenesStmt->execute()      or die $dropValidGenesStmt->errstr;
$dropMaxRankSpeciesStmt = $dbh->prepare('DROP TABLE rnaSeqMaxSpecies');
$dropMaxRankSpeciesStmt->execute()  or die $dropMaxRankSpeciesStmt->errstr;
if ( $auto == 0 ){
    $dbh->commit()  or die('Failed commit');
}

# Now, close and exit
$dbh->disconnect();

exit 0;

