#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

#############################
# Insert EST data for miRNA
#############################

use Getopt::Long;
use File::Slurp;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

$| = 1;    # stdout not out in memory buffer

# Define arguments & their default value
my ($species, $bgee_connector) = ('', '');
my ($smiRNAdbFile, $miRNA)     = ('', '');
my ($sex_info)                 = ('');
my ($Aport, $Sport)            = (0, 0);
my %opts = ('species=i'    => \$species,            # speciesCommonName from TSV for or Bgee db
            'bgee=s'       => \$bgee_connector,     # Bgee connector string
            'smiRNAdb=s'   => \$smiRNAdbFile,       # smiRNAdb libs mapping file
            'miRNA=s'      => \$miRNA,              # miRNA.dat file
            'sex_info=s'   => \$sex_info,           # generated_files/uberon/uberon_sex_info.tsv
            'Aport=i'      => \$Aport,              # ID MAPPING socket port
            'Sport=i'      => \$Sport,              # ID MAPPING stage socket port
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $species eq '' || $bgee_connector eq '' || $sex_info eq '' || $Aport == 0 || $Sport == 0 ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -species=9606 -bgee=\$(BGEECMD) -smiRNAdb=... -miRNA=miRNA.dat -sex_info=\$(UBERON_SEX_INFO_FILE_PATH) -Aport=\$(IDMAPPINGPORT) -Sport=\$(STGMAPPINGPORT)
\t-species   speciesId from Bgee db
\t-bgee      Bgee connector string
\t-smiRNAdb  smiRNAdb libs mapping file
\t-miRNA     miRNA.dat file
\t-sex_info  file containing sex-related info about anatomical terms
\t-Aport     ID MAPPING socket port
\t-Sport     ID MAPPING stage socket port
\n";
    exit 1;
}

# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);


# Get gene_mapping to bgeeGeneId
my $gene_mapping = Utils::query_bgeeGene($dbh, $species);


## Defining species aliases
my $spDB = $dbh->prepare('SELECT speciesId, CONCAT(SUBSTR(genus, 1, 1), SUBSTR(species, 1, 2)) FROM species WHERE speciesId = ?');
$spDB->execute($species)  or die $spDB->errstr;
my ($speciesCommonName, $speciesAlias) = $spDB->fetchrow_array;
# Adapt $speciesCommonName to the file nomenclature
# (to this point of the script, $speciesCommonName actually stored the speciesId...)
$speciesCommonName = $speciesCommonName ==  7955  ? 'Fish'
                   : $speciesCommonName ==  7227  ? 'Fly'
                   : $speciesCommonName == 10090  ? 'Mouse'
                   : $speciesCommonName ==  9606  ? 'Human'
                   :                                '';
if ( $speciesCommonName eq '' ){
    print "\tThis species is not annotated for miRNA EST [$species-$speciesAlias]\n";
    #NOTE Kill the mappers before the next species insertion
    Utils::get_anatomy_mapping([], $Aport, 0);
    Utils::get_anatomy_mapping([], $Sport, 0);
    exit 0;
}


# Check if data source exists
my $dsDB = $dbh->prepare('SELECT dataSourceId FROM dataSource WHERE dataSourceName = ? AND category = ?');
$dsDB->execute('smiRNAdb', 'EST data source')  or die $dsDB->errstr;
my $dataSourceId = $dsDB->fetchrow_array;
die "Data source smiRNAdb not found\n"  if ( !defined $dataSourceId );


# Ontology mapping
my @Anat;
my @Stg;
for my $line ( read_file("$smiRNAdbFile", chomp => 1) ){
    if ( $line =~ /\t$species\t[MF]?$/){
        my @tmp = split(/\t/, $line);
        push @Anat, $tmp[1];
        push @Stg,  $tmp[2];
    }
}
my $doneAnat = Utils::get_anatomy_mapping(\@Anat, $Aport, 0)  if ( exists $Anat[0] );
my $doneStg  = Utils::get_anatomy_mapping(\@Stg,  $Sport, 0)  if ( exists $Stg[0] );  # Same as anatomy mapping, not in between stages!


# Get already known conditions
my $conditions         = Utils::query_conditions($dbh);

# Get simpler (upper level) stage equivalences
my $stage_equivalences = Utils::get_stage_equivalences($dbh);

# Get sex-related information needed for sub 'insert_get_condition'
my $anatSexInfo    = Utils::get_anat_sex_info($sex_info);
my $speciesSexInfo = Utils::get_species_sex_info($dbh);

## Insert normal and annotated libraries into Bgee
my $estLibraryDB = $dbh->prepare('INSERT INTO estLibrary (estLibraryId, estLibraryName, conditionId, dataSourceId, estLibraryDescription) VALUES (?, ?, ?, ?, "")');
my $count = 0;
my $estLibraryId;
my %inserted_libs;
my $annotated_species = 0;
for my $line ( read_file("$smiRNAdbFile", chomp => 1) ){
    if ( $line =~ /\t$species\t[MF]?$/){
        my @tmp = split(/\t/, $line);
        if ( !exists $doneAnat->{$tmp[1]} || $doneAnat->{$tmp[1]} eq '' ){
            warn "[$tmp[1]] unmapped organ id\n";
            next;
        }
        if ( !exists $doneStg->{$tmp[2]}  || $doneStg->{$tmp[2]} eq '' ){
            warn "[$tmp[2]] unmapped stage id\n";
            next;
        }

        # Insert and get conditionId instead of anatEntityId/stageId to insert in estLibrary table
        # NOTE: sex and strain are currently NOT annotated for ESTs
        $count++;
        $estLibraryId = (lc $speciesAlias).'-mi'.$count;
        my $condKeyMap;
        ($condKeyMap, $conditions) = Utils::insert_get_condition($dbh, $conditions, $stage_equivalences,
                                                                 $doneAnat->{ $tmp[1] }, $doneStg->{ $tmp[2] },
                                                                 $species,
                                                                 $Utils::NOT_ANNOTATED_SEX, $Utils::NOT_ANNOTATED_STRAIN, # vars from Util.pm
                                                                 $anatSexInfo, $speciesSexInfo,
                                                                 $estLibraryId, '');

        $estLibraryDB->execute($estLibraryId, $tmp[0], $condKeyMap->{'conditionId'}, $dataSourceId)  or warn $estLibraryDB->errstr;
        $inserted_libs{$tmp[0]} = $estLibraryId; # used when parsing the other files
        $annotated_species = 1;
    }
}
if ( $annotated_species == 0 ){
    print "\tThis species is not annotated for miRNA EST [$species-$speciesAlias]\n";
    #NOTE Kill the mappers before the next species insertion
    Utils::get_anatomy_mapping([], $Aport, 0);
    Utils::get_anatomy_mapping([], $Sport, 0);
    exit 0;
}


## Insertion of lib keywords / description
my $upEstLib     = $dbh->prepare('UPDATE estLibrary SET estLibraryDescription = ? WHERE estLibraryId = ?');
my $selKeywordId = $dbh->prepare('SELECT keywordId FROM keyword WHERE keyword = ?');
my $insEstLib2K  = $dbh->prepare('INSERT INTO estLibraryToKeyword (estLibraryId, keywordId) VALUES (?, ?)');
my $insKeyword   = $dbh->prepare('INSERT INTO keyword (keyword) VALUES (?)');
my ($libName, $libId, $libDescription);
my @keywords;
my $keywordIdInBgee;
my $countLine = 0;
for my $line ( read_file('../../source_files/ESTs/smirnadb/S.csv', chomp => 1) ){
    $countLine++;
    my @tmp = split(/\t/, $line);
    if ( defined $tmp[1] && exists $inserted_libs{$tmp[1]} ){
        $libName        = $tmp[1];
        $libId          = $inserted_libs{$libName};
        $libDescription = $tmp[2];
        $upEstLib->execute($libDescription, $libId)  or die $upEstLib->errstr;

        # Keywords
        @keywords = split(/, /, $tmp[4]);
        push(@keywords, $tmp[3])  if ( $tmp[3] ne '' );

        # Check if the keywords are already defined in Bgee
        my %keywords;
        for my $keyword ( @keywords ){
            $keyword =~ s{^\s+}{};
            $keyword =~ s{\s+$}{};

            unless ( exists $keywords{lc($keyword)} ){
                $selKeywordId->execute( lc $keyword )  or die $selKeywordId->errstr;
                $keywordIdInBgee = $selKeywordId->fetchrow_array;
                if ( defined $keywordIdInBgee ){
                    $insEstLib2K->execute($libId, $keywordIdInBgee)  or warn $insEstLib2K->errstr;
                }
                else {
                    $insKeyword->execute( lc $keyword )  or warn $insKeyword->errstr;
                    # Retrieve the autoincremented keywordId
                    $keywordIdInBgee = $dbh->{mysql_insertid};
                    $insEstLib2K->execute($libId, $keywordIdInBgee)  or warn $insEstLib2K->errstr;
                }
                $keywords{ lc $keyword }++;
            }
        }
    }
    elsif ( defined $tmp[1] ){
        #warn "Warning, mirzId not found in mapping_libs_smiRNAdb.csv: $tmp[1] - line $countLine of S.csv: $line\n";
    }
}


## Insert ESTs into Bgee
# Retrieve miRNA geneBioTypeId
my $selGeneBioTypeId = $dbh->prepare('SELECT geneBioTypeId FROM geneBioType WHERE geneBioTypeName IN (?, ?)');
$selGeneBioTypeId->execute('miRNA', 'pre_miRNA')  or die $selGeneBioTypeId->errstr;
my $gene_biotype_id_ref = $selGeneBioTypeId->fetchall_arrayref;

# Retrieve miRBase dataSourceId
my $selDataSourceId = $dbh->prepare('SELECT dataSourceId FROM dataSource WHERE dataSourceName = ? AND category = ?');
$selDataSourceId->execute('miRBase', 'Genomics database')  or die $selDataSourceId->errstr;
my $data_source_id = $selDataSourceId->fetchrow_array;
die "Data source miRBase not found\n"  if ( !defined $data_source_id );

# Retrieve the mapping of mirbase names to Ensembl IDs
my %EnsemblToMiRBaseAcc;            # miRBase ID => geneId
my $selGeneXref = $dbh->prepare('SELECT t1.geneId, t2.XRefId FROM gene AS t1 INNER JOIN geneXRef AS t2 ON t1.bgeeGeneId = t2.bgeeGeneId WHERE geneBioTypeId IN (?, ?) AND dataSourceId = ? AND speciesId = ?');
$selGeneXref->execute($gene_biotype_id_ref->[0]->[0], $gene_biotype_id_ref->[1]->[0], $data_source_id, $species)  or die $selGeneXref->errstr;
while ( my @data = $selGeneXref->fetchrow_array ){
    $EnsemblToMiRBaseAcc{$data[1]}->{$data[0]}++; # miRBase ID => geneId
}

# Read miRNA.data for mapping miRBase Acc -> name
my ($miRNAName, $miRNAId) = ('', '');
my %miRBaseAccToName;
for my $line ( read_file("$miRNA", chomp => 1) ){
    if ( $line =~ m/^ID\s+(\w{3}-\S+)/ ){
        $miRNAName = $1;
    }
    elsif ( $line =~ m/^AC\s+(\S+);/ && $miRNAName ne '' ){
        $miRNAId = $1;
        $miRBaseAccToName{$miRNAId} = lc($miRNAName);
    }
    elsif ( $line =~ m/^\/\// ){
        $miRNAName = '';
        $miRNAId = '';
    }
}

# Make mapping mirbase name -> Ensembl ID
my %EnsemblToMiRBaseName;
my %EnsemblToMiRBaseName2;
for my $Acc ( keys %EnsemblToMiRBaseAcc ){
    for my $Id ( keys %{$EnsemblToMiRBaseAcc{$Acc}} ){
        if ( exists $miRBaseAccToName{$Acc} ){
            $EnsemblToMiRBaseName{$miRBaseAccToName{$Acc}}->{$Id}++;
            $EnsemblToMiRBaseName2{$Id} = $miRBaseAccToName{$Acc};
        }
    }
}


## Parse the smiRNAdb file
my $file = 'Report_'.$speciesCommonName.'.csv';
my %libIdToColumn;      # libId => column number in the file
my %geneIdToExpression; # geneName => expressionLib1, expressionLib2 .. expressionLibn
for my $line ( read_file("../../source_files/ESTs/smirnadb/$file", chomp => 1) ){
    # The first line contains the names of the libraries. We record the columns to read
    if ( $line =~ /^Sequence/gi ){
        my @tmp = split(/\t/, $line);
        # Check if the library is in Bgee
        for my $column ( 0..$#tmp ){
            my $libName = $tmp[$column];
            if ( defined $inserted_libs{$libName} ){
                # Map the id to a column in the file
                my $libId = $inserted_libs{$libName};
                $libIdToColumn{$libId} = $column;
            }
        }
    }
    # The other lines contain expression data
    elsif ( $line =~ /^[AUGC]+/gi ){
        # the sequence matches a miRNA sequence
        if ( lc($line) =~ /$speciesAlias/gi ){
            my @tmp = split(/\t/, $line);
            # the sequence should not match more than one miRNA (some lines contain more than one gene; see Report_Human.csv line 16)
            unless ( $tmp[1] =~ /\// ){
                my $fileGeneName = lc($tmp[1]);
                # if the gene is mapped to Ensembl
                if ( exists $EnsemblToMiRBaseName{$fileGeneName} ){
                    for my $EnsemblId ( keys %{$EnsemblToMiRBaseName{$fileGeneName}} ){
                        # If the gene name appears more than once in the file, add its expression data
                        if ( exists $geneIdToExpression{$EnsemblId} ){
                            $geneIdToExpression{$EnsemblId} = vectorSum($geneIdToExpression{$EnsemblId}, [@tmp[2..$#tmp]]);
                        }
                        else {
                            $geneIdToExpression{$EnsemblId} = [@tmp[2..$#tmp]];
                        }
                    }
                }
            }
        }
    }
}


## Interpret the data & insert them into bgee
my $insEST = $dbh->prepare('INSERT INTO expressedSequenceTag (estId, estLibraryId, bgeeGeneId) VALUES (?, ?, ?)');
for my $EnsemblId ( keys %geneIdToExpression ){
    for my $libId ( sort keys %libIdToColumn ){
        my $pos         = $libIdToColumn{$libId} - 2;  # The sequence and annotation fields have been removed
        my $countClones =  $geneIdToExpression{$EnsemblId}[$pos];
        unless ( $countClones eq 0 ){
            for my $c ( 1..$countClones ){
                $libId =~ /-(mi\d+)/;   # retrieve the number of the library to build the estId
                my $estId = $EnsemblToMiRBaseName2{$EnsemblId}.'_'.$EnsemblId.'_'.$1.'.'.$c;
                $insEST->execute($estId, $libId, $gene_mapping->{ $EnsemblId })  or warn $insEST->errstr;
            }
        }
    }
}



print "\tThis species IS     annotated for miRNA EST [$species-$speciesAlias]\n";
#$dbh->disconnect;
exit 0;


# add two vectors
sub vectorSum {
# References to the arrays
    my ($v1_ref, $v2_ref) = @_;

# Obtain the arrays
    my @v1 = @{$v1_ref};
    my @v2 = @{$v2_ref};

    my @ret;
    if ( scalar(@v1) != scalar(@v2) ){
        die "Error: trying to add two arrays of different length\n";
    }
    else {
        for ( my $i = 0; $i < scalar(@v1); $i++ ){
            push(@ret, ($v1[$i] + $v2[$i]));
        }
        return([@ret]);
    }
}

