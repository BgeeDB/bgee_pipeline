#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use File::Slurp;

use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;

$| = 1; # stdout not put in memory buffer


# Define arguments & their default value
my ($species, $bgee_connector, $ensembl_connector) = ('', '', '');
my ($ESTspeciesFile, $mappingFile, $libReport)     = ('', '', '');
my ($sex_info)                                     = ('');
my ($Aport, $Sport)                                = (0, 0);
my %opts = ('species=i'    => \$species,            # speciesId from TSV for or Bgee db
            'bgee=s'       => \$bgee_connector,     # Bgee connector string
            'ensembl=s'    => \$ensembl_connector,  # Ensembl connector string
            'ESTspecies=s' => \$ESTspeciesFile,     # file with Species EST annotations
            'mapping=s'    => \$mappingFile,        # mapping file between unigene & sequences, if no mapping in Ensembl
            'libReport=s'  => \$libReport,          # library.report file
            'sex_info=s'   => \$sex_info,           # generated_files/uberon/uberon_sex_info.tsv
            'Aport=i'      => \$Aport,              # ID MAPPING socket port
            'Sport=i'      => \$Sport,              # ID MAPPING stage socket port
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $species eq '' || $bgee_connector eq '' || $ensembl_connector eq '' || $sex_info eq '' || $Aport == 0 || $Sport == 0 ){
    print "\n\tInvalid or missing argument:
\te.g. $0 -species=9606 -bgee=\$(BGEECMD) -ensembl=\$(ENSCMD) -ESTspecies=... -libReport=library.report [-mapping=...] -sex_info=\$(UBERON_SEX_INFO_FILE_PATH) -Aport=\$(IDMAPPINGPORT) -Sport=\$(STGMAPPINGPORT)
\t-species    speciesId from Bgee db
\t-bgee       Bgee connector string
\t-ensembl    Ensembl connector string
\t-ESTspecies file with Species EST annotations
\t-mapping    mapping file between unigene & sequences, if no mapping in Ensembl
\t-libReport  library.report file
\t-sex_info   file containing sex-related info about anatomical terms
\t-Aport      ID MAPPING socket port
\t-Sport      ID MAPPING stage socket port
\n";
    exit 1;
}

# Ensembl connection via Ensembl API/Registry
my $reg = Utils::connect_ensembl_registry($ensembl_connector);
# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);


# Get full species ids
my $speciesDB = $dbh->prepare('SELECT genus, species, speciesCommonName FROM species WHERE speciesId=?');
$speciesDB->execute($species)  or die $speciesDB->errstr;
my $species_ref = $speciesDB->fetchall_arrayref;


# Get gene_mapping to bgeeGeneId
my $gene_mapping = Utils::query_bgeeGene($dbh, $species);


# Check UniGene dataSource
my $unigeneDB = $dbh->prepare('SELECT dataSourceId FROM dataSource WHERE dataSourceName=? AND category=?');
$unigeneDB->execute('UniGene', 'EST data source')  or die $unigeneDB->errstr;
my $unigene = $unigeneDB->fetchall_arrayref;
die "UniGene data source not found\n"  if ( !exists $unigene->[0]->[0] );


## Insert normal and annotated libraries into Bgee
#TODO merge annotation_libs_SPECIES.txt files AND pass it as script argument
my $annot_file = '../../source_files/ESTs/unigene/annotation_libs_'.$species_ref->[0]->[0].'_'.$species_ref->[0]->[1].'.txt';
# quit if no or empty annotation file
if ( !-e "$annot_file" || -z "$annot_file" ){
    print "\tNo UniGene Bgee annotation for $species_ref->[0]->[0] $species_ref->[0]->[1] ($species)\n";
    #NOTE Kill the mappers before the next species insertion
    Utils::get_anatomy_mapping([], $Aport, 0);
    Utils::get_anatomy_mapping([], $Sport, 0);
    exit 0;
}

# Read annotation_libs_SPECIES.txt file and put columns in hash according to header
#                                 file_name, separator, parser, sheet number
# XXX: is the quote mode '"' when launching read_spreadsheet?
my %tsv = %{ Utils::read_spreadsheet("$annot_file", "\t", 'csv', '', 1) };


# Ontology mapping
my @Anat = @{$tsv{'anatEntityId'}};
my @Stg  = @{$tsv{'stageId'}};
my $doneAnat = Utils::get_anatomy_mapping(\@Anat, $Aport, 0);
my $doneStg  = Utils::get_anatomy_mapping(\@Stg,  $Sport, 0); # Same as anatomy mapping, not in between stages!


# Get already known conditions
my $conditions         = Utils::query_conditions($dbh);

# Get simpler (upper level) stage equivalences
my $stage_equivalences = Utils::get_stage_equivalences($dbh);

# Get sex-related information needed for sub 'insert_get_condition'
my $anatSexInfo    = Utils::get_anat_sex_info($sex_info);
my $speciesSexInfo = Utils::get_species_sex_info($dbh);


# Insert estLibrary
my $insUniGeneDB = $dbh->prepare('INSERT INTO estLibrary (estLibraryId, estLibraryName, conditionId, dataSourceId, estLibraryDescription) VALUES (?, ?, ?, ?, "")');
my %inserted_libs;
for my $line ( 0..$#{$tsv{'estLibraryId'}} ){
#    warn "[$tsv{'estLibraryId'}[$line]]\t[$tsv{'estLibraryName'}[$line]]\t[$tsv{'anatEntityId'}[$line]]\t[$tsv{'stageId'}[$line]]\t[$unigene->[0]->[0]]\n";
    if ( !exists $doneAnat->{ $tsv{'anatEntityId'}[$line] } || $doneAnat->{ $tsv{'anatEntityId'}[$line] } eq '' ){
        warn "[$tsv{'anatEntityId'}[$line]] unmapped organ id\n";
        next;
    }
    if ( !exists $doneStg->{ $tsv{'stageId'}[$line] } || $doneStg->{ $tsv{'stageId'}[$line] } eq '' ){
        warn "[$tsv{'stageId'}[$line]] unmapped stage id\n";
        next;
    }

    # Insert and get conditionId instead of anatEntityId/stageId to insert in estLibrary table
    # NOTE: sex and strain are currently NOT annotated for ESTs
    my $condKeyMap;
    ($condKeyMap, $conditions) = Utils::insert_get_condition($dbh, $conditions, $stage_equivalences,
                                                             $doneAnat->{ $tsv{'anatEntityId'}[$line] }, $doneStg->{ $tsv{'stageId'}[$line] },
                                                             $species,
                                                             $Utils::NOT_ANNOTATED_SEX, $Utils::NOT_ANNOTATED_STRAIN, # vars from Util.pm
                                                             $anatSexInfo, $speciesSexInfo,
                                                             $tsv{'estLibraryId'}[$line], '');

    $insUniGeneDB->execute($tsv{'estLibraryId'}[$line], $tsv{'estLibraryName'}[$line],
                           $condKeyMap->{'conditionId'}, $unigene->[0]->[0])
        or warn $insUniGeneDB->errstr;
    $inserted_libs{ $tsv{'estLibraryId'}[$line] }++;
}


## Insert EST lib keywords / description
my $insKeywordDB  = $dbh->prepare('INSERT INTO keyword (keyword) VALUES (?)');
my $insEST2KeywDB = $dbh->prepare('INSERT INTO estLibraryToKeyword (estLibraryId, keywordId) VALUES (?, ?)');
my $upESTDB       = $dbh->prepare('UPDATE estLibrary SET estLibraryDescription = ? WHERE estLibraryId = ?');
my $keywordDB     = $dbh->prepare('SELECT keywordId FROM keyword WHERE keyword = ?');

# Read library.report
my $flag = 0;
my ($lib_id, $lib_description) = ('', '');
my %keywords;
for my $line ( read_file("$libReport", chomp => 1) ){
    # this library is in Bgee
    if ( $line =~ /^dbEST lib id:\s*(.+?)\s*$/ && exists $inserted_libs{$1} ){
        $lib_id = $1;
        $flag = 1;
    }

    # Get lib description
    if ( $line =~ /^Description:\s*(.+?)\s*$/ && $flag == 1 ){
        $lib_description = $1;
        $lib_description =~ s{  +}{ }g;
    }

    # Get keywords
    if ( $line =~ /^Keyword:\s*(.+?)\s*$/ && $flag == 1 ){
        $keywords{$1}++;
    }

    # Avoid insertion of library without any sequences
    if ( $line eq 'Sequences generated to date: 0' ){
        $flag = 0;
    }

    # End of description of the lib
    if ( $line =~ /^> .+$/ && $flag == 1 ){
        KEYWORD:
        for my $word ( sort keys %keywords ){
            $keywordDB->execute($word)  or die $keywordDB->errstr;
            my $k_id = $keywordDB->fetchrow_array;
            # Check if the keyword already exists in the table
            if ( defined $k_id ){
                $insEST2KeywDB->execute($lib_id, $k_id)  or warn $insEST2KeywDB->errstr;
            }
            else { # Or insert the new keyword
                $insKeywordDB->execute($word)  or warn $insKeywordDB->errstr;
                # Retrieve the auto_incremented keywordId:
                # MySQL has the ability to choose unique key values automatically. If this
                # happened, the new ID will be stored in this attribute. An alternative
                # way for accessing this attribute is via $dbh->{'mysql_insertid'}
                $k_id = $dbh->{'mysql_insertid'};
                $insEST2KeywDB->execute($lib_id, $k_id)  or warn  $insEST2KeywDB->errstr;
            }
        }
        if ( $lib_description ne '' ){
            $upESTDB->execute($lib_description, $lib_id)  or die $upESTDB->errstr;
        }
    }

    if ( $line =~ /^\> .+$/ ){
        $flag            = 0;
        %keywords        = ();
        $lib_description = '';
        $lib_id          = '';
    }
}


## Insert ESTs into Bgee
# Retrieve the file where the EST Unigene cluster infos are
my $est_file;
for my $line ( read_file("$ESTspeciesFile", chomp => 1 ) ){
    my @tmp = split(/\t/, $line);
    my $db_species = $species_ref->[0]->[2];
    # Fix species name
    $db_species = 'fruitfly'  if ( $db_species eq 'fruit fly' );
    $db_species = 'xenopus'   if ( $db_species eq 'western clawed frog' );
    if ( lc $tmp[0] eq lc $db_species ){
        $est_file = $tmp[1];
        last;
    }
}

## Retrieve mapping from Biomart
my %mapping;
# if a mapping file is given
if ( $mappingFile ne '' ){
    for my $line ( read_file("$mappingFile", chomp => 1 ) ){
        my @tmp = split(/\t/, $line);
        push(@{$mapping{$tmp[0]}}, $tmp[1]);
    }
}
else {
    # Get a slice adaptor for the  $species  core database
    my $gene_adaptor = $reg->get_adaptor( $species, 'Core', 'Gene' );

    # Fetch all clones from a slice adaptor (returns a list reference)
    my @genes = @{$gene_adaptor->fetch_all()};
    GENE:
    for my $gene ( sort {$a->stable_id cmp $b->stable_id} (@genes) ){ #Sort to always get the same order
        my @unigene = map  { $_->display_id }
                      grep { lc $_->dbname() eq 'unigene' } # only UniGene terms
                      @{$gene->get_all_xrefs()};

        if ( !exists $unigene[0] ){
#            warn "\t\t", 'UniGene mapping was not retrieved for [', $gene->stable_id, "]\n";
            next GENE;
        }
        for my $unig ( sort @unigene ){
            push(@{$mapping{ $unig }}, $gene->stable_id);
        }
    }
}


## Filling BGee
# Check UniGene EST file availability and gunzip it if required
if ( -e "$est_file" ){
    # Fine
}
elsif ( -e "$est_file.gz" ){
    system("gunzip --to-stdout $est_file.gz > $est_file")==0  or die "Cannot gunzip [$est_file.gz]\n";
}
else {
    die "Cannot find [$est_file]\n";
}
# Use bgeeGeneId instead of geneId through $gene_mapping
my $estDB = $dbh->prepare('INSERT INTO expressedSequenceTag (estId, estId2, estLibraryId, bgeeGeneId, UniGeneClusterId) VALUES (?, ?, ?, ?, ?)');
my $unigene_ID;
my (@EST_ACC, @EST_GI, @LID);
# Read UniGene EST file for insertion
for my $line ( read_file("$est_file", chomp => 1 ) ){
    if ( $line =~ /^ID\s+(\w+\.\d+)/ ){
        $unigene_ID = $1; # Retrieve Unigene cluster ID
    }

    # Retrieve ESTs IDs and infos
    if ( $line =~ /^SEQUENCE\s+ACC=(\w+\.?\w*);\s+NID=(\w+);.+?LID=(\w+);.+?SEQTYPE=EST.*/ && exists $inserted_libs{$3} ){
        push(@EST_ACC, $1);
        push(@EST_GI,  $2);
        push(@LID,     $3);
    }

    # End of a UniGene record -> insertion into Bgee
    if ( $line =~ /^\/\// ){
        if ( exists $EST_ACC[0] ){
            # If at least one EST to insert in this unigene cluster -> Retrieving EnsemblID corresponding to Unigene cluster ID
            # We consider only genes with a one to one relation btwn UniGene and Ensembl
            if ( $#{$mapping{$unigene_ID}} == 0 ){
                for (my $i=0; $i<=$#EST_ACC; $i++){
                    $estDB->execute($EST_ACC[$i], $EST_GI[$i], $LID[$i], $gene_mapping->{ ${$mapping{$unigene_ID}}[0] }, $unigene_ID)  or warn $estDB->errstr;
                }
            }
        }

        $unigene_ID = undef;
        @EST_ACC    = ();
        @EST_GI     = ();
        @LID        = ();
    }
}
# Delete gunzip file to keep only gzip ones
unlink "$est_file";


# Close db connections
$reg->clear();
#$dbh->disconnect;

exit 0;

