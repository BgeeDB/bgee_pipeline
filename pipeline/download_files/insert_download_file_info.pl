#!/usr/bin/env perl

# Julien Wollbrett, May 2021

# This script is reponsible for inserting download files information

#############################################################

use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
use Data::Dumper;

$| = 1;

# Define arguments & their default value
my ($bgee_connector)                = ('');
my ($processed_dir)                 = ('');
my ($calls_dir)                     = ('');
my ($h5ad_dir)                      = ('');

my %opts = ('bgee=s'                 => \$bgee_connector,     # Bgee connector string
            'processed_dir=s'   => \$processed_dir,
            'calls_dir=s'       => \$calls_dir,
            'h5ad_dir=s'        => \$h5ad_dir
           );

my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $processed_dir eq '' ||
     $calls_dir eq '' || $h5ad_dir eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEE_CMD) -processed_dir=\$(PROCESSED_DIR) -calls_dir=\$(CALLS_DIR)
\t-bgee                 Bgee    connector string
\t-processed_dir        The path where processed expression download files have been generated
\t-calls_dir            The path where calls expression download files have been generated
\t-h5ad_dir             The path where h5ad download files have been generated
\n";
    exit 1;
}

my $calls_relative_path     = "calls/expr_calls/";
my $processed_relative_path = "processed_expr_values/";
my $h5ad_relative_path      = "h5ad/";

# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);

# retrieve species information (notably to generate file names)
# $species{speciesId}{'name'} = genus . species
my %species;

my $speciesSql = 'SELECT t1.speciesId, t1.genus, t1.species, t2.speciesDataGroupId FROM species 
                  as t1 inner join speciesToDataGroup as t2 ON t1.speciesId = t2.speciesId;';

my $speciesStmt = $dbh->prepare($speciesSql);
$speciesStmt->execute()  or die $speciesStmt->errstr;
while ( my @data = $speciesStmt->fetchrow_array ){
    my $speciesName = $data[1].' '.$data[2];
    $speciesName =~ s/ /_/g;
    $species{$speciesName}{'speciesId'} = $data[0];
    $species{$speciesName}{'dataGroupId'} = $data[3];
}

my %calls_files = detect_calls_files($calls_dir, $calls_relative_path, %species);


my %processed_files = detect_expression_files($processed_dir, $processed_relative_path, %species);

my %h5ad_files = detect_h5ad_files($h5ad_dir, $h5ad_relative_path, %species);

my %all_files = (%calls_files, %processed_files, %h5ad_files);
#print Dumper(\%all_files);

insert_into_database($dbh, \%all_files);


## Detect processed expression files

# all call files are generated in the same directory. The name of these files always
# finish with the extension .tsv.gz
# this function return the name of all calls files without the path to their directory and
# add a prefix before the name. This prefix correspond to the relative path from the download_files dir
# in the bgee_vXX dir of the FTP to the actual directory where the files are stored 
# (e.g /calls/expr_calls/)
sub detect_calls_files {
    my $directory = shift;
    my $prefix = shift;
    my $species = shift;
    my %all_calls_files_info;
    my $advanced_file_name = "_expr_advanced(.*).tsv.gz";
    my $simple_file_name = "_expr_simple(.*).tsv.gz";
    opendir(DIR, $directory) or die $!;
    while (my $calls_file = readdir(DIR)) {
        next unless (-f "$directory/$calls_file");
        next unless ($calls_file =~ m/\.tsv.gz$/);
        my $species_name = "";
        my $file_type = "";
        my $conditions = "";
        if ($calls_file =~ m/$advanced_file_name$/) {
            $calls_file =~ /(.*)_expr_advanced(.*).tsv.gz$/;
            $species_name = $1;
            $file_type = "expr_complete";
            $conditions = "anatomicalEntity";
            #TODO: for the moment we only generate files with anat entity or all conditions. Detection of conditions
            # has to be updated if we add other combinations of condition parameters.
            if ( defined($2) && $2 ne "") {
                $conditions = "anatomicalEntity,developmentalStage,sex,strain";
            }
        } elsif ($calls_file =~ m/$simple_file_name$/) {
            $calls_file =~ /(.*)_expr_simple(.*).tsv.gz$/;
            $species_name = $1;
            $file_type = "expr_simple";
            $conditions = "anatomicalEntity";
            #TODO: for the moment we only generate files with anat entity or all conditions. Detection of conditions
            # has to be updated if we add other combinations of condition parameters.
            if ( defined($2) && $2 ne "" ) {
                $conditions = "anatomicalEntity,developmentalStage,sex,strain";
            }
        }
        #print "$species_name\n";
        if (defined($species_name) && $species_name ne "") {
            $all_calls_files_info{$calls_file}{'conditions'} = $conditions;
            $all_calls_files_info{$calls_file}{'rel_path'} = "$prefix$calls_file";
            $all_calls_files_info{$calls_file}{'size'} = -s "$directory$calls_file";
            $all_calls_files_info{$calls_file}{'file_type'} = $file_type;
            $all_calls_files_info{$calls_file}{'data_group'} = $species{$species_name}{'dataGroupId'};;
        }
    }
    return %all_calls_files_info;
}

# expression files are generated in the different directories. The name of these files always
# finish with the extension .tar.gz
# this function should be improved by using perl modules to detect files in order to remove the ugly hardcoded 
# loops in directories.
# this function return the name of all expression files without the path to their directory and
# add a prefix before the name. This prefix correspond to the relative path from the download_files dir
# in the bgee_vXX dir of the FTP to the actual directory where the files are stored 
# (e.g processed_expr_values/rna_seq/)
sub detect_expression_files {
    my $directory = shift;
    my $prefix = shift;
    my $species = shift;
    my %all_calls_files_info;
    my $rnaseq_metadata_suffix = "_experiments_libraries.tar.gz";
    my $rnaseq_data_suffix = "read_counts_TPM.tar.gz";
    my $droplet_based_data_suffix = "read_counts_CPM.tar.gz";
    my $affy_metadata_suffix = "Affymetrix_experiments_chips.tar.gz";
    my $affy_data_suffix = "Affymetrix_probesets.tar.gz";
    # open directory containing datatypes dir
    opendir(DIR1, $directory) or die $!;
    while (my $datatypes_dir = readdir(DIR1)) {
        my $full_dt_dir = "$directory/$datatypes_dir";
        next unless (-d $full_dt_dir);
        if ($datatypes_dir =~ m/^affymetrix$/) {
            # open directory containing species dirs
            opendir(DIR2, $full_dt_dir) or die $!;
            while (my $species_dir = readdir(DIR2)) {
                my $full_species_dir = "$full_dt_dir/$species_dir";
                next unless (-d $full_species_dir);
                next unless (exists($species{$species_dir}));
                # open dir containing data
                opendir(DIR3, $full_species_dir) or die $!;
                while (my $data_file = readdir(DIR3)) {
                    #print "$data_file\n";
                    my $full_data_file = "$full_species_dir/$data_file";
                    next unless (-f $full_data_file);
                    if ($data_file =~ m/$affy_metadata_suffix$/) {
                        $all_calls_files_info{$data_file}{'conditions'} = undef;
                        $all_calls_files_info{$data_file}{'rel_path'} = "$prefix$datatypes_dir/$species_dir/$data_file";
                        $all_calls_files_info{$data_file}{'size'} = -s $full_data_file;
                        $all_calls_files_info{$data_file}{'data_group'} = $species{$species_dir}{'dataGroupId'};
                        $all_calls_files_info{$data_file}{'file_type'} = 'affy_annot';
                    } elsif ($data_file =~ m/$affy_data_suffix$/) {
                        $all_calls_files_info{$data_file}{'conditions'} = undef;
                        $all_calls_files_info{$data_file}{'rel_path'} = "$prefix$datatypes_dir/$species_dir/$data_file";
                        $all_calls_files_info{$data_file}{'size'} = -s $full_data_file;
                        $all_calls_files_info{$data_file}{'data_group'} = $species{$species_dir}{'dataGroupId'};
                        $all_calls_files_info{$data_file}{'file_type'} = 'affy_data';
                    }
                }
            }
        } elsif ($datatypes_dir =~ m/^rna_seq_processed$/) {
            # open directory containing species dirs
            opendir(DIR2, $full_dt_dir) or die $!;
            while (my $species_dir = readdir(DIR2)) {
                my $full_species_dir = "$full_dt_dir/$species_dir";
                next unless (-d $full_species_dir);
                next unless (exists($species{$species_dir}));
                # open dir containing data
                opendir(DIR3, $full_species_dir) or die $!;
                while (my $data_file = readdir(DIR3)) {
                    my $full_data_file = "$full_species_dir/$data_file";
                    next unless (-f $full_data_file);
                    if ($data_file =~ m/$rnaseq_metadata_suffix$/) {
                        $all_calls_files_info{$data_file}{'conditions'} = undef;
                        $all_calls_files_info{$data_file}{'rel_path'} = $prefix."rna_seq/$species_dir/$data_file";
                        $all_calls_files_info{$data_file}{'size'} = -s $full_data_file;
                        $all_calls_files_info{$data_file}{'data_group'} = $species{$species_dir}{'dataGroupId'};
                        $all_calls_files_info{$data_file}{'file_type'} = 'rnaseq_annot';

                    } elsif ($data_file =~ m/$rnaseq_data_suffix$/) {
                        $all_calls_files_info{$data_file}{'conditions'} = undef;
                        $all_calls_files_info{$data_file}{'rel_path'} = $prefix."rna_seq/$species_dir/$data_file";
                        $all_calls_files_info{$data_file}{'size'} = -s $full_data_file;
                        $all_calls_files_info{$data_file}{'data_group'} = $species{$species_dir}{'dataGroupId'};
                        $all_calls_files_info{$data_file}{'file_type'} = 'rnaseq_data';
                    }
                }
            }
        } elsif ($datatypes_dir =~ m/^sc_rna_seq_fl_processed$/) {
            # open directory containing species dirs
            opendir(DIR2, $full_dt_dir) or die $!;
            while (my $species_dir = readdir(DIR2)) {
                my $full_species_dir = "$full_dt_dir/$species_dir";
                next unless (-d $full_species_dir);
                next unless (exists($species{$species_dir}));
                # open dir containing data
                opendir(DIR3, $full_species_dir) or die $!;
                while (my $data_file = readdir(DIR3)) {
                    my $full_data_file = "$full_species_dir/$data_file";
                    next unless (-f $full_data_file);
                    if ($data_file =~ m/$rnaseq_metadata_suffix$/) {
                        $all_calls_files_info{$data_file}{'conditions'} = undef;
                        $all_calls_files_info{$data_file}{'rel_path'} = $prefix."full_length/$species_dir/$data_file";
                        $all_calls_files_info{$data_file}{'size'} = -s $full_data_file;
                        $all_calls_files_info{$data_file}{'data_group'} = $species{$species_dir}{'dataGroupId'};
                        $all_calls_files_info{$data_file}{'file_type'} = 'full_length_annot';
                    } elsif ($data_file =~ m/$rnaseq_data_suffix$/) {
                        $all_calls_files_info{$data_file}{'conditions'} = undef;
                        $all_calls_files_info{$data_file}{'rel_path'} = $prefix."full_length/$species_dir/$data_file";
                        $all_calls_files_info{$data_file}{'size'} = -s $full_data_file;
                        $all_calls_files_info{$data_file}{'data_group'} = $species{$species_dir}{'dataGroupId'};
                        $all_calls_files_info{$data_file}{'file_type'} = 'full_length_data';
                    }
                }
            }
        } elsif ($datatypes_dir =~ m/^sc_rna_seq_db_processed$/) {
            # open directory containing species dirs
            opendir(DIR2, $full_dt_dir) or die $!;
            while (my $species_dir = readdir(DIR2)) {
                my $full_species_dir = "$full_dt_dir/$species_dir";
                next unless (-d $full_species_dir);
                next unless (exists($species{$species_dir}));
                # open dir containing data
                opendir(DIR3, $full_species_dir) or die $!;
                while (my $data_file = readdir(DIR3)) {
                    my $full_data_file = "$full_species_dir/$data_file";
                    next unless (-f $full_data_file);
                    if ($data_file =~ m/$rnaseq_metadata_suffix$/) {
                        $all_calls_files_info{$data_file}{'conditions'} = undef;
                        $all_calls_files_info{$data_file}{'rel_path'} = $prefix."droplet_based/$species_dir/$data_file";
                        $all_calls_files_info{$data_file}{'size'} = -s $full_data_file;
                        $all_calls_files_info{$data_file}{'data_group'} = $species{$species_dir}{'dataGroupId'};
                        $all_calls_files_info{$data_file}{'file_type'} = 'droplet_based_annot';
                    } elsif ($data_file =~ m/$droplet_based_data_suffix$/) {
                        $all_calls_files_info{$data_file}{'conditions'} = undef;
                        $all_calls_files_info{$data_file}{'rel_path'} = $prefix."sdroplet_based/$species_dir/$data_file";
                        $all_calls_files_info{$data_file}{'size'} = -s $full_data_file;
                        $all_calls_files_info{$data_file}{'data_group'} = $species{$species_dir}{'dataGroupId'};
                        $all_calls_files_info{$data_file}{'file_type'} = 'droplet_based_data';
                    }
                }
            }
        }
    }
    return %all_calls_files_info;
}

sub detect_h5ad_files {
    my $directory = shift;
    my $prefix = shift;
    my $species = shift;
    my %h5ad_files_info;
    my $h5ad_suffix = "_h5ad.tar.gz";
    # open directory containing datatypes dir
    opendir(DIR1, $directory) or die $!;
    while (my $datatypes_dir = readdir(DIR1)) {
        my $full_dt_dir = "$directory/$datatypes_dir";
        print "dt dir : $datatypes_dir\n";
        next unless (-d $full_dt_dir);
        if ($datatypes_dir =~ m/^droplet_based$/) {
            # open directory containing species dirs
            opendir(DIR2, $full_dt_dir) or die $!;
            while (my $species_dir = readdir(DIR2)) {
                my $full_species_dir = "$full_dt_dir/$species_dir";
                print "spe dir : $full_species_dir\n";
                next unless (-d $full_species_dir);
                next unless (exists($species{$species_dir}));
                # open dir containing data
                opendir(DIR3, $full_species_dir) or die $!;
                while (my $data_file = readdir(DIR3)) {
                    my $full_data_file = "$full_species_dir/$data_file";
                    next unless (-f $full_data_file);
                    if ($data_file =~ m/$h5ad_suffix$/) {
                        print "data file : $full_data_file\n";
                        $h5ad_files_info{$data_file}{'conditions'} = undef;
                        $h5ad_files_info{$data_file}{'rel_path'} = $prefix."droplet_based/$species_dir/$data_file";
                        $h5ad_files_info{$data_file}{'size'} = -s $full_data_file;
                        $h5ad_files_info{$data_file}{'data_group'} = $species{$species_dir}{'dataGroupId'};
                        $h5ad_files_info{$data_file}{'file_type'} = 'droplet_based_h5ad';

                    } 
                }
            }
        } elsif ($datatypes_dir =~ m/^full_length$/) {
            # open directory containing species dirs
            opendir(DIR2, $full_dt_dir) or die $!;
            while (my $species_dir = readdir(DIR2)) {
                my $full_species_dir = "$full_dt_dir/$species_dir";
                next unless (-d $full_species_dir);
                next unless (exists($species{$species_dir}));
                # open dir containing data
                opendir(DIR3, $full_species_dir) or die $!;
                while (my $data_file = readdir(DIR3)) {
                    my $full_data_file = "$full_species_dir/$data_file";
                    next unless (-f $full_data_file);
                    if ($data_file =~ m/$h5ad_suffix$/) {
                        $h5ad_files_info{$data_file}{'conditions'} = undef;
                        $h5ad_files_info{$data_file}{'rel_path'} = $prefix."full_length/$species_dir/$data_file";
                        $h5ad_files_info{$data_file}{'size'} = -s $full_data_file;
                        $h5ad_files_info{$data_file}{'data_group'} = $species{$species_dir}{'dataGroupId'};
                        $h5ad_files_info{$data_file}{'file_type'} = 'full_length_h5ad';

                    } 
                }
            }
        }
    } 
    return %h5ad_files_info;
}

sub insert_into_database {
    my $bgee = shift;
    my $data_to_insert = shift;
    print "Insert download files info...\n";
    #TODO: add a check that downloadFile table is empty. Otherwise it will be possible to insert 
    # information about the same file twice. We could also retrieve file names and check that a file is not already present
    my $insert = $dbh->prepare('INSERT INTO downloadFile (downloadFileRelativePath, downloadFileName, downloadFileDescription, downloadFileCategory,
        speciesDataGroupId, downloadFileSize, downloadFileConditionParameters) VALUES (?, ?, ?, ?, ?, ?, ?)');
    for my $file_name ( keys %{$data_to_insert} ){
        my $cond = $data_to_insert->{$file_name}->{'conditions'};
        my $file_path = $data_to_insert->{$file_name}->{'rel_path'};
        my $size = $data_to_insert->{$file_name}->{'size'};
        my $data_group = $data_to_insert->{$file_name}->{'data_group'};
        my $file_type = $data_to_insert->{$file_name}->{'file_type'};
        print "INSERT INTO downloadFile (downloadFileRelativePath, downloadFileName, downloadFileDescription, downloadFileCategory,".
          " speciesDataGroupId, downloadFileSize, downloadFileConditionParameters) VALUES ($file_path, $file_name, \"\", $file_type, $data_group, $size, $cond)\n";
        $insert->execute($file_path, $file_name, "", $file_type, $data_group, $size, $cond)  or die $insert->errstr;
    }
    $insert->finish();
}


