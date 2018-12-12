#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

# Frederic Bastian, created November 2012
# Julien Roux, updated March 2016, October 2016
# list of subs useful for analyzing RNA-Seq data
#############################################################

use List::MoreUtils qw(all any);
use LWP::Simple;
use File::Slurp;
use Digest::SHA;
$| = 1;


my $floatingPointRegex   = '^[-+]?\d*\.?\d+([eE][-+]?\d+)?$';
my $integerRegex         = '^[-+]?\d+$';
my $positiveIntegerRegex = '^\d+$';

sub bgeeTrim {
    my ($stringToTrim) = @_;
    if ( defined $stringToTrim ){
        $stringToTrim =~ s/^\s+|\s+$//g ;
    }
    return $stringToTrim;
}

# Retrieve the library ID (SRX...) from the sample ID (GSMxxx) by requesting the ncbi website
# and parsing the webpage
sub retrieveLibraryId {
    my ($sampleId) = @_;
    # to not overload the server
    sleep 10;
    my $fileName = 'temp_ncbi_download';
    # var to check that the downloaded file has the expected format:
    # line where the "Relations" title is reached
    my $relationsLine = 0;
    # line where the "SRA" column is reached
    my $sraLine       = 0;
    my $libId = undef;

    my $statusCode = getstore('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='.$sampleId, $fileName);
    if ( $statusCode != 200 ){
        die 'Warning: could not get https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=', $sampleId, ': error ', $statusCode, "\n";
    }

    my $lineCount = 0;
    for my $line ( read_file($fileName, chomp => 1) ){
        if ( $line =~ /<strong>Relations<\/strong>/ ){
            $relationsLine = $lineCount;
        }
        elsif ( $line =~ /<td>SRA<\/td>/ ){
            $sraLine = $lineCount;
        }
        # OK, now it should be the line from which to extract the SRX ID
        elsif ( $lineCount == $relationsLine + 2 && $lineCount == $sraLine + 1 && $line =~ /<a.+?>([SEDC]RX.+?)<\/a>/ ){
            $libId = bgeeTrim($1);
            last;
        }
        $lineCount++;
    }
    unlink($fileName)  or die("Could not remove [$fileName]\n");
    if ( !defined $libId || $libId !~ /^[SEDC]RX/ ){
        print("Could not retrieve the library ID for sampleId [$sampleId]\n");
    }
    return $libId;
}
# TODO is this function still used in pipeline? Isn't it redundant with get_SRA.pl?

# Generate SHA512 checksum for a given file
sub generateCheckSum {
    my ($file, $removeFirstLine) = @_;
    # $removeFirstLine is optional argument
    # to define whether the first line of the file should be removed before generating the checksum
    # it is the case for filtered processed mas5, where lines are ordered by probeset IDs, but header is not standardized.
    # better chance to detect duplicates by removing it
    my $sha = Digest::SHA->new(512);

    if ( defined $removeFirstLine && $removeFirstLine ){
        open(my $TEMPOUT, '>', 'firtLineRemovedTempFile')  or die "could not open temp file\n";
        open(my $IN, '<', $file)                           or die "could not read [$file]\n";
        my $line = <$IN>;#first line removed
        while ( defined ($line = <$IN>) ){
            print {$TEMPOUT} $line;
        }
        close $IN;
        close $TEMPOUT;

        $sha->addfile('firtLineRemovedTempFile', 'p');
        unlink('firtLineRemovedTempFile');
    }
    else {
        $sha->addfile($file, 'p');
    }

    return $sha->hexdigest();
}
# TODO still used in pipeline? Keep?

# Extract experiments information from the experiment annotation file
sub getAllAnnotatedExperiments2 {
    my ($TSV) = @_;
    my %tsv = %$TSV;
    # $experiments{experimentId}->{'name'}
    # $experiments{experimentId}->{'description'}
    # $experiments{experimentId}->{'source'}
    # $experiments{experimentId}->{'status'}
    # $experiments{experimentId}->{'commented'}
    my %experiments = ();

    my $totalLineCount = 0;
    for my $line ( 0..$#{$tsv{'experimentId'}} ){
        my $commented = 0;
        if ( ($tsv{'experimentId'}[$line] =~ /^#(.+)/) or ($tsv{'experimentId'}[$line] =~ /^\"#(.+)/) ){
            $tsv{'experimentId'}[$line] = $1;
            $commented = 1;
        }

        if ( length($tsv{'experimentId'}[$line]) < 255 && $tsv{'experimentId'}[$line] !~ /\s/ ){
            # Initialize the hash for this expId
            $experiments{ $tsv{'experimentId'}[$line] }->{'name'}        = '';
            $experiments{ $tsv{'experimentId'}[$line] }->{'description'} = '';
            $experiments{ $tsv{'experimentId'}[$line] }->{'source'}      = '';
            $experiments{ $tsv{'experimentId'}[$line] }->{'status'}      = '';
            $experiments{ $tsv{'experimentId'}[$line] }->{'commented'}   = '';
        }
        else {
            warn "Badly formatted RNAseqExperiment file at line [$line]: [invalid experimentId]\n";
            next;
        }

        $experiments{ $tsv{'experimentId'}[$line] }->{'name'}        = $tsv{'experimentName'}[$line];
        $experiments{ $tsv{'experimentId'}[$line] }->{'description'} = $tsv{'experimentDescription'}[$line];
        $experiments{ $tsv{'experimentId'}[$line] }->{'source'}      = $tsv{'experimentSource'}[$line];
        $experiments{ $tsv{'experimentId'}[$line] }->{'status'}      = $tsv{'experimentStatus'}[$line];
        $experiments{ $tsv{'experimentId'}[$line] }->{'commented'}   = $commented;
    }

    return %experiments;
}

## TODO remove? This is the old format for RNASeqExperiment.tsv
sub getAllAnnotatedExperiments {
    my ($rnaSeqExpAnnotationFile) = @_;
    # $experiments{experimentId}->{'name'}
    # $experiments{experimentId}->{'description'}
    # $experiments{experimentId}->{'source'}
    # $experiments{experimentId}->{'status'}
    my %experiments = ();

    my $lineCount      = 0;
    my $totalLineCount = 0;
    my $error          = '';
    my $expId          = undef;
    for my $line ( read_file($rnaSeqExpAnnotationFile, chomp => 1) ){
        # file format:
        # first line: experimentId\texpId
        # second line: experimentName\texpName
        # third line: experimentDescription\tdescr.
        # fourth line: experimentSource\tsource
        # fifth line: experimentStatus\status
        # 6th line: either optional COMMENT, or experiment separator //

        my @tmp;
        ($tmp[0], $tmp[1]) = map { bgeeTrim($_) } split(/\t/, $line);
        # first line describing an experiment
        if ( $lineCount == 0 && !defined $expId ){
            if ( $tmp[0] eq 'experimentId' && length($tmp[1]) < 255 && $tmp[1] !~ /\s/ ){
                $expId = $tmp[1];
                # Initialize the hash for this expId
                $experiments{$expId}->{'name'}        = '';
                $experiments{$expId}->{'description'} = '';
                $experiments{$expId}->{'source'}      = '';
                $experiments{$expId}->{'status'}      = '';
            }
            else {
                $error = 'missing experimentId';
            }
        }
        elsif ( $lineCount == 1 && defined $expId ){
            if ( $tmp[0] eq 'experimentName' ){
                $experiments{$expId}->{'name'} = $tmp[1];
            }
            else {
                $error = 'missing experimentName';
            }
        }
        elsif ( $lineCount == 2 && defined $expId ){
            if ( $tmp[0] eq 'experimentDescription' ){
                $experiments{$expId}->{'description'} = $tmp[1];
            }
            else {
                $error = 'missing experimentDescription';
            }
        }
        elsif ( $lineCount == 3 && defined $expId ){
            if ( $tmp[0] eq 'experimentSource' ){
                $experiments{$expId}->{'source'} = $tmp[1];
            }
            else {
                $error = 'missing experimentSource';
            }
        }
        elsif ( $lineCount == 4 && defined $expId ){
            if ( $tmp[0] eq 'experimentStatus' ){
                $experiments{$expId}->{'status'} = $tmp[1];
            }
            else {
                $error = 'missing experimentStatus';
            }
        } elsif ( $lineCount == 5 && defined $expId ){
            if ( $tmp[0] eq 'COMMENT' ){
                # optionnal comment line, everything's OK
                # $line-- to match "//" during next iteration
                $lineCount--;
            }
            elsif ( $line =~ /^\/\/\s*$/ ){
                #everything's OK
            }
            else {
                $error = 'missing COMMENT or end of experiment //';
            }
        }
        else {
            warn "Error, no experiment ID already defined when reaching line [$totalLineCount]\n";
        }
        $totalLineCount++;
        if ( $error ne '' ){
            warn "Badly formatted rnaSeaExperiment file at line [$totalLineCount]: [$error]\n";
            $error = '';
        }
        if ( $line =~ /^\/\/\s*$/ ){
            $expId     = undef;
            $lineCount = -1;
            $error     = '';
        }

        $lineCount++;
    }
    return %experiments;
}

## Problem with this: the # at beginning of lines marks header by also commented libraries
sub getAllRnaSeqAnnotations2 {
    my ($TSV) = @_;
    my %tsv = %$TSV;
    # $rnaSeqAnnotations{expId}->{libraryId (SRXxxx)}->{'platform'}   = platform
    # $rnaSeqAnnotations{expId}->{libraryId (SRXxxx)}->{'uberonId'}   = uberonId
    # $rnaSeqAnnotations{expId}->{libraryId (SRXxxx)}->{'stageId'}    = stageId
    # $rnaSeqAnnotations{expId}->{libraryId (SRXxxx)}->{'sex'}        = sex
    # $rnaSeqAnnotations{expId}->{libraryId (SRXxxx)}->{'strain'}     = strain
    # $rnaSeqAnnotations{expId}->{libraryId (SRXxxx)}->{'commented'}  = commented
    my %rnaSeqAnnotations;

    for my $line ( 0..$#{$tsv{'libraryId'}} ){
        my $commented = 0;
        if ( ($tsv{'libraryId'}[$line] =~ /^#(.+)/) or ($tsv{'libraryId'}[$line] =~ /^\"#(.+)/)){
            $tsv{'libraryId'}[$line] = $1;
            $commented = 1;
        }

        # file format: libraryId experimentId platform ... organId stageId
# #libraryId        experimentId   platform      organId organName       uberonId        uberonName      stageId stageName       infoOrgan       infoStage       libraryTitle     librarySource    libraryDescription       libraryCharacteristics    organAnnotationStatus   organBiologicalStatus   stageAnnotationStatus   stageBiologicalStatus   sex     strain  speciesId   comment annotatorId     lastModificationDate
        # file format: libraryId experimentId platform ... organId stageId
        my $libraryId    = $tsv{'libraryId'}[$line];
        my $experimentId = $tsv{'experimentId'}[$line];
        my $platform     = $tsv{'platform'}[$line];
        my $uberonId     = $tsv{'uberonId'}[$line];
        my $stageId      = $tsv{'stageId'}[$line];
        my $sex          = $tsv{'sex'}[$line];
        my $strain       = $tsv{'strain'}[$line];
        my $speciesId    = $tsv{'speciesId'}[$line];

        if ( !defined $rnaSeqAnnotations{$experimentId}->{$libraryId} ){
            $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'commented'} = $commented;
            # platform
            if ( $platform ne '' ){
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'platform'} = $platform;
            }
            else {
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'platform'} = '';
                warn "Warning: no platform specified for [$experimentId--$libraryId]. Commented: $commented\n";
            }
            # uberonId
            if ( $uberonId ne '' ){
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'uberonId'} = $uberonId;
            }
            else {
                warn "Warning: no uberonId specified for [$experimentId--$libraryId]. Commented: $commented\n";
            }
            # stageId
            if ( $stageId ne '' ){
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'stageId'} = $stageId;
            }
            else {
                warn "Warning: no stageId specified for [$experimentId--$libraryId]. Commented: $commented\n";
            }
            # sex
            if ( $sex ne '' ){
              # Normalize sex info
              my $norm_sex = $sex eq 'F'         ? 'female'
                           : $sex eq 'M'         ? 'male'
                           : $sex eq 'H'         ? 'hermaphrodite'
                           : $sex eq 'U'         ? 'not annotated'
                           : $sex eq 'mixed'     ? 'mixed'
                           : $sex =~ /^[Mm]ixed/ ? 'mixed' # Mixed 1:1
                           : $sex eq 'NA'        ? 'NA'
                           : 'NA';

                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'sex'} = $norm_sex;
            }
            else {
                warn "Warning: no sex specified for [$experimentId--$libraryId]. Commented: $commented\n";
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'sex'} = 'NA';
            }
            # strain
            if ( $strain ne '' ){
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'strain'} = $strain;
            }
            else {
                warn "Warning: no strain specified for [$experimentId--$libraryId]. Commented: $commented\n";
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'strain'} = 'NA';
            }
            # species
            if ( $speciesId ne '' ){
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'speciesId'} = $speciesId;
            }
            else {
                warn "Warning: no species specified for [$experimentId--$libraryId]. Commented: $commented\n";
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'speciesId'} = 'NA';
            }
        }
        else {
            warn 'Warning: library present several times in the annotation file: experiment: ',
            $experimentId, ' - library: ', $libraryId, "\n";
        }
    }
    return %rnaSeqAnnotations;
}

## TODO remove? This doesn't take care of removing the quotes
sub getAllRnaSeqAnnotations {
    my ($rnaSeqAnnotationFile) = @_;
    # $rnaSeqAnnotations{expId}->{libraryId (SRXxxx)}->{'platform'}   = platform
    # $rnaSeqAnnotations{expId}->{libraryId (SRXxxx)}->{'uberonId'}   = uberonId
    # $rnaSeqAnnotations{expId}->{libraryId (SRXxxx)}->{'stageId'}    = stageId
    # $rnaSeqAnnotations{expId}->{libraryId (SRXxxx)}->{'sex'}        = sex
    # $rnaSeqAnnotations{expId}->{libraryId (SRXxxx)}->{'strain'}     = strain
    # $rnaSeqAnnotations{expId}->{libraryId (SRXxxx)}->{'commented'}  = commented
    my %rnaSeqAnnotations;

    open(my $IN, '<', $rnaSeqAnnotationFile)  or die "Could not read file [$rnaSeqAnnotationFile]\n";
    my $line = <$IN>; #header
    while ( defined ($line = <$IN>) ){
        chomp $line;
        my $commented = 0;
        if ( ($line =~ /^#(.+)/) or ($line =~ /^\"#(.+)/) ){
            $line      = $1;
            $commented = 1;
        }

        # file format: libraryId experimentId platform ... organId stageId
# #libraryId        experimentId   platform      organId organName       uberonId        uberonName      stageId stageName       infoOrgan       infoStage       libraryTitle     librarySource    libraryDescription       libraryCharacteristics    organAnnotationStatus   organBiologicalStatus   stageAnnotationStatus   stageBiologicalStatus   sex     strain  speciesId   comment annotatorId     lastModificationDate
        my @tmp = map { bgeeTrim($_) } split(/\t/, $line);
        # remove quotes
        for my $i (0 .. $#tmp){
          $tmp[$i] =~ s/^\"//;
          $tmp[$i] =~ s/\"$//;
        }

        my $libraryId    = $tmp[0];
        my $experimentId = $tmp[1];
        my $platform     = $tmp[2];
        my $uberonId     = $tmp[5];
        my $stageId      = $tmp[7];
        my $sex          = $tmp[19];
        my $strain       = $tmp[20];
        my $speciesId    = $tmp[21];

        if ( !defined $rnaSeqAnnotations{$experimentId}->{$libraryId} ){
            $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'commented'} = $commented;
            # platform
            if ( $platform ne '' ){
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'platform'} = $platform;
            }
            else {
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'platform'} = '';
                warn "Warning: no platform specified for [$experimentId--$libraryId]. Commented: $commented\n";
            }
            # uberonId
            if ( $uberonId ne '' ){
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'uberonId'} = $uberonId;
            }
            else {
                warn "Warning: no uberonId specified for [$experimentId--$libraryId]. Commented: $commented\n";
            }
            # stageId
            if ( $stageId ne '' ){
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'stageId'} = $stageId;
            }
            else {
                warn "Warning: no stageId specified for [$experimentId--$libraryId]. Commented: $commented\n";
            }
            # sex
            if ( $sex ne '' ){
              # Normalize sex info
              my $norm_sex = $sex eq 'F'         ? 'female'
                           : $sex eq 'M'         ? 'male'
                           : $sex eq 'H'         ? 'hermaphrodite'
                           : $sex eq 'U'         ? 'not annotated'
                           : $sex eq 'mixed'     ? 'mixed'
                           : $sex =~ /^[Mm]ixed/ ? 'mixed' # Mixed 1:1
                           : $sex eq 'NA'        ? 'NA'
                           : 'NA';

                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'sex'} = $norm_sex;
            }
            else {
                warn "Warning: no sex specified for [$experimentId--$libraryId]. Commented: $commented\n";
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'sex'} = 'NA';
            }
            # strain
            if ( $strain ne '' ){
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'strain'} = $strain;
            }
            else {
                warn "Warning: no strain specified for [$experimentId--$libraryId]. Commented: $commented\n";
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'strain'} = 'NA';
            }
            # species
            if ( $speciesId ne '' ){
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'speciesId'} = $speciesId;
            }
            else {
                warn "Warning: no species specified for [$experimentId--$libraryId]. Commented: $commented\n";
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'speciesId'} = 'NA';
            }
        }
        else {
            warn 'Warning: library present several times in the annotation file: experiment: ',
            $experimentId, ' - library: ', $libraryId, "\n";
        }
    }
    close $IN;

    return %rnaSeqAnnotations;
}

# Retrieve information on libraries
sub getAllRnaSeqLibrariesInfo {
    # Updated to be compatible with new Bgee v14 RNA-seq pipeline. This is reading rna_seq_sample_info.txt file
    my ($rnaSeqLibraryFile) = @_;
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'speciesId'} = ...
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'organism'} = ...
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'genomeFilePath'} = ...
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'database'} = ...
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'platform'} = ...
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'libraryType'} = ...
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'libraryInfo'} = ...
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'readLength'} = ...
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'runIds'}->{$runId} = ()

    my @valid_platforms = ('Illumina', 'NextSeq');

    my %rnaSeqLibrary;
    for my $line ( read_file("$rnaSeqLibraryFile", chomp=>1) ){
        next  if ( $line =~ /^#/ or $line =~ /^\"#/ );
        #libraryId      experimentId    speciesId       organism        genomeFilePath  database        platform        libraryType     libraryInfo     readLength      runIds
        my @tmp = map { bgeeTrim($_) } split(/\t/, $line);

        my $experimentId                    = $tmp[1];
        my $libraryId                       = $tmp[0];
        my $speciesId                       = $tmp[2];
        my $organism                        = $tmp[3];
        my $genomeFilePath                  = $tmp[4];
        my $database                        = $tmp[5];
        my $platform                        = $tmp[6];
        my $libraryType                     = $tmp[7];
        my $libraryInfo                     = $tmp[8];
        my $readLength                      = $tmp[9];
        my $runIds                          = $tmp[10];
        die "tsv field number problem [$line]\n"  if ( scalar @tmp != 11 );

        if ( !defined $rnaSeqLibrary{$experimentId}->{$libraryId} ){
            # Perform format checks
            my $discarded = 0;
            if ( $experimentId eq '' ){
                warn "Warning, wrong format for experimentId [$experimentId]\n";
                $discarded = 1;
            }
            if ( $libraryId eq '' ){
                warn "Warning, wrong format for libraryId [$libraryId]\n";
                $discarded = 1;
            }
            if ( $speciesId eq '' ){
                warn "Warning, wrong format for speciesId [$speciesId]\n";
                $discarded = 1;
            }
            if ( $organism eq '' ){
                warn "Warning, wrong format for organism [$organism]\n";
                $discarded = 1;
            }
            if ( $genomeFilePath eq '' ){
                warn "Warning, wrong format for genomeFilePath [$genomeFilePath]\n";
                $discarded = 1;
            }
            if ( $database ne 'Ensembl' && $database ne 'EnsemblMetazoa' ){
                warn "Warning, wrong format for database [$database]\n";
                $discarded = 1;
            }
            if ( $platform eq '' || all { $platform !~ /^$_/ } @valid_platforms ){
                warn "Warning, wrong format for platform [$platform]\n";
                $discarded = 1;
            }
            if ( $libraryType ne 'SINGLE' && $libraryType ne 'PAIRED' ){
                warn "Warning, wrong format for libraryType [$libraryType]\n";
                $discarded = 1;
            }
            if ( $readLength ne '' && ($readLength !~ /$integerRegex/ || $readLength < 0 ) ){
                warn "Warning, wrong format for readLength [$readLength]\n";
                $discarded = 1;
            }

            my @runs = split(/,/, $runIds);
            my %runs;
            foreach my $runId ( @runs ){
                if ( $runId eq '' ){
                    warn "Warning, wrong format for one runId [$runId]\n";
                    $discarded = 1;
                }
                $runs{$runId} = ();
            }

            if ( $discarded ){
                warn ' for experiment: ', $experimentId, ' - sample: ', $libraryId, ", sample discarded!\n";
            }
            else {
                $rnaSeqLibrary{$experimentId}->{$libraryId}->{'speciesId'}      = $speciesId;
                $rnaSeqLibrary{$experimentId}->{$libraryId}->{'organism'}       = $organism;
                $rnaSeqLibrary{$experimentId}->{$libraryId}->{'genomeFilePath'} = $genomeFilePath;
                $rnaSeqLibrary{$experimentId}->{$libraryId}->{'database'}       = $database;
                $rnaSeqLibrary{$experimentId}->{$libraryId}->{'platform'}       = $platform;
                $rnaSeqLibrary{$experimentId}->{$libraryId}->{'libraryType'}    = lc($libraryType);
                $rnaSeqLibrary{$experimentId}->{$libraryId}->{'libraryInfo'}    = $libraryInfo;
                $rnaSeqLibrary{$experimentId}->{$libraryId}->{'readLength'}     = $readLength;
                foreach my $runId ( keys %runs ){
                    $rnaSeqLibrary{$experimentId}->{$libraryId}->{'runIds'}->{$runId} = ();
                }
            }
        }
        else {
            warn 'Warning: sample present several times in the library file: experiment: ', $experimentId, ' - sample: ', $libraryId, "\n";
        }
    }

    return %rnaSeqLibrary;
}

sub getExcludedLibraries {
    # record libraries excluded because of low mapping quality or problem during mapping
    my ($excludedLibraryFile) = @_;
    my %excludedLibraries;

    for my $line ( read_file("$excludedLibraryFile", chomp=>1) ){
        next  if ( $line =~ /^#/ or $line =~ /^\"#/ );
        my @tmp = map { bgeeTrim($_) } split(/\t/, $line);
        if ( $tmp[1] eq 'TRUE' ){
            $excludedLibraries{$tmp[0]} = ();
        }
    }

    return %excludedLibraries;
}

## function to read cutoff and percent present infos
sub getAllRnaSeqLibrariesStats {
    # Updated to be compatible with new Bgee v14 RNA-seq pipeline. This is reading presence_absence_all_samples.txt file
    my ($rnaSeqLibraryFile) = @_;
    my %rnaSeqLibraries;

    for my $line ( read_file("$rnaSeqLibraryFile", chomp=>1) ){
        next  if ( $line =~ /^#/ or $line =~ /^\"#/ );
        #libraryId      max_intergenic  cutoffTPM       cutoffGenicTPM  cutoffFPKM      cutoffGenicFPKM proportionGenicPresent  numberGenicPresent      numberGenic     proportionCodingPresent numberPresentCoding     numberCoding    proportionIntergenicPresent     numberIntergenicPresent numberIntergenic        ratioIntergenicCodingPresent    species organism

        my @tmp = map { bgeeTrim($_) } split(/\t/, $line);

        my $libraryId                       = $tmp[0];
        my $cutoffTPM                       = $tmp[3];
        my $cutoffFPKM                      = $tmp[5];
        my $allGenesPercentPresent          = $tmp[6];
        my $proteinCodingPercentPresent     = $tmp[9];
        my $intergenicRegionsPercentPresent = $tmp[12];
        my $ratioIntergenicCodingPresent    = $tmp[15];
        my $speciesId                       = $tmp[16];
        my $organism                        = $tmp[17];

        die "tsv field number problem [$line]\n"  if ( scalar @tmp != 18 );

        if ( !defined $rnaSeqLibraries{$libraryId} ){
            # Perform format checks
            my $discarded = 0;
            if ( $libraryId eq '' ){
                warn "Warning, wrong format for libraryId [$libraryId]\n";
                $discarded = 1;
            }
            if ( $cutoffTPM !~ /$floatingPointRegex/ ){
                warn "Warning, wrong format for cutoffTPM [$cutoffTPM]\n";
                $discarded = 1;
            }
            if ( $cutoffFPKM !~ /$floatingPointRegex/ ){
                warn "Warning, wrong format for cutoffFPKM [$cutoffFPKM]\n";
                $discarded = 1;
            }
            if ( $allGenesPercentPresent !~ /$floatingPointRegex/ || $allGenesPercentPresent < 0 || $allGenesPercentPresent > 100 ){
                warn "Warning, wrong format for allGenesPercentPresent [$allGenesPercentPresent]\n";
                $discarded = 1;
            }
            if ( $proteinCodingPercentPresent !~ /$floatingPointRegex/ || $proteinCodingPercentPresent < 0 || $proteinCodingPercentPresent > 100 ){
                warn "Warning, wrong format for proteinCodingPercentPresent [$proteinCodingPercentPresent]\n";
                $discarded = 1;
            }
            if ( $intergenicRegionsPercentPresent !~ /$floatingPointRegex/ || $intergenicRegionsPercentPresent < 0 || $intergenicRegionsPercentPresent > 100 ){
                warn "Warning, wrong format for intergenicRegionsPercentPresent [$intergenicRegionsPercentPresent]\n";
                $discarded = 1;
            }
            if ( $ratioIntergenicCodingPresent !~ /$floatingPointRegex/ || $ratioIntergenicCodingPresent < 0 || $ratioIntergenicCodingPresent > 1 ){
                warn "Warning, wrong format for ratioIntergenicCodingPresent [$ratioIntergenicCodingPresent]\n";
                $discarded = 1;
            }
            if ( $speciesId eq '' ){
                warn "Warning, wrong format for speciesId [$speciesId]\n";
                $discarded = 1;
            }
            if ( $organism eq '' ){
                warn "Warning, wrong format for organism [$organism]\n";
                $discarded = 1;
            }

            if ( $discarded ){
                warn 'Sample ', $libraryId, " discarded!\n";
            }
            else {
                $rnaSeqLibraries{$libraryId}->{'cutoffTPM'}                       = $cutoffTPM;
                $rnaSeqLibraries{$libraryId}->{'cutoffFPKM'}                      = $cutoffFPKM;
                $rnaSeqLibraries{$libraryId}->{'allGenesPercentPresent'}          = $allGenesPercentPresent;
                $rnaSeqLibraries{$libraryId}->{'proteinCodingPercentPresent'}     = $proteinCodingPercentPresent;
                $rnaSeqLibraries{$libraryId}->{'intergenicRegionsPercentPresent'} = $intergenicRegionsPercentPresent;
                # needs to be between 0 and 100, rounded to 2 decimals
                $rnaSeqLibraries{$libraryId}->{'ratioIntergenicCodingPresent'}    = int($ratioIntergenicCodingPresent * 10000 + 0.5) / 100;
                $rnaSeqLibraries{$libraryId}->{'speciesId'}                       = $speciesId;
                $rnaSeqLibraries{$libraryId}->{'organism'}                        = $organism;
            }
        }
        else {
            warn 'Warning: sample present several times in the file: ', $libraryId, "\n";
        }
    }

    return %rnaSeqLibraries;
}

## function to read percent mapping and read length infos
sub getAllRnaSeqReportInfo {
    my ($rnaSeqReportFile) = @_;
    my %rnaSeqLibraries;

    for my $line ( read_file("$rnaSeqReportFile", chomp=>1) ){
        next  if ( $line =~ /^#/ or $line =~ /^\"#/ ); # this includes the header

        my @tmp = map { bgeeTrim($_) } split(/\t/, $line);

        # columns: libraryId, allReadsCount, mappedReadsCount, minReadLength, maxReadLength

        my $libraryId                       = $tmp[0];
        my $allReadsCount                   = $tmp[1];
        my $mappedReadsCount                = $tmp[2];
        my $minReadLength                   = $tmp[3];
        my $maxReadLength                   = $tmp[4];

        die "tsv field number problem [$line]\n"  if ( scalar @tmp != 5 );

        if ( !defined $rnaSeqLibraries{$libraryId} ){
            # Perform format checks
            my $discarded = 0;
            if ( $libraryId eq '' ){
                warn "Warning, wrong format for libraryId [$libraryId]\n";
                $discarded = 1;
            }
            if ( $allReadsCount !~ /$positiveIntegerRegex/ and $allReadsCount ne 'NULL' ){
                warn "Warning, wrong format for allReadsCount [$allReadsCount]\n";
                $discarded = 1;
            }
            if ( $mappedReadsCount !~ /$positiveIntegerRegex/ and $mappedReadsCount ne 'NULL' ){
                warn "Warning, wrong format for mappedReadsCount [$mappedReadsCount]\n";
                $discarded = 1;
            }
            if ( $minReadLength !~ /$positiveIntegerRegex/ and $minReadLength ne 'NULL' ){
                warn "Warning, wrong format for minReadLength [$minReadLength]\n";
                $discarded = 1;
            }
            if ( $maxReadLength !~ /$positiveIntegerRegex/ and $maxReadLength ne 'NULL' ){
                warn "Warning, wrong format for maxReadLength [$maxReadLength]\n";
                $discarded = 1;
            }

            if ( $discarded ){
                warn 'Sample ', $libraryId, " discarded!\n";
            }
            else {
                $rnaSeqLibraries{$libraryId}->{'allReadsCount'}    = $allReadsCount;
                $rnaSeqLibraries{$libraryId}->{'mappedReadsCount'} = $mappedReadsCount;
                if ( $minReadLength eq 'NULL' ){
                    $rnaSeqLibraries{$libraryId}->{'minReadLength'}  = 0;
                }
                else {
                    $rnaSeqLibraries{$libraryId}->{'minReadLength'}  = $minReadLength;
                }
                if ( $maxReadLength eq 'NULL' ){
                    $rnaSeqLibraries{$libraryId}->{'maxReadLength'}  = 0;
                }
                else {
                    $rnaSeqLibraries{$libraryId}->{'maxReadLength'}  = $maxReadLength;
                }
            }
        }
        else {
            warn 'Warning: sample present several times in the file: ', $libraryId, "\n";
        }
    }
    return %rnaSeqLibraries;
}

# Retrieve for a sample the expression calls, TPM and FPKM values for all genes
sub getGenesResults {
    my ($sampleExpCallsFile) = @_;
    # $expressionCalls{geneId}->{'estimatedCount'} = ...
    # $expressionCalls{geneId}->{'FPKM'} = ...
    # $expressionCalls{geneId}->{'TPM'} = ...
    # $expressionCalls{geneId}->{'biotype'} = ...
    # $expressionCalls{geneId}->{'expressionCall'} = ...
    my %expressionCalls;

    open(my $IN, '<', $sampleExpCallsFile)  or die "Could not read file [$sampleExpCallsFile]\n";
    my $line = <$IN>;    #header
    while ( defined ($line = <$IN>) ){
        chomp $line;
        # file format: geneId, estimatedCount, FPKM, TPM, biotype, expression call
        my @tmp = map { bgeeTrim($_) } split(/\t/, $line);
        my $geneId         = $tmp[0];
        my $estimatedCount = $tmp[1];
        my $TPM            = $tmp[2];
        my $FPKM           = $tmp[3];
        my $biotype        = $tmp[4];
        my $expressionCall = $tmp[5];

        if ( !defined $expressionCalls{$geneId} ){
            # Perform format checks
            my $discarded = 0;
            if ( $geneId eq '' || $geneId =~ /\s/ ){
                warn "Warning, wrong format for geneId [$geneId]\n";
                $discarded = 1;
            }
            if ( $estimatedCount !~ /$floatingPointRegex/ || $estimatedCount < 0 ){
                warn "Warning, wrong format for estimatedCount [$estimatedCount]\n";
                $discarded = 1;
            }
            if ( $FPKM !~ /$floatingPointRegex/ ){
                warn "Warning, wrong format for FPKM [$FPKM]\n";
                $discarded = 1;
            }
            if ( $TPM !~ /$floatingPointRegex/ ){
                warn "Warning, wrong format for TPM [$TPM]\n";
                $discarded = 1;
            }
            if ( $biotype eq '' ){
                warn "Warning, wrong format for biotype [$biotype]\n";
                $discarded = 1;
            }
            if ( $expressionCall ne 'absent' && $expressionCall ne 'present' ){
                warn "Warning, wrong format for expressionCall [$expressionCall]\n";
                $discarded = 1;
            }

            if ( $discarded ){
                warn ' for gene: ', $geneId, "\n";
            }
            else {
                $expressionCalls{$geneId}->{'estimatedCount'} = $estimatedCount;
                $expressionCalls{$geneId}->{'FPKM'}           = $FPKM;
                $expressionCalls{$geneId}->{'TPM'}            = $TPM;
                $expressionCalls{$geneId}->{'expressionCall'} = $expressionCall;
                $expressionCalls{$geneId}->{'biotype'}        = $biotype;
            }
        }
        else {
            warn "Warning: expression call defined several times for geneId [$geneId]\n";
        }
    }
    close $IN;

    return %expressionCalls;
}

## TODO keep this?
# Julien Roux, Oct 2015
# Getting feature length used for FPKM calculation
# We retrieve feature length for genes of any species and any biotype
sub getFeatureLength {
    my ($featureLengthFolder, $ensRelease) = @_;
    my %featureLength;
    for my $file ( glob($featureLengthFolder.'/*'.$ensRelease.'.gtf_all_length') ){
        $file = basename($file);
        open(my $IN, '<', $featureLengthFolder.'/'.$file)  or die "Could not read file [$file]\n";
        while ( defined (my $line = <$IN>) ) {
            chomp $line;
            # file format: geneId, featureLength, biotype
            my @tmp = map { bgeeTrim($_) } split(/\t/, $line);
            my $geneId         = $tmp[0];
            my $featureLength  = $tmp[1];
            my $biotype        = $tmp[2];

            $featureLength{$geneId}->{'featureLength'} = $featureLength;
            $featureLength{$geneId}->{'biotype'}       = $biotype;
        }
        close $IN;
    }
    return %featureLength;
}

1;

