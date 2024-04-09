#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;
use FindBin;
use lib "$FindBin::Bin/.."; # Get lib path for Utils.pm
use Utils;
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
        my $experimentId = $tsv{'experimentId'}[$line];
        if ( ($tsv{'experimentId'}[$line] =~ /^#(.+)/) or ($tsv{'experimentId'}[$line] =~ /^\"#(.+)/) ){
            $tsv{'experimentId'}[$line] = $1;
            $commented = 1;
        }

        if ( length($tsv{'experimentId'}[$line]) < 255 && $tsv{'experimentId'}[$line] !~ /\s/ ){
            # Initialize the hash for this expId
            $experiments{$experimentId}->{'name'}        = '';
            $experiments{$experimentId}->{'description'} = '';
            $experiments{$experimentId}->{'source'}      = '';
            $experiments{$experimentId}->{'status'}      = '';
            $experiments{$experimentId}->{'commented'}   = '';
        }
        else {
            warn "Badly formatted RNAseqExperiment file at line [$line]: [invalid experimentId]\n";
            next;
        }
        if ($tsv{'experimentName'}[$line] eq '') {
            warn "experiment Name is emtpy for $experimentId";
        }
        if ($tsv{'experimentDescription'}[$line] eq '') {
            warn "experiment description is emtpy for $experimentId";
        }
        if ($tsv{'experimentSource'}[$line] eq '') {
            warn "experiment source is emtpy for $experimentId";
        }
        if ($tsv{'experimentStatus'}[$line] eq '') {
            warn "experiment status is emtpy for $experimentId";
        }
        $experiments{$experimentId}->{'name'}        = $tsv{'experimentName'}[$line];
        $experiments{$experimentId}->{'description'} = $tsv{'experimentDescription'}[$line];
        $experiments{$experimentId}->{'source'}      = $tsv{'experimentSource'}[$line];
        $experiments{$experimentId}->{'status'}      = $tsv{'experimentStatus'}[$line];
        $experiments{$experimentId}->{'commented'}   = $commented;
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
        ($tmp[0], $tmp[1]) = map { bgeeTrim($_) } map { s/^\"//; s/\"$//; $_ } split(/\t/, $line);
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
    # $rnaSeqAnnotations{expId}->{libraryId (SRXxxx)}->{'anatId'}     = anatId
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
# #libraryId        experimentId   platform      organId organName      anatId        anatName      stageId stageName       infoOrgan       infoStage       libraryTitle     librarySource    libraryDescription       libraryCharacteristics    organAnnotationStatus   organBiologicalStatus   stageAnnotationStatus   stageBiologicalStatus   sex     strain  speciesId   comment annotatorId     lastModificationDate
        # file format: libraryId experimentId platform ... organId stageId
        my $libraryId    = $tsv{'libraryId'}[$line];
        my $experimentId = $tsv{'experimentId'}[$line];
        my $platform     = $tsv{'platform'}[$line];
        my $anatId     = $tsv{'anatId'}[$line];
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
            # anatId
            if ( $anatId ne '' ){
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'anatId'} = $anatId;
            }
            else {
                warn "Warning: no anatId specified for [$experimentId--$libraryId]. Commented: $commented\n";
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
            $experimentId, ' - library: ', $libraryId, ". Commented: $commented\n";
        }
    }
    return %rnaSeqAnnotations;
}

## TODO remove? This doesn't take care of removing the quotes
sub getAllRnaSeqAnnotations {
    my ($rnaSeqAnnotationFile) = @_;
    # $rnaSeqAnnotations{expId}->{libraryId (SRXxxx)}->{'platform'}   = platform
    # $rnaSeqAnnotations{expId}->{libraryId (SRXxxx)}->{'anatId'}   = anatId
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

        # file format: 
        # #libraryId    experimentId  platform  organId   organName anatId  uberonName    stageId   stageName infoOrgan infoStage sampleTitle   sampleSource  sampleDescription sampleCharacteristics organAnnotationStatus organBiologicalStatus stageAnnotationStatus stageBiologicalStatus sex   strain    speciesId comment   annotatorId   lastModificationDate  replicate infoReplicate SRSId tags  RNASeqProtocol    physiological status  globin_reduction  PATOid    PATOname
        my @tmp = map { bgeeTrim($_) }
                  map { s/^\"//; s/\"$//; $_ } # remove quotes
                  split(/\t/, $line);

        my $libraryId     = $tmp[0];
        my $experimentId  = $tmp[1];
        my $platform      = $tmp[2];
        my $anatId        = $tmp[4];
        my $stageId       = $tmp[6];
        my $freeTextOrgan = $tmp[9];
        my $freeTextStage = $tmp[10];
        my $sex           = $tmp[14];
        my $strain        = $tmp[15];
        my $genotype      = $tmp[16];
        my $speciesId     = $tmp[17];
        my $protocol      = $tmp[18];
        my $protocolType  = $tmp[19];
        my $popCapture    = $tmp[20];

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
            # anatId
            if ( $anatId ne '' ){
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'anatId'} = $anatId;
            }
            else {
                warn "Warning: no anatId specified for [$experimentId--$libraryId]. Commented: $commented\n";
            }
            # stageId
            if ( $stageId ne '' ){
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'stageId'} = $stageId;
            }
            else {
                warn "Warning: no stageId specified for [$experimentId--$libraryId]. Commented: $commented\n";
            }
            # free text author annotation. no warning throws if empty
            $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'freeTextOrgan'} = $freeTextOrgan;
            $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'freeTextStage'} = $freeTextStage;
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
            # genotype information. Do not throw a warning if empty
            $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'genotype'} = $genotype;
            # species
            if ( $speciesId ne '' ){
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'speciesId'} = $speciesId;
            }
            else {
                warn "Warning: no species specified for [$experimentId--$libraryId]. Commented: $commented\n";
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'speciesId'} = 'NA';
            }
            # protocol and protocol type are often empty. Do not throw a warning if empty
            $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'protocol'} = $protocol;
            $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'protocolType'} = $protocolType;
            # population capture
            if ( !defined($popCapture) || $popCapture ne '' ){
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'populationCapture'} = $popCapture;
            }
            else {
                $rnaSeqAnnotations{$experimentId}->{$libraryId}->{'populationCapture'} = '';
                warn "Warning: no RNA population captured specified for [$experimentId--$libraryId]. Commented: $commented\n";
            }
        }
        else {
            warn 'Warning: library present several times in the annotation file: experiment: ',
            $experimentId, ' - library: ', $libraryId, ". Commented: $commented\n";
        }
    }
    close $IN;

    return %rnaSeqAnnotations;
}

# Retrieve information on libraries
sub getAllRnaSeqLibrariesInfo {
    # Updated to be compatible with new Bgee v15 RNA-seq pipeline. This is reading rna_seq_sample_info.txt file
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
        my @tmp = map { bgeeTrim($_) } map { s/^\"//; s/\"$//; $_ } split(/\t/, $line);

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
            if ( $database ne 'Ensembl' && $database ne 'EnsemblMetazoa' && $database ne 'RefSeq'){
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
            #NOTE SRA may return float for read length because of read length variability
            if ( $readLength ne '' && ($readLength !~ /$floatingPointRegex/ || $readLength < 0 ) ){
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
        next  if ( $line =~ /^#/ or $line =~ /^\"#/ or $line =~ /^$/ );
        my @tmp = map { bgeeTrim($_) } map { s/^\"//; s/\"$//; $_ } split(/\t/, $line);
        if ( $tmp[1] eq 'TRUE' ){
            $excludedLibraries{$tmp[0]} = $tmp[2];
        }
    }

    return %excludedLibraries;
}

## function to read cutoff and percent present infos
sub getAllRnaSeqLibrariesStats {
    # Updated to be compatible with new Bgee v15 RNA-seq pipeline. This is reading presence_absence_all_samples.txt file
    my ($rnaSeqLibraryFile) = @_;
    my %rnaSeqLibraries;

    for my $line ( read_file("$rnaSeqLibraryFile", chomp=>1) ){
        next  if ( $line =~ /^#/ or $line =~ /^\"#/ or $line =~ /^libraryId/);
        #libraryId    cutoffTPM    proportionGenicPresent  numberGenicPresent      numberGenic     proportionCodingPresent numberPresentCoding     numberCoding    proportionIntergenicPresent     numberIntergenicPresent numberIntergenic    pValueCutoff    meanIntergenic  sdIntergenic    speciesId

        my @tmp = map { bgeeTrim($_) } map { s/^\"//; s/\"$//; $_ } split(/\t/, $line);

        my $libraryId                       = $tmp[0];
        my $cutoffTPM                       = $tmp[1];
        my $allGenesPercentPresent          = $tmp[2];
        my $proteinCodingPercentPresent     = $tmp[5];
        my $intergenicRegionsPercentPresent = $tmp[8];
        my $pValueThreshold                 = $tmp[11];
        my $meanIntergenic                  = $tmp[12];
        my $sdIntergenic                    = $tmp[13];
        my $speciesId                       = $tmp[14];

        die "tsv field number problem [$line]\n"  if ( scalar @tmp != 15 );

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
            if ( $pValueThreshold !~ /$floatingPointRegex/ || $pValueThreshold < 0 || $pValueThreshold > 1 ){
                warn "Warning, wrong format for pValueThreshold [$pValueThreshold]\n";
                $discarded = 1;
            }
            if ( $meanIntergenic !~ /$floatingPointRegex/ || $meanIntergenic < 0){
                warn "Warning, wrong format for meanIntergenic [$meanIntergenic]\n";
                $discarded = 1;
            }
            if ( $sdIntergenic !~ /$floatingPointRegex/ || $sdIntergenic < 0){
                warn "Warning, wrong format for sdIntergenic [$sdIntergenic]\n";
                $discarded = 1;
            }
            if ( $speciesId eq '' ){
                warn "Warning, wrong format for speciesId [$speciesId]\n";
                $discarded = 1;
            }

            if ( $discarded ){
                warn 'Sample ', $libraryId, " discarded!\n";
            }
            else {
                $rnaSeqLibraries{$libraryId}->{'cutoffTPM'}                       = $cutoffTPM;
                $rnaSeqLibraries{$libraryId}->{'allGenesPercentPresent'}          = $allGenesPercentPresent;
                $rnaSeqLibraries{$libraryId}->{'proteinCodingPercentPresent'}     = $proteinCodingPercentPresent;
                $rnaSeqLibraries{$libraryId}->{'intergenicRegionsPercentPresent'} = $intergenicRegionsPercentPresent;
                $rnaSeqLibraries{$libraryId}->{'pValueThreshold'}                 = $pValueThreshold;
                $rnaSeqLibraries{$libraryId}->{'meanIntergenic'}                  = $meanIntergenic;
                $rnaSeqLibraries{$libraryId}->{'sdIntergenic'}                    = $sdIntergenic;
                $rnaSeqLibraries{$libraryId}->{'speciesId'}                       = $speciesId;

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
        next  if ( $line =~ /^#/ or $line =~ /^\"#/ or $line =~ /^libraryId/); # this includes the header

        my @tmp = map { bgeeTrim($_) } map { s/^\"//; s/\"$//; $_ } split(/\t/, $line);

        # columns: libraryId    reads   min_read_size   max_read_size   number_aligned  number_unique_aligned   prop_aligned    prop_unique_aligned number_targets  start_time  kallisto_version

        my $libraryId                       = $tmp[0];
        my $allReadsCount                   = $tmp[1];
        #transform from potential scientific notation to integer
        if ( $allReadsCount =~ /e[+|-][0-9]+$/ ) {
            $allReadsCount = sprintf("%.0f", $tmp[1]);
        }
        my $minReadLength                   = $tmp[2];
        my $maxReadLength                   = $tmp[3];
        my $mappedReadsCount                = $tmp[4];
        #transform from potential scientific notation to integer
        if ( $allReadsCount =~ /^[0-9]+e[+|-][0-9]+$/ ) {
            $allReadsCount = sprintf("%.0f", $tmp[4]);
        }


        die "tsv field number problem [$line]\n"  if ( scalar @tmp != 11 );

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

# Retrieve for a sample the expression calls, TPM values for all genes
sub getGenesResults {
    my ($sampleExpCallsFile) = @_;
    # $expressionCalls{geneId}->{'estimatedCount'} = ...
    # $expressionCalls{geneId}->{'FPKM'} = ...
    # $expressionCalls{geneId}->{'TPM'} = ...
    # $expressionCalls{$geneId}->{'zscore'} = ...
    # $expressionCalls{$geneId}->{'pValue'} = ...
    # $expressionCalls{geneId}->{'biotype'} = ...
    # $expressionCalls{geneId}->{'expressionCall'} = ...
    my %expressionCalls;

    open(my $IN, '<', $sampleExpCallsFile)  or die "Could not read file [$sampleExpCallsFile]\n";
    my $line = <$IN>;    #header
    while ( defined ($line = <$IN>) ){
        chomp $line;
        # file format: geneId, tpm, counts, length, biotype, type, zscore, pvalue, expression call, fpkm

        my @tmp = map { bgeeTrim($_) } map { s/^\"//; s/\"$//; $_ } split(/\t/, $line);
        my $geneId         = $tmp[0];
        my $TPM            = $tmp[1];
        my $estimatedCount = $tmp[2];
        my $biotype        = $tmp[4];
        my $zscore         = $tmp[6];
        my $pValue         = $tmp[7];
        my $expressionCall = $tmp[8];

        if ( !defined $expressionCalls{$geneId} ){
            # Perform format checks
            my $discarded = 0;
            if ( $geneId eq '' || $geneId =~ /\s/ ){
                warn "Warning, wrong format for geneId [$geneId]\n";
                $discarded = 1;
            }
            if ( $TPM !~ /$floatingPointRegex/ || $TPM < 0 || $TPM > 1000000){
                warn "Warning, wrong format for TPM [$TPM]\n";
                $discarded = 1;
            }
            if ( $estimatedCount !~ /$floatingPointRegex/ || $estimatedCount < 0 ){
                warn "Warning, wrong format for estimatedCount [$estimatedCount]\n";
                $discarded = 1;
            }
            if ( $biotype eq '' ){
                warn "Warning, wrong format for biotype [$biotype]\n";
                $discarded = 1;
            }
             if ( !($zscore =~ /$floatingPointRegex/ || $zscore eq 'NA') ){
                warn "Warning, wrong format for zscore [$zscore]\n";
                $discarded = 1;
            }
            if ( !( $pValue eq 'NA' || ($pValue =~ /$floatingPointRegex/ && $pValue >= 0 && $pValue <= 1) ) ){
                warn "Warning, wrong format for pValue [$pValue]\n";
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
                $expressionCalls{$geneId}->{'TPM'}            = $TPM;
                $expressionCalls{$geneId}->{'estimatedCount'} = $estimatedCount;
                $expressionCalls{$geneId}->{'biotype'}        = $biotype;
                $expressionCalls{$geneId}->{'zscore'}         = $zscore;
                $expressionCalls{$geneId}->{'pValue'}         = $pValue;
                $expressionCalls{$geneId}->{'expressionCall'} = $expressionCall;
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
            my @tmp = map { bgeeTrim($_) } map { s/^\"//; s/\"$//; $_ } split(/\t/, $line);
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

# Julien Wollbrett, Feb. 2021
# Read file containing mapping between RNA-Seq protocols and gene biotypes
# for which absent calls does not have to be created
# return a hash with a protocol name as key and an array of biotypes as value
sub retrieveProtocolsToBiotypeExcludeAbsentCalls {
    my %protocolToBiotypes = ();
    my ($mappingProtocolToBiotypesFile) = @_;
    for my $line ( read_file("$mappingProtocolToBiotypesFile", chomp=>1) ){
        #skip the header
        next  if ( $line =~ /^RNASeqProtocol/ or $line =~ /^\"RNASeqProtocol/ );
        # RNASeqProtocol    biotypes_excluded_for_absent_calls
        my @columns = map { bgeeTrim($_) } map { s/^\"//; s/\"$//; $_ } split(/\t/, $line);
        my @biotypes = split(/,/, $columns[1]);

        $protocolToBiotypes{$columns[0]} = \@biotypes;
    }
    return %protocolToBiotypes;
}

######################################################################################################
######################################################################################################
################################### FULL LENGTH scRNASEQ FUNCTIONS ###################################
######################################################################################################
######################################################################################################

# TODO: homogenize columns of sample_info file from bulk and single cell in order to get libraries info from the same function
# Retrieve information on libraries of full length RNASeq
sub getAllFullLengthScRnaSeqLibrariesInfo {
    my ($fullLengthScRnaSeqLibraryFile) = @_;
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'speciesId'} = ...
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'organism'} = ...
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'platform'} = ...
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'libraryType'} = ...
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'libraryInfo'} = ...
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'readLength'} = ...
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'anatId'} = ...
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'stageId'} = ...
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'cellTypeId'} = ...
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'sex'} = ...
    # $rnaSeqLibrary{experimentId}->{libraryId (SRX...)}->{'strain'} = ...

    my @valid_platforms = ('Illumina HiSeq 2500', 'Illumina HiSeq 2000', 'Illumina MiSeq', 'NextSeq 500');

    my @valid_protocols = ('Adapted_SMART_seq2', 'Fluidigm C1 instrument and Nextera XT protocol',
        'NEBNext® Ultra™ DNA Library Prep Kit for Illumina®', 'SMARTer Ultra Low', 'SMART_seq', 'SMART_seq2');


    my %fullLengthScRnaSeqLibrary;
    for my $line ( read_file("$fullLengthScRnaSeqLibraryFile", chomp=>1) ){
        next  if ( $line =~ /^#/ or $line =~ /^libraryId/ );
        #libraryId  experimentId    cellTypeName    cellTypeId  speciesId   platform    protocol    protocolType    libraryType infoOrgan   stageId uberonId    sex strain  readLength  organism
        my @tmp = map { bgeeTrim($_) } map { s/^\"//; s/\"$//; $_ } split(/\t/, $line);

        my $libraryId                       = $tmp[0];
        my $experimentId                    = $tmp[1];
        my $cellTypeId                      = $tmp[3];
        my $speciesId                       = $tmp[4];
        my $platform                        = $tmp[5];
        my $protocol                        = $tmp[6];
        my $libraryType                     = $tmp[8];
        my $stageId                         = $tmp[10];
        my $uberonId                        = $tmp[11];
        my $sex                             = $tmp[12];
        my $strain                          = $tmp[13];
        my $readLength                      = $tmp[14];
        my $organism                        = $tmp[15];
        my $genotype                        = $tmp[16];

        #TODO: change the annotation to fit authorized sex values in the DB ('not annotated','hermaphrodite','female','male','mixed','NA')
        # it is ugly to have to manually modify the values in this script
        if($sex eq "(Missing)") {
            $sex = $Utils::NA_SEX;
        }
        if($sex eq "M") {
            $sex = $Utils::MALE_SEX;
        }

        die "tsv field number problem [$line]\n"  if ( scalar @tmp != 16 );

        if ( !defined $fullLengthScRnaSeqLibrary{$experimentId}->{$libraryId} ){
            # Perform format checks
            my $discarded = 0;
            if ( $libraryId eq '' ){
                warn "Warning, wrong format for libraryId [$libraryId]\n";
                $discarded = 1;
            }
            if ( $experimentId eq '' ){
                warn "Warning, wrong format for experimentId [$experimentId]\n";
                $discarded = 1;
            }
            if ( $speciesId eq '' ){
                warn "Warning, wrong format for speciesId [$speciesId]\n";
                $discarded = 1;
            }
            if ( $platform eq '' || all { $platform !~ /^$_/ } @valid_platforms ){
                warn "Warning, wrong format for platform [$platform]\n";
                $discarded = 1;
            }
            if ( $protocol eq '' || all { $protocol !~ /^$_/ } @valid_protocols ){
                warn "Warning, wrong format for protocol [$protocol]\n";
                $discarded = 1;
            }
            if ( $libraryType ne 'SINGLE' && $libraryType ne 'PAIRED' ){
                warn "Warning, wrong format for libraryType [$libraryType]\n";
                $discarded = 1;
            }
            #NOTE SRA may return float for read length because of read length variability
            if ( $readLength ne '' && ($readLength !~ /$floatingPointRegex/ || $readLength < 0 ) ){
                warn "Warning, wrong format for readLength [$readLength]\n";
                $discarded = 1;
            }
            if ( $organism eq '' ){
                warn "Warning, wrong format for organism [$organism]\n";
                $discarded = 1;
            }
            if ( $cellTypeId eq '' ){
                warn "Warning, wrong format for cellTypeId [$cellTypeId]\n";
                $discarded = 1;
            }
            if ( $cellTypeId eq '' ){
                warn "Warning, wrong format for cellTypeId [$cellTypeId]\n";
                $discarded = 1;
            }
            if ( $stageId eq '' ){
                warn "Warning, wrong format for stageId [$stageId]\n";
                $discarded = 1;
            }
            if ( $uberonId eq '' ){
                warn "Warning, wrong format for uberonId [$uberonId]\n";
                $discarded = 1;
            }
            if ( $sex eq '' ){
                warn "Warning, wrong format for sex [$sex]\n";
                $discarded = 1;
            }
            if ( $strain eq '' ){
                warn "Warning, wrong format for strain [$strain]\n";
                $discarded = 1;
            }
            # genotype can be null if no information provided
            if ( $genotype eq '' ){
                $genotype = undef;
            }

            if ( $discarded ){
                warn ' for experiment: ', $experimentId, ' - sample: ', $libraryId, ", sample discarded!\n";
            }
            else {
                $fullLengthScRnaSeqLibrary{$experimentId}->{$libraryId}->{'speciesId'}      = $speciesId;
                $fullLengthScRnaSeqLibrary{$experimentId}->{$libraryId}->{'organism'}       = $organism;
                $fullLengthScRnaSeqLibrary{$experimentId}->{$libraryId}->{'platform'}       = $platform;
                $fullLengthScRnaSeqLibrary{$experimentId}->{$libraryId}->{'protocol'}       = $protocol;
                $fullLengthScRnaSeqLibrary{$experimentId}->{$libraryId}->{'libraryType'}    = lc($libraryType);
                $fullLengthScRnaSeqLibrary{$experimentId}->{$libraryId}->{'readLength'}     = $readLength;
                $fullLengthScRnaSeqLibrary{$experimentId}->{$libraryId}->{'cellTypeId'}     = $cellTypeId;
                $fullLengthScRnaSeqLibrary{$experimentId}->{$libraryId}->{'stageId'}        = $stageId;
                $fullLengthScRnaSeqLibrary{$experimentId}->{$libraryId}->{'uberonId'}       = $uberonId;
                $fullLengthScRnaSeqLibrary{$experimentId}->{$libraryId}->{'sex'}            = $sex;
                $fullLengthScRnaSeqLibrary{$experimentId}->{$libraryId}->{'strain'}         = $strain;
                $fullLengthScRnaSeqLibrary{$experimentId}->{$libraryId}->{'genotype'}       = $genotype;
            }
        }
        else {
            warn 'Warning: sample present several times in the library file: experiment: ', $experimentId, ' - sample: ', $libraryId, "\n";
        }
    }

    return %fullLengthScRnaSeqLibrary;
}

1;

