#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;
use Data::Dumper;

# Retrieve information from target base library annotation files manually created by curators. 
# This file also contains condition parameters.
##XXX The cell type ID information in this file corresponds to the parent cell type to all cell type
##    detected by the authors. The cell type inserted in the bgee database is not the one from this file.
sub getTargetBaseCuratedLibrariesAnnotation {
    my ($targetBaseLibraryFile) = @_;
    # $targetBaseLibrary{experimentId}->{libraryId (SRX...)}->{'platform'} = ...
    # $targetBaseLibrary{experimentId}->{libraryId (SRX...)}->{'anatEntityId'} = ...
    # $targetBaseLibrary{experimentId}->{libraryId (SRX...)}->{'cellTypeId'} = ...
    # $targetBaseLibrary{experimentId}->{libraryId (SRX...)}->{'stageId'} = ...
    # $targetBaseLibrary{experimentId}->{libraryId (SRX...)}->{'sex'} = ...
    # $targetBaseLibrary{experimentId}->{libraryId (SRX...)}->{'strain'} = ...
    # $targetBaseLibrary{experimentId}->{libraryId (SRX...)}->{'genotype'} = ...
    # $targetBaseLibrary{experimentId}->{libraryId (SRX...)}->{'speciesId'} = ...
    # $targetBaseLibrary{experimentId}->{libraryId (SRX...)}->{'protocol'} = ...
    # $targetBaseLibrary{experimentId}->{libraryId (SRX...)}->{'sequencedTranscriptPart'} = ...
    # $targetBaseLibrary{experimentId}->{libraryId (SRX...)}->{'whiteList'} = ...
    # $targetBaseLibrary{experimentId}->{libraryId (SRX...)}->{'sampleName'} = ...

    my @valid_platforms = ('BGISEQ-500', 'HiSeq X Ten', 'Illumina HiSeq 2000',
        'Illumina HiSeq 2500', 'Illumina HiSeq 3000', 'Illumina HiSeq 4000',
        'Illumina HiSeq X Ten', 'Illumina NovaSeq 6000', 'NextSeq 500');

    my @valid_protocols = ('10X Genomics', 'CEL-seq2', 'Drop-seq',
        'inDrop', 'MARS-Seq', 'Microwell-seq');

    my @validCellCompartments = ("scRNA-seq", "Sn-scRNA-seq");

    my %targetBaseLibrary;
    for my $line ( read_file("$targetBaseLibraryFile", chomp=>1) ){
        next  if ( $line =~ /^#/ or $line =~ /^\"#/ );
        # there is currently 32 columns in the target base library annotation file
        my @tmp = map { bgeeTrim($_) } map { s/^\"//; s/\"$//; $_ } split(/\t/, $line);

        my $libraryId                       = $tmp[0];
        my $experimentId                    = $tmp[1];
        my $platform                        = $tmp[2];
        my $anatEntityId                    = $tmp[3];
        my $cellTypeId                      = $tmp[5];
        my $stageId                         = $tmp[8];
        my $sex                             = $tmp[17];
        my $strain                          = $tmp[18];
        my $genotype                        = $tmp[19];
        my $speciesId                       = $tmp[20];
        my $cellCompartment                 = $tmp[27];
        my $protocol                        = $tmp[28];
        my $sequencedTranscriptPart         = $tmp[29];
        my $whiteList                       = $tmp[30];
        my $sampleName                      = $tmp[31];

        die "tsv field number problem [$line]\n"  if ( scalar @tmp != 32 );

        if ( !defined $targetBaseLibrary{$experimentId}->{$libraryId} ){
            # Perform format checks
            # do not check format for genotype and whiteList as they can be empty
            my $discarded = 0;
            if ($libraryId eq '' ){
                warn "Warning, wrong format for libraryId [$libraryId]\n";
                $discarded = 1;
            }
            if ($experimentId eq '' ){
                warn "Warning, wrong format for experimentId [$experimentId]\n";
                $discarded = 1;
            }
            if ($platform eq '' || all { $platform !~ /^$_/ } @valid_platforms ){
                warn "Warning, wrong format for platform [$platform]\n";
                $discarded = 1;
            }
            if ($anatEntityId eq '' ){
                warn "Warning, wrong format for anatEntityId [$anatEntityId]\n";
                $discarded = 1;
            }
            if ($cellTypeId eq '' ){
                warn "Warning, wrong format for cellTypeId [$cellTypeId]\n";
                $discarded = 1;
            }
            if ($stageId eq '' ){
                warn "Warning, wrong format for stageId [$stageId]\n";
                $discarded = 1;
            }
            if ($sex eq '' ){
                warn "Warning, wrong format for sex [$sex]\n";
                $discarded = 1;
            #TODO create hash with all potential sexes
            } else {
                if ($sex eq "F") {
                    $sex = "female";
                } elsif ($sex eq "M") {
                    $sex = "male";
                } else {
                    warn "Warning, wrong format for sex [$sex]\n";
                    $discarded = 1;
                }
            }
            if ($strain eq '' ){
                warn "Warning, wrong format for strain [$strain]\n";
                $discarded = 1;
            }
            if ($speciesId eq '' ){
                warn "Warning, wrong format for speciesId [$speciesId]\n";
                $discarded = 1;
            }
            if ($cellCompartment eq '' || all { $cellCompartment !~ /^$_/ }
                @validCellCompartments ){
                warn "Warning, wrong format for cellCompartment [$cellCompartment]\n";
                $discarded = 1;
            } elsif ($cellCompartment eq "Sn-scRNA-seq") {
                $cellCompartment = "nucleus";
            } elsif ($cellCompartment eq "scRNA-seq") {
                $cellCompartment = "cell";
            }
            if ($protocol eq '' || all { $protocol !~ /^$_/ } @valid_protocols ){
                warn "Warning, wrong format for protocol [$protocol]\n";
                $discarded = 1;
            }
            if ($sequencedTranscriptPart ne '3\'end'){
                warn "Warning, wrong format for sequencedTranscriptPart",
                    "[$sequencedTranscriptPart]\n";
                $discarded = 1;
            } else {
                $sequencedTranscriptPart = '3prime';
            }
            if ($sampleName eq '' ){
                warn "Warning, wrong format for sampleName [$sampleName]\n";
                $discarded = 1;
            }

            if ( $discarded ){
                warn ' for experiment: ', $experimentId, ' - library: ', $libraryId,
                    ", library discarded!\n";
            }
            else {
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'platform'} = $platform;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'anatEntityId'} = $anatEntityId;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'cellTypeId'} = $cellTypeId;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'stageId'} = $stageId;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'sex'} = $sex;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'strain'} = $strain;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'genotype'} = $genotype;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'speciesId'} = $speciesId;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'cellCompartment'} = $cellCompartment;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'protocol'} = $protocol;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'sequencedTranscriptPart'} =
                    $sequencedTranscriptPart;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'whiteList'} = $whiteList;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'sampleName'} = $sampleName;
            }
        }
        else {
            warn 'Warning: sample present several times in the library file: experiment: ',
            $experimentId, ' - sample: ', $libraryId, "\n";
        }
    }

    return %targetBaseLibrary;
}

sub getSingleCellExperiments {
    my ($targetBaseExperimentFile, @experimentTypes) = @_;
    # $experiments{experimentId}->{'name'}
    # $experiments{experimentId}->{'description'}
    # $experiments{experimentId}->{'source'}
    # $experiments{experimentId}->{'status'}
    # $experiments{experimentId}->{'commented'}
    # $experiments{experimentId}->{'protocol'}
    # $experiments{experimentId}->{'experimentType'}

    #@allowedExpTypes = ('3\'end', 'Full-length', 'Full-length and 3\'end'];

    my %experiments;
    for my $line ( read_file("$targetBaseExperimentFile", chomp=>1) ){
        next  if ( $line =~ /^#/ or $line =~ /^\"#/ );
        # there is currently 32 columns in the target base library annotation file
        my @tmp = map { bgeeTrim($_) } map { s/^\"//; s/\"$//; $_ } split(/\t/, $line);

        my $experimentId                    = $tmp[0];
        my $name                            = $tmp[1];
        my $description                     = $tmp[2];
        my $source                          = $tmp[3];
        my $status                          = $tmp[4];
        my $comment                         = $tmp[5];
        my $protocol                        = $tmp[8];
        my $experimentType                  = $tmp[9];


        die "tsv field number problem [$line]\n"  if ( scalar @tmp != 10 );

        if ( !defined $experiments{$experimentId} ){
            # Perform format checks
            # do not check format for genotype and whiteList as they can be empty
            my $discarded = 0;
            if ($experimentId eq '' ){
                warn "Warning, wrong format for experimentId [$experimentId]\n";
                $discarded = 1;
            }
            if ($source eq '' ){
                warn "Warning, wrong format for source [$source]\n";
                $discarded = 1;
            }
            if ($protocol eq '' ){
                warn "Warning, wrong format for protocol [$protocol]\n";
                $discarded = 1;
            }
            if ($experimentType eq ''){
                warn "Warning, wrong format for experiment type [$experimentType]\n";
                $discarded = 1;
            }
            if ( $discarded ){
                warn ' experiment: ', $experimentId, " discarded!\n";
            }
            else {
                if (grep(/^$experimentType$/, @experimentTypes)) {
                    $experiments{$experimentId}->{'name'} = $name;
                    $experiments{$experimentId}->{'description'} = $description;
                    $experiments{$experimentId}->{'source'} = $source;
                    $experiments{$experimentId}->{'status'} = $status;
                    $experiments{$experimentId}->{'comment'} = $comment;
                    $experiments{$experimentId}->{'protocol'} = $protocol;
                    $experiments{$experimentId}->{'experimentType'} = $experimentType;
                }
            }
        } else {
            warn 'Warning: experiment present several times in the experiment file: experiment: ',
            $experimentId, "\n";
        }
    }
    return %experiments;
}


# Extract library information from the annotation generated by the bgee pipeline (metadata_info_10X.txt)
sub get_processed_libraries_info {
    my ($targetBaseLibraryFile) = @_;
    # $experiments{experimentId}->{libraryId}->{runId}->{'readCount'} = ...
    # $experiments{experimentId}->{libraryId}->{'libraryLayout'} = ...

    my %libInfos = ();
    for my $line ( read_file("$targetBaseLibraryFile", chomp=>1) ){
        next  if ( $line =~ /^#/ or $line =~ /^\"#/ );
        # there is currently 10 columns in the metadata_info_10X.txt file
        my @tmp = map { bgeeTrim($_) } map { s/^\"//; s/\"$//; $_ } split(/\t/, $line);

        my $experimentId                    = $tmp[1];
        my $libraryId                       = $tmp[2];
        my $libraryType                     = $tmp[8];

        die "tsv field number problem [$line]\n"  if ( scalar @tmp != 11 );

        if (!defined $libInfos{$experimentId}->{$libraryId}) {
            # Perform format checks
            # do not check format for genotype and whiteList as they can be empty
            my $discarded = 0;
            if ($libraryId eq '' ){
                warn "Warning, wrong format for libraryId [$libraryId]\n";
                $discarded = 1;
            }
            if ($experimentId eq '' ){
                warn "Warning, wrong format for experimentId [$experimentId]\n";
                $discarded = 1;
            }
            else {
                $libInfos{$experimentId}->{$libraryId}->{'libraryType'} =
                    $libraryType;
            }
        } else {
            warn 'Warning: library present several times in the library file: experiment: ',
            $experimentId, ' - library: ', $libraryId, "\n";
        }
    }

    return %libInfos;
}


# Extract the mapping from barcodes to cell types. Also extract cell type Ids
sub getBarcodeToCellType {
    my ($barcodeToCellTypeFile) = @_;
    # $barcodes{libraryId}->{'barcodes'}->{barcode}->{'cellTypeId'} = ...
    # $barcodes{libraryId}->{'barcodes'}->{barcode}->{'annotationStatus'} = ...
    # $barcodes{libraryId}->{cellType} = ()

    my %barcodes = ();
    my %celltypes = ();
    for my $line ( read_file("$barcodeToCellTypeFile", chomp=>1) ){
        next  if ( $line =~ /^barcode/ or $line =~ /^\"barcode/ );
        # there is currently 10 columns in the metadata_info_10X.txt file
        my @tmp = map { bgeeTrim($_) } map { s/^\"//; s/\"$//; $_ } split(/\t/, $line);

        my $barcode                         = $tmp[0];
        my $libraryId                       = $tmp[2];
        my $cellTypeId                      = $tmp[7];
        my $cellTypeAnnotationStatus        = $tmp[9];

        die "tsv field number problem [$line]\n"  if (scalar @tmp != 12);

        if (!defined $barcodes{$libraryId}->{'barcodes'}->{$barcode}) {
            # Perform format checks
            # do not check format for genotype and whiteList as they can be empty
            my $discarded = 0;
            if ($libraryId eq '' ){
                warn "Warning, wrong format for libraryId [$libraryId]\n";
                $discarded = 1;
            }
            if ($barcode eq '' ){
                warn "Warning, wrong format for barcode [$barcode]\n";
                $discarded = 1;
            }
            if ($cellTypeId eq '' ){
                warn "Warning, wrong format for cellTypeId [$cellTypeId]\n";
                $discarded = 1;
            }
            if ($cellTypeAnnotationStatus eq '' ){
                warn "Warning, wrong format for cellTypeId [$cellTypeAnnotationStatus]\n";
                $discarded = 1;
            }
            else {
                $barcodes{$libraryId}->{'barcodes'}->{$barcode}->{'cellTypeId'} =
                    $cellTypeId;
                $barcodes{$libraryId}->{'barcodes'}->{$barcode}->{'annotationStatus'} =
                    $cellTypeAnnotationStatus;
                $barcodes{$libraryId}->{'cellTypes'}->{$cellTypeId} = ();
            }
        } else {
            warn 'Warning: barcode described several times in the library : library: ',
            $libraryId, ' - barcode: ', $barcode, "\n";
        }
    }

    return %barcodes;
}
# insert one rnaseq library annotated sample and retrieve the value
# of internal autoincrement rnaSeqLibraryAnnotatedSampleId
sub insert_get_annotated_sample {
    my ($insAnnotatedSample, $selectAnnotatedSampleId, $conditionId,
        $libraryId, $debug) = @_;
    my $annotatedSampleId;
    #insert annotated sample
    if ($debug) {
        print 'INSERT INTO rnaSeqLibraryAnnotatedSample: ', $conditionId, ' - ',
                    $libraryId, "\n";
    } else {
        $insAnnotatedSample->execute($libraryId, $conditionId)
            or die $insAnnotatedSample->errstr;
    }
    #retrieve annotated sample ID
    $selectAnnotatedSampleId->execute($conditionId, $libraryId)
        or die $selectAnnotatedSampleId->errstr;
    while ( my @data = $selectAnnotatedSampleId->fetchrow_array ){
        $annotatedSampleId = $data[0];
    }
    return $annotatedSampleId;
}

# insert one rnaseq library individual sample and retrieve the value
# of internal autoincrement rnaSeqLibraryIndividualSampleId
sub insert_get_individual_sample {
    my ($insIndividualSample, $selectIndividualSampleId, $annotatedSampleId,
        $barcode, $sampleName, $debug) = @_;
    my $indiviudalSampleId;
    # insert individual sample
    if ($debug) {
        print 'INSERT INTO rnaSeqLibraryIndividualSample: ', $annotatedSampleId, ' - ', 
                    $barcode, ' - ', $sampleName, "\n";
    } else {
        $insIndividualSample->execute($annotatedSampleId, $barcode, $sampleName)
            or die $insIndividualSample->errstr;
    }
    # retrieve individual sample
    $selectIndividualSampleId->execute($annotatedSampleId, $barcode, $sampleName)
        or die $selectIndividualSampleId->errstr;
    while ( my @data = $selectIndividualSampleId->fetchrow_array ){
        $indiviudalSampleId = $data[0];
    }
    return $indiviudalSampleId;
}

sub read_sparse_matrix {
    my ($directory, $name) = @_;
    # read gene index file
    my $geneIndexFile = "$directory/$name.genes.txt";
    my %indexToGene = read_bustools_index_matrix($geneIndexFile);

    # read barcode index file
    my $barcodeIndexFile = "$directory/$name.barcodes.txt";
    my %indexToBarcode = read_bustools_index_matrix($barcodeIndexFile);
    # read coordinates file with counts
    my $coordinatesCountFile = "$directory/$name.mtx";
    my %coordinates_count = read_bustools_coordinates_matrix($coordinatesCountFile,
        \%indexToGene, \%indexToBarcode);
    return %coordinates_count;

}

sub read_bustools_index_matrix {
    my ($indexFile) = @_;
    my %indexHash;
    my $index = 1;
    for my $line ( read_file("$indexFile", chomp=>1) ){
        $indexHash{$index} = $line;
        $index++;
    }
    return %indexHash;
}

sub read_bustools_coordinates_matrix {
    my ($coordinatesFile, $indexToGeneRef, $indexToBarcodeRef) = (@_);
    my %indexToGene = %$indexToGeneRef;
    my %indexToBarcode = %$indexToBarcodeRef;
    my $foundHeader = 0;
    my %sparseMatrix;
    for my $line ( read_file("$coordinatesFile", chomp=>1) ){
        next  if ( $line =~ /^%/);
        # do not parse the line containing dimension of the matrix
        if(!$foundHeader) {
            $foundHeader = 1;
            next;
        }
        my @tmp = split(/\s+/, $line);
        if (scalar @tmp ne 3) {
            die "Lines of bustools coordinates file should contain 3 columns. @tmp";
        }
        #$sparseMatrix{barcode}->{indexToGene{gene} = ... (UMI count or cpm depending on the sparse matrix)
        # 
        if (!exists $indexToBarcode{$tmp[0]}) {
            warn "barcode : $tmp[1] - $tmp[0] - $tmp[2]\n";
        }
        # if gene does not exists it means it was an intergenic region removed in pre filtering of cpm
        elsif (!exists $indexToGene{$tmp[1]}) {
            print "gene :$tmp[1]-$tmp[0]-$tmp[2]\n";
        }

        $sparseMatrix{$indexToBarcode{$tmp[0]}}->{$indexToGene{$tmp[1]}} = $tmp[2];
    }
    return %sparseMatrix;
}

sub insert_individual_sample_gene_result {
    my ($insIndividualSampleGeneResult, $sparseMatrixCountRef, $sparseMatrixCpmRef,
        $genesResultsRef, $individualSampleId, $barcode, $debug) = @_;
    my %sparseMatrixCount = %$sparseMatrixCountRef;
    my %sparseMatrixCpm = %$sparseMatrixCpmRef;
    my %genesResults = %$genesResultsRef;

    for my $geneId (sort keys %{$sparseMatrixCount{$barcode}} ) {
        # check that the gene is present in the database. It is both a
        # security check and a way to remove intergenic regions
        next if (! exists $genesResults{$geneId} ||
            ! exists $sparseMatrixCount{$barcode}{$geneId});
        my $bgeeGeneId = $genesResults{$geneId};
        if (! exists $sparseMatrixCpm{$barcode}{$geneId}) {
            warn "Warning, gene $geneId has count for barcode $barcode but",
                " no abundance was generated";
            next;
        }
        if ($debug) {
            print 'INSERT INTO rnaSeqLibraryIndividualSampleGeneResult: ',
            $individualSampleId, ' - ', $bgeeGeneId, ' - ', "cpm", ' - ',
            $sparseMatrixCpm{$barcode}{$geneId}, ' - ', 0, ' - ',
            $sparseMatrixCount{$barcode}{$geneId}, ' - ',
            "high quality", ' - ', 'not excluded', "\n";
        } else {
            $insIndividualSampleGeneResult->execute($individualSampleId, $bgeeGeneId,
                'cpm', $sparseMatrixCpm{$barcode}{$geneId}, 0,
                $sparseMatrixCount{$barcode}{$geneId}, 'high quality', 'not excluded')
                    or die $insIndividualSampleGeneResult->errstr;
        }
    }
}

1;