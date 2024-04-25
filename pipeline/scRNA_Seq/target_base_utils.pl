#!/usr/bin/env perl

# Perl core modules
use strict;
use warnings;
use diagnostics;

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
    # $targetBaseLibrary{experimentId}->{libraryId (SRX...)}->{'authorAnatEntityAnnotation'} = ...
    # $targetBaseLibrary{experimentId}->{libraryId (SRX...)}->{'authorStageAnnotation'} = ...
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
        'Illumina HiSeq X Ten', 'Illumina NovaSeq 6000', 'NextSeq 500', 'NextSeq 550');

    my @valid_protocols = ('10X Genomics V2', '10X Genomics V3');

    my @validCellCompartments = ("scRNA-seq", "Sn-scRNA-seq");
    
    my %targetBaseLibrary;
    for my $line ( read_file("$targetBaseLibraryFile", chomp=>1) ){
        next  if ( $line =~ /^#libraryId/ or $line =~ /^\"#libraryId/ or $line =~ /^libraryId/);
        # there is currently 32 columns in the target base library annotation file
        my @tmp = map { bgeeTrim($_) } map { s/^\"//; s/\"$//; $_ } split(/\t/, $line, -1);

        my $libraryId                       = $tmp[0];
        my $experimentId                    = $tmp[1];
        my $platform                        = $tmp[2];
        my $anatEntityId                    = $tmp[4];
        my $cellTypeId                      = $tmp[6];
        my $stageId                         = $tmp[8];
        my $authorAnatEntityAnnotation      = $tmp[11];
        my $authorStageAnnotation           = $tmp[16];
        my $sex                             = $tmp[20];
        my $strain                          = $tmp[21];
        my $genotype                        = $tmp[22];
        my $speciesId                       = $tmp[23];
        my $cellCompartment                 = $tmp[24];
        my $protocol                        = $tmp[25];
        my $sequencedTranscriptPart         = $tmp[26];
        my $sampleName                      = $tmp[27];

        die "tsv field number problem [$line]\n"  if ( scalar @tmp != 33 );

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
                $sex = "NA";
            #TODO create hash with all potential sexes
            } else {
                if ($sex eq "F") {
                    $sex = "female";
                } elsif ($sex eq "M") {
                    $sex = "male";
                } elsif ( !grep(/^$sex$/, ("NA", "not annotated", "hermaphrodite", "mixed" ))) {
                    warn "Warning, wrong format for sex [$sex]\n";
                    $discarded = 1;
                }
            }
            if ($strain eq '' ){
                $strain= "NA";
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
            if ($sequencedTranscriptPart ne '3\'end' && $sequencedTranscriptPart ne '5\'end'){
                warn "Warning, wrong format for sequencedTranscriptPart",
                    "[$sequencedTranscriptPart]\n";
                $discarded = 1;
            } else {
                $sequencedTranscriptPart = '3prime';
            }

            if ( $discarded ){
                warn ' for experiment: ', $experimentId, ' - library: ', $libraryId,
                    ", library discarded!\n";
            } else {
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'platform'} = $platform;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'anatEntityId'} = $anatEntityId;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'cellTypeId'} = $cellTypeId;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'stageId'} = $stageId;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'authorAnatEntityAnnotation'} =
                    $authorAnatEntityAnnotation;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'authorStageAnnotation'} =
                    $authorStageAnnotation;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'sex'} = $sex;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'strain'} = $strain;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'genotype'} = $genotype;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'speciesId'} = $speciesId;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'cellCompartment'} = $cellCompartment;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'protocol'} = $protocol;
                $targetBaseLibrary{$experimentId}->{$libraryId}->{'sequencedTranscriptPart'} =
                    $sequencedTranscriptPart;
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
        # is it useful???
        my $status                          = $tmp[4];
        my $protocol                        = $tmp[7];
        my $protocolType                    = $tmp[8];
        my $comment                         = $tmp[13];


        die "tsv field number problem [$line]\n"  if ( scalar @tmp != 14 );

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
            if ($source eq 'HCA' ){
                warn "Warning, we remap HCA sources to SRA to be able to reuse the SRA experimentId\n";
                $source = 'SRA';
            }
            if ($source eq 'EBI' ){
                warn "Warning, we remap EBI source to ENA as it is the name of the source in the database\n";
                $source = 'ENA';
            }
            if ($protocol eq '' ){
                warn "Warning, wrong format for protocol [$protocol]\n";
                $discarded = 1;
            }
            if ($protocolType eq ''){
                warn "Warning, wrong format for experiment type [$protocolType]\n";
                $discarded = 1;
            }
            
            ## in order not to have "uninitialized value... warning, remove potentially
            # undef value per empty string
            $name = "" unless defined $name;
            $description = "" unless defined $description;
            $comment = "" unless defined $comment;
            
            if ( $discarded ){
                warn ' experiment: ', $experimentId, " discarded!\n";
            } else {
                if (grep(/^$protocolType$/, @experimentTypes)) {
                    $experiments{$experimentId}->{'name'} = $name;
                    $experiments{$experimentId}->{'description'} = $description;
                    $experiments{$experimentId}->{'source'} = $source;
                    $experiments{$experimentId}->{'status'} = $status;
                    $experiments{$experimentId}->{'comment'} = $comment;
                    $experiments{$experimentId}->{'protocol'} = $protocol;
                    $experiments{$experimentId}->{'experimentType'} = $protocolType;
                }
            }
        } else {
            warn 'Warning: experiment present several times in the experiment file: experiment: ',
            $experimentId, "\n";
        }
    }
    return %experiments;
}

##TODO: the mappedUMI info is not yet implemented in the pipeline but once it is
## we just need to uncomment corresponding lines in the function
sub getCallsSummaryAtLibraryAnnotatedLevel {
    my ($pipelineCallsSummary) = @_;
    # $callsSummary{libraryId}{cellTypeId}->{authorCellTypeAnnotation}->{'abundanceThreshold'}
    # $callsSummary{libraryId}{cellTypeId}->{authorCellTypeAnnotation}->{'allGenesPercentPresent'}
    # $callsSummary{libraryId}{cellTypeId}->{authorCellTypeAnnotation}->{'proteinCodingGenesPercentPresent'}
    # $callsSummary{libraryId}{cellTypeId}->{authorCellTypeAnnotation}->{'intergenicRegionsPercentPresent'}
    # $callsSummary{libraryId}{cellTypeId}->{authorCellTypeAnnotation}->{'pValueThreshold'}
    # $callsSummary{libraryId}{cellTypeId}->{authorCellTypeAnnotation}->{'meanRefIntergenic'}
    # $callsSummary{libraryId}{cellTypeId}->{authorCellTypeAnnotation}->{'sdRefIntergenic'}

    my %callsSummary;
    for my $line ( read_file("$pipelineCallsSummary", chomp=>1) ){
        next  if ( $line =~ /^libraryId/ );
        # there is currently 18 columns in the target base library annotation file
        my @tmp = map { bgeeTrim($_) } map { s/^\"//; s/\"$//; $_ } split(/\t/, $line);

        my $libraryId                        = $tmp[0];
        my $cellTypeId                       = $tmp[1];
        my $authorCellTypeAnnotation         = $tmp[2];
        my $abundanceThreshold               = $tmp[3];
        my $allGenesPercentPresent           = $tmp[6];
        my $proteinCodingGenesPercentPresent = $tmp[9];
        my $intergenicRegionsPercentPresent  = $tmp[12];
        my $pValueThreshold                  = $tmp[13];
        my $meanRefIntergenic                = $tmp[14];
        my $sdRefIntergenic                  = $tmp[15];


        die "tsv field number problem [$line]\n"  if ( scalar @tmp != 18 );

        if ( !defined $callsSummary{$libraryId}{$cellTypeId}{$authorCellTypeAnnotation} ){
            # Perform format checks
            my $discarded = 0;
            if ($libraryId eq '' ){
                warn "Warning, wrong format for libraryId [$libraryId]\n";
                $discarded = 1;
            }
            if ($cellTypeId eq '' ){
                warn "Warning, wrong format for cellTypeId [$cellTypeId]\n";
                $discarded = 1;
            }
            if ($authorCellTypeAnnotation eq '' ){
                warn "authorCellTypeAnnotation is empty for library $libraryId and celltype $cellTypeId. The library/celltype is not discarded.\n";
            }
            if ($abundanceThreshold eq '' ){
                warn "Warning, wrong format for abundanceThreshold [$abundanceThreshold]\n";
                $discarded = 1;
            }
            if ($allGenesPercentPresent eq ''){
                warn "Warning, wrong format for allGenesPercentPresent type [$allGenesPercentPresent]\n";
                $discarded = 1;
            }
            if ($proteinCodingGenesPercentPresent eq ''){
                warn "Warning, wrong format for proteinCodingGenesPercentPresent type [$proteinCodingGenesPercentPresent]\n";
                $discarded = 1;
            }
            if ($intergenicRegionsPercentPresent eq ''){
                warn "Warning, wrong format for intergenicRegionsPercentPresent type [$intergenicRegionsPercentPresent]\n";
                $discarded = 1;
            }
            if ($pValueThreshold eq ''){
                warn "Warning, wrong format for pValueThreshold type [$pValueThreshold]\n";
                $discarded = 1;
            }
            if ($meanRefIntergenic eq ''){
                warn "Warning, wrong format for meanRefIntergenic type [$meanRefIntergenic]\n";
                $discarded = 1;
            }
            if ($sdRefIntergenic eq ''){
                warn "Warning, wrong format for sdRefIntergenic type [$sdRefIntergenic]\n";
                $discarded = 1;
            }
            # if ($mappedUMIs eq ''){
            #     warn "Warning, wrong format for mappedUMIs type [$mappedUMIs]\n";
            #     $discarded = 1;
            # }
            if ( $discarded ){
                warn ' libraryId: ', $libraryId, ", cellTypeId: ", $cellTypeId, " discarded!\n";
            } else {
                $callsSummary{$libraryId}->{$cellTypeId}->{$authorCellTypeAnnotation}->{'abundanceThreshold'} = $abundanceThreshold;
                $callsSummary{$libraryId}->{$cellTypeId}->{$authorCellTypeAnnotation}->{'allGenesPercentPresent'} = $allGenesPercentPresent;
                $callsSummary{$libraryId}->{$cellTypeId}->{$authorCellTypeAnnotation}->{'proteinCodingGenesPercentPresent'} = $proteinCodingGenesPercentPresent;
                $callsSummary{$libraryId}->{$cellTypeId}->{$authorCellTypeAnnotation}->{'intergenicRegionsPercentPresent'} = $intergenicRegionsPercentPresent;
                $callsSummary{$libraryId}->{$cellTypeId}->{$authorCellTypeAnnotation}->{'pValueThreshold'} = $pValueThreshold;
                $callsSummary{$libraryId}->{$cellTypeId}->{$authorCellTypeAnnotation}->{'meanRefIntergenic'} = $meanRefIntergenic;
                $callsSummary{$libraryId}->{$cellTypeId}->{$authorCellTypeAnnotation}->{'sdRefIntergenic'} = $sdRefIntergenic;
            }
        } else {
            warn 'Warning: couple libraryId/cellTypeId present several times in the file: libraryId: ',
            $libraryId, ", cellTypeId : ", $cellTypeId, "\n";
        }
    }
    return %callsSummary;
}

# Extract library information from the annotation generated by the bgee pipeline for pipeline full length and target based
# added an argument to check if the metadata file comes from FL or TB as the target based file has one more column
sub get_processed_libraries_info {
    my ($targetBaseLibraryFile, $isTargetBased) = @_;
    # $experiments{experimentId}->{libraryId}->{runId}->{'readCount'} = ...
    # $experiments{experimentId}->{libraryId}->{'libraryLayout'} = ...
    
    my @validLibraryType = ('paired', 'single');

    my %libInfos = ();
    for my $line ( read_file("$targetBaseLibraryFile", chomp=>1) ){
        next  if ( $line =~ /^sample_accession/ or $line =~ /^\"sample_accession/ );
        my @tmp = map { bgeeTrim($_) } map { s/^\"//; s/\"$//; $_ } split(/\t/, $line, -1);
        if ($isTargetBased) {
            # there is currently 12 columns in the metadata_info_10X.txt from target based pipeline
            die "tsv field number problem [$line]\n"  if ( scalar @tmp != 12 );
        } else {
            # there is currently 11 columns in the metadata_info.txt file from full length pipeline
            die "tsv field number problem [$line]\n"  if ( scalar @tmp != 11 );
        }
        my $runId                           = $tmp[3];
        my $experimentId                    = $tmp[1];
        my $libraryId                       = $tmp[2];
        my $speciesId                       = $tmp[5];
        my $speciesName                     = $tmp[6];
        my $libraryType                     = $tmp[8];
        my $fastqFTP                        = $tmp[9];
        my $submittedFTP                    = $tmp[10];
        my $downloadSource = "";
        if ($isTargetBased) {
            $downloadSource                 = $tmp[11];
        }

        if (!defined $libInfos{$experimentId}->{$libraryId}->{'runIds'}->{$runId}) {
            # Perform format checks
            # do not check format for genotype and whiteList as they can be empty
            my $discarded = 0;
            if ($runId eq '' ){
                warn "Warning, wrong format for runId [$runId]\n";
                $discarded = 1;
            }
            if ($libraryId eq '' ){
                warn "Warning, wrong format for libraryId [$libraryId]\n";
                $discarded = 1;
            }
            if ($experimentId eq '' ){
                warn "Warning, wrong format for experimentId [$experimentId]\n";
                $discarded = 1;
            }
            if ($libraryType eq 'PAIRED') {
                $libraryType = 'paired';
            } elsif ($libraryType eq 'SINGLE') {
                $libraryType = 'single';
            } elsif (all { $libraryType !~ /^$_/ } @validLibraryType) {
                warn "Warning, wrong format for libraryType [$libraryType]\n";
                $discarded = 1;
            }
            if ($speciesId eq '' ){
                warn "Warning, wrong format for speciesId [$speciesId]\n";
                $discarded = 1;
            }
            if ($speciesName eq '' ){
                warn "Warning, wrong format for speciesId [$speciesName]\n";
                $discarded = 1;
            }
            if ($isTargetBased && $downloadSource eq '' ){
                warn "Warning, wrong format for downloadSource [$downloadSource]\n";
                $discarded = 1;
            }
            if ( $discarded ){
                warn ' experimentId: ', $experimentId, ", libraryId: ", $libraryId, " discarded!\n";
            } else {
                $libInfos{$experimentId}->{$libraryId}->{'runIds'}->{$runId}->{'libraryType'} =
                    $libraryType;
                if ($isTargetBased) {
                    $libInfos{$experimentId}->{$libraryId}->{'runIds'}->{$runId}->{'downloadSource'} =
                        $downloadSource;
                }
                $libInfos{$experimentId}->{$libraryId}->{'runIds'}->{$runId}->{'submittedFTP'} =
                    $submittedFTP;
                $libInfos{$experimentId}->{$libraryId}->{'runIds'}->{$runId}->{'fastqFTP'} =
                    $fastqFTP;
                $libInfos{$experimentId}->{$libraryId}->{'speciesId'} = $speciesId;
                $libInfos{$experimentId}->{$libraryId}->{'speciesName'} = $speciesName;
            }
        } else {
            warn 'Warning: run present several times in the metadata file: experiment: ',
            $experimentId, ' - library: ', $libraryId, 'run: ', $runId, "\n";
        }
    }

    return %libInfos;
}

sub getCallsInfoPerLibrary {
    my ($callsInfoFile) = @_;
    # $callsInfo{$geneId}->{'cpm'}
    # $callsInfo{$geneId}->{'sumUMI'}
    # $callsInfo{$geneId}->{'zScore'}
    # $callsInfo{$geneId}->{'pValue'}

    my %callsInfo = ();
    for my $line ( read_file("$callsInfoFile", chomp=>1) ){
        next  if ( $line =~ /^gene_id/);
        # there is currently 10 columns in the Calls_LIBRARY_ID.tsv file
        my @tmp = map { bgeeTrim($_) } map { s/^\"//; s/\"$//; $_ } split(/\t/, $line, -1);

        my $geneId                    = $tmp[0];
        my $sumUMI                    = $tmp[1];
        my $cpm                       = $tmp[2];
        my $zScore                    = $tmp[5];
        my $pValue                    = $tmp[6];
        my $cellTypeId                = $tmp[8];
        my $authorCellTypeAnnotation  = $tmp[9];

        die "tsv field number problem [$line]\n"  if ( scalar @tmp != 10 );

        # Perform format checks
        my $discarded = 0;
        if ($geneId eq '' ){
            warn "Warning, wrong format for geneId [$geneId]\n";
            $discarded = 1;
        }
        if ($sumUMI eq '' ){
            warn "Warning, wrong format for sumUMI [$sumUMI]\n";
            $discarded = 1;
        }
        if ($cpm eq '' ){
            warn "Warning, wrong format for cpm [$cpm]\n";
            $discarded = 1;
        }
        if ($zScore eq '' ){
            warn "Warning, wrong format for zScore [$zScore]\n";
            $discarded = 1;
        # non numerical values of zscore are stored as NULL in the database
        # dbi iquivalent to null is undef so if value of zscore is NA we modified
        # it to be inserted as NULL in the database
        } elsif ($zScore eq "NA") {
            $zScore = undef;
        }
        if ($pValue eq '' ){
            warn "Warning, wrong format for pValue [$pValue]\n";
            $discarded = 1;
        } elsif ($pValue eq "NA") {
            $pValue = 1;
        }
        if ($discarded) {
            warn "gene: $geneId, cellTypeId $cellTypeId, cell-type author annotation $authorCellTypeAnnotation discarded!\n";
        } else {
            $callsInfo{$geneId}->{'cpm'} = $cpm;
            $callsInfo{$geneId}->{'sumUMI'} = $sumUMI;
            $callsInfo{$geneId}->{'zScore'} = $zScore;
            $callsInfo{$geneId}->{'pValue'} = $pValue;
        }
    }

    return %callsInfo;
}

# Extract the mapping from barcodes to cell types. Also extract cell type Ids
sub getBarcodeToCellType {
    my ($barcodeToCellTypeFile) = @_;
    # $barcodes{libraryId}->{'barcodes'}->{barcode}->{'clusterId'} = ...
    # $barcodes{libraryId}->{'barcodes'}->{barcode}->{'authorCellTypeAnnotation'} = ...
    # $barcodes{libraryId}->{'barcodes'}->{barcode}->{'annotationStatus'} = ...
    # $barcodes{libraryId}->{'celltypes'}->{clusterId}->{'cellTypeId'} = ()
    # $barcodes{libraryId}->{'celltypes'}->{clusterId}->{'authorCellTypeAnnotation'} = ()

    my %barcodes = ();
    for my $line ( read_file("$barcodeToCellTypeFile", chomp=>1) ){
        next  if ( $line =~ /^barcode/ or $line =~ /^\"barcode/ );
        # there is currently 10 columns in the metadata_info_10X.txt file
        my @tmp = map { bgeeTrim($_) } map { s/^\"//; s/\"$//; $_ } split(/\t/, $line);

        my $barcode                         = $tmp[0];
        my $libraryId                       = $tmp[2];
        my $authorCellTypeAnnotation        = $tmp[5];
        my $cellTypeId                      = $tmp[9];
        my $cellTypeAnnotationStatus        = $tmp[11];
        my $bgeeClusterId                   = $tmp[14];

        die "tsv field number problem [$line]\n"  if (scalar @tmp != 15);

        if (!defined $barcodes{$libraryId}->{'barcodes'}->{$barcode}) {
            # Perform format checks
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
            } else {
                #it is possible to have more than one id in this column.
                #For now, as a temporary solution, if there are more than one term
                #and exactly 2 terms, then we take the second term (this is because,
                #after discussing with Anne, it looks like the first term is a more specific
                #anat. entity term
                #TODO: modify this part of the pipeline for Bgee 16 and use the ontology to find
                #the closest common ancestor
                my @cellTypeIds = split( /,\s+/, $cellTypeId);
                if(scalar @cellTypeIds > 2) {
                    die "a maximum of 2 cellTypeIds can be present in the column cellTypeId";
                } elsif (scalar @cellTypeIds == 2) {
                    $cellTypeId = $cellTypeIds[1];
                }
               
            }
            if ($cellTypeAnnotationStatus eq '' ){
                warn "Warning, wrong format for cellTypeAnnotationStatus [$cellTypeAnnotationStatus]\n";
                $discarded = 1;
            }
            if ($bgeeClusterId eq '' ){
                warn "Warning, wrong format for bgeeClusterId [$bgeeClusterId]\n";
                $discarded = 1;
            }
            if ($authorCellTypeAnnotation eq '' ){
                warn "Warning, authorCellTypeAnnotation is empty. library : $libraryId celltype : $cellTypeId. The cluster has not been discarded\n";
            }
            if ($discarded) {
                warn 'libraryId: ', $libraryId, ', barcode ', $barcode, ', bgeeClusterId ', $bgeeClusterId,
                ', authorCellTypeAnnotation ', $authorCellTypeAnnotation,
                ', cellTypeId ', $cellTypeId, " discarded!\n";
            } else {
                $barcodes{$libraryId}->{'barcodes'}->{$barcode}->{'clusterId'} =
                    $bgeeClusterId;
                $barcodes{$libraryId}->{'barcodes'}->{$barcode}->{'authorCellTypeAnnotation'} =
                    $authorCellTypeAnnotation;
                $barcodes{$libraryId}->{'barcodes'}->{$barcode}->{'cellTypeId'} =
                    $cellTypeId;
                $barcodes{$libraryId}->{'barcodes'}->{$barcode}->{'annotationStatus'} =
                    $cellTypeAnnotationStatus;
                $barcodes{$libraryId}->{'clusters'}->{$bgeeClusterId}->{'cellTypeId'} = $cellTypeId;
                $barcodes{$libraryId}->{'clusters'}->{$bgeeClusterId}->{'authorCellTypeAnnotation'} =
                    $authorCellTypeAnnotation;
                $barcodes{$libraryId}->{'cellTypeIds'}->{$cellTypeId} = ();
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
    my ($libraryId, $conditionId, $cellTypeAuthorAnnotation, $anatEntityAuthorAnnotation,
        $stageAuthorAnnotation, $abundanceUnit, $meanRefIntergenic, $sdRefIntergenic,
        $abundanceThreshold, $allGenesPercentPresent, $proteinCodingGenesPercentPresent,
        $intergenicRegionsPercentPresent, $pValueThreshold, $allUMIsCount, $mappedUMIsCount,
        $isDropletBased, $barcode, $time, $timeUnit, $freeTextAnnotation, $insAnnotatedSample,
        $selectAnnotatedSampleId, $debug) = @_;

    #insert annotated sample
    if ($debug) {
        my $timeNullable = $time;
        if (!defined $time) {
            $timeNullable = 'undef';
        }
        print 'INSERT INTO rnaSeqLibraryAnnotatedSample: ',
                    $libraryId,                        ' - ',
                    'conditionId',                      ' - ',
                    $cellTypeAuthorAnnotation,         ' - ',
                    $anatEntityAuthorAnnotation,       ' - ',
                    $stageAuthorAnnotation,            ' - ',
                    $abundanceUnit,                    ' - ',
                    $meanRefIntergenic,                ' - ',
                    $sdRefIntergenic,                  ' - ',
                    $abundanceThreshold,               ' - ',
                    $allGenesPercentPresent,           ' - ',
                    $proteinCodingGenesPercentPresent, ' - ',
                    $intergenicRegionsPercentPresent,  ' - ',
                    $pValueThreshold,                  ' - ',
                    $allUMIsCount,                     ' - ',
                    $mappedUMIsCount,                  ' - ',
                    $isDropletBased,                   ' - ',
                    $barcode,                          ' - ',
                    $timeNullable,                             ' - ',
                    $timeUnit,                         ' - ',
                    $freeTextAnnotation,
                    "\n";
    } else {
        $insAnnotatedSample->execute($libraryId, $conditionId, $cellTypeAuthorAnnotation,
        $anatEntityAuthorAnnotation, $stageAuthorAnnotation, $abundanceUnit, $meanRefIntergenic,
        $sdRefIntergenic, $abundanceThreshold, $allGenesPercentPresent,
        $proteinCodingGenesPercentPresent, $intergenicRegionsPercentPresent, $pValueThreshold,
        $allUMIsCount, $mappedUMIsCount, $isDropletBased, $barcode, $time, $timeUnit,
        $freeTextAnnotation)
            or die $insAnnotatedSample->errstr;
    }
    
    #retrieve annotated sample ID
    my $annotatedSampleId = ();
    $selectAnnotatedSampleId->execute($libraryId, $conditionId, $cellTypeAuthorAnnotation)
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
        print 'INSERT INTO rnaSeqLibraryIndividualSample: ', 'annotatedSampleId', ' - ', 
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

# function loading both matrices in the same hash to reduce memory usage
sub read_count_and_cpm_matrices {
    my ($countDirectory, $countName, $cpmDirectory, $cpmName) = @_;
    my %sparseMatrixCount = read_sparse_matrix($countDirectory, $countName);
    my %sparseMatrixCpm = read_sparse_matrix($cpmDirectory, $cpmName);
    my %sparseMatrixCombined;
    for my $barcode (keys %sparseMatrixCount) {
        for my $geneId (keys %{$sparseMatrixCount{$barcode}}) {
            $sparseMatrixCombined{$barcode}->{$geneId}->{'count'} =
                $sparseMatrixCount{$barcode}{$geneId};
            if (exists $sparseMatrixCpm{$barcode}{$geneId}) {
                $sparseMatrixCombined{$barcode}->{$geneId}->{'cpm'} =
                $sparseMatrixCpm{$barcode}{$geneId};
            }
        }
    }
    return %sparseMatrixCombined;
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

1;
