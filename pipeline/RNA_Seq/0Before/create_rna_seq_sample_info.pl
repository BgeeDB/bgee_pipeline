#!/usr/bin/env perl
# Julien Roux

# This scripts parses the annotation file, and retrieves information on the species (from Bgee database) and the SRA IDs, usign NCBI e-utils

# Dec 17, 2015:
# - merged with get_sra_id.pl to simplify pipeline and create less intermediate info files.
# - Added output file as argument.
# - The script now issues warnings if the organism and platform information in the annotation don't match the information on the SRA record (before, TRUE or FALSE flags columsn were present in output file)
# - The script checks annotation file RNASeqLibraryPlatformChecks.tsv, which documents previous investigations on the warnings from this script (previous pipeline), so that no new warning is issued if there is really a mistake on the SRA record

# Mar 14, 2016
# - Now should work on out annotation file merged with the Wormbase annotation
# - An external file is checked for libraries which were manually checked (with decision to include them or not).
# - warnings issues for ncRNA-seq
# - try to detect CAGE-seq ("Library selection" == CAGE flag in XML), SAGE-seq, DeepSage and issues warning
# - issue warning if read length too short (<36nt): miRNA or special library?

# Mar 29, 2016
# - Added if the species comes from ensembl or ensembl metazoa, in rna_seq_sample_info.txt
#
# TODO a lot of information could be extracted from SRA for the annotators, who will just have to check. We could make a new script that is run for each SRX ID and outputs the annotation line to copy paste in annotation spreadsheet
# We could get tissue source to compare to annotation. Problem: info not present in most SRA records
# if ( $info =~ /<TAG>OrganismPart<\/TAG><VALUE>([^<]+?)<\/VALUE>/ ) {
#   $source = $1;
# }
# We could get experiment (GSE ID) to compare to annotation. Problem: info not present for GTEx samples
# if ( $info =~ // ){
#   $series = $1;
# }
# Any other information we could get to help or to compare to annotation?

# Perl core modules
use strict;
use warnings;
use diagnostics;

use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../.."; # Get lib path for Utils.pm
use Utils;
use File::Slurp;
use List::MoreUtils qw(all any);
use LWP::Simple;

# Define arguments & their default value
my ($bgee_connector)         = ('');
my ($RNAseqLib)              = ('');
my ($RNAseqLibChecks)        = ('');
my ($RNAseqLibWormExclusion) = ('');
my ($extraMapping)           = ('');
my ($outFile)                = ('');
my ($debug)                  = (0);
my ($sample)                 = '';
my %opts = ('bgee=s'                   => \$bgee_connector,     # Bgee connector string
            'RNAseqLib=s'              => \$RNAseqLib,
            'RNAseqLibChecks=s'        => \$RNAseqLibChecks,
            'RNAseqLibWormExclusion=s' => \$RNAseqLibWormExclusion,
            'extraMapping=s'           => \$extraMapping,
            'outFile=s'                => \$outFile,
            'debug'                    => \$debug,
            'sample=s'                 => \$sample,
           );

# Check arguments
my $test_options = Getopt::Long::GetOptions(%opts);
if ( !$test_options || $bgee_connector eq '' || $RNAseqLib eq '' || $RNAseqLibChecks eq '' || $RNAseqLibWormExclusion eq '' || $outFile eq ''){
    print "\n\tInvalid or missing argument:
\te.g. $0  -bgee=\$(BGEECMD) -RNAseqLib=\$(RNASEQ_LIB_FILEPATH_FULL) -RNAseqLibChecks=\$(RNASEQ_LIB_CHECKS_FILEPATH_FULL) -RNAseqLibWormExclusion=\$(RNASEQ_LIB_EXCLUSION_FILEPATH_WORM) -outFile=\$(RNASEQ_SAMPINFO_FILEPATH) -extraMapping=\$(EXTRAMAPPING_FILEPATH)
\t-bgee                   Bgee connector string
\t-RNAseqLib              RNAseq Libraries from annotation file
\t-RNAseqLibChecks        RNAseq Libraries previously checked for platform annotation errors
\t-RNAseqLibWormExclusion RNAseq Libraries excluded from WormBase annotation
\t-outFile                Output file: TSV with all species and SRA information from NCBI
\t-extraMapping           Extra mapping file
\t-debug                  more verbose output
\t-sample                 Analyse this specific sample ONLY
\n";
    exit 1;
}

print 'Connecting to Bgee to retrieve species/stage/organ information... ';
# Bgee db connection
my $dbh = Utils::connect_bgee_db($bgee_connector);

# Get species list in bgee db
my %species;
my $selSpecies = $dbh->prepare('SELECT species.speciesId, CONCAT(species.genus, " ", species.species), species.genomeFilePath, dataSource.dataSourceName FROM species INNER JOIN dataSource ON species.dataSourceId = dataSource.dataSourceId');
$selSpecies->execute()  or die $selSpecies->errstr;
while ( my @data = $selSpecies->fetchrow_array ){
    $species{$data[0]}->{'organism'}       = $data[1];
    $species{$data[0]}->{'genomeFilePath'} = $data[2];
    $species{$data[0]}->{'database'}       = $data[3];
}
$selSpecies->finish;

# Retrieve all organs/stages from Bgee
my %organs = %{ Utils::getBgeedbOrgans($dbh) };
my %stages = %{ Utils::getBgeedbStages($dbh) };

# Parse extra mapping info for currently too up-to-date annotations
##UnmappedId    UnmappedName    UberonID    UberonName    Comment
my %extra = map  { my @tmp = split(/\t/, $_, -1); if ( $tmp[2] ne '' && $tmp[0] ne '' ){ $tmp[0] => $tmp[2] } else { 'nonono' => 'nonono' } }
            grep { !/^#/ }
            read_file("$extraMapping", chomp => 1);

$dbh->disconnect;
print "Done\n";


# Read annotation files
print "Reading annotation file $RNAseqLib... ";
die "[$RNAseqLib] file is empty\n"  if ( -z $RNAseqLib );
my %tsv = %{ Utils::read_spreadsheet("$RNAseqLib", "\t", 'csv', '"', 1) };
# See modified headers, issue #90 on gitHub:
# libraryId    experimentId    chipTypeId    organId    organName    uberonId    uberonName    stageId    stageName    infoOrgan    infoStage    sampleTitle    sampleSource    sampleDescription    sampleCharacteristics    organAnnotationStatus    organBiologicalStatus    stageAnnotationStatus    stageBiologicalStatus    sex    strain    speciesId    comment    annotatorId    lastModificationDate    replicate    infoReplicate    SRSId    tags
# SRX081869    GSE30352    Illumina Genome Analyzer IIx            UBERON:0000955    brain    GgalDv:0000008    1-year-old chicken stage    Brain    ~1 year, adult    gga br F 1                perfect match    not documented    other    partial sampling    F    Red Junglefowl    9031    Chicken    ANN    2013-08-23    1    GEO - [...] we generated RNA-Seq data [...] of brain (cerebral cortex or whole brain without cerebellum), cerebellum, heart, kidney, liver and testis (usually from one male and one female per somatic tissue and two males for testis)[...]
print "Done\n";

print "Reading annotation file $RNAseqLibChecks... ";
warn "Warning: [$RNAseqLibChecks] file is empty.\n"  if ( -z $RNAseqLibChecks );
my %tsv_checks = %{ Utils::read_spreadsheet("$RNAseqLibChecks", "\t", 'csv', '"', 1) };
#libraryId      platformOld     platformNew     platformSRA     platformGEO     platformPaper   comment annotatorId     lastModificationDate
my %checked_libraries;
for my $library (@{$tsv_checks{'libraryId'}}){
    # print "\t$library platform was previously checked, and problems should be solved.\n";
    $checked_libraries{$library} = ();
}
print "Done\n";

print "Reading annotation file $RNAseqLibWormExclusion... ";
warn "Warning: [$RNAseqLibWormExclusion] file is empty.\n"  if ( -z $RNAseqLibWormExclusion );
%tsv_checks = %{ Utils::read_spreadsheet("$RNAseqLibWormExclusion", "\t", 'csv', '"', 1) };
#libraryId    excluded    comment    annotatorId    lastModificationDate
my %checked_libraries_worm;
for my $i ( 0..$#{$tsv_checks{'libraryId'}} ) {
    $checked_libraries_worm{$tsv_checks{'libraryId'}[$i]} = $tsv_checks{'excluded'}[$i];
}
print "Done\n";

print "For each library, retrieving SRA record information and writing in output file...\n";
open (my $OUT, '>', "$outFile")  or die "Cannot write [$outFile]\n";
# output file header
print {$OUT} join("\t", '#libraryId', 'experimentId', 'speciesId', 'organism', 'genomeFilePath', 'database', 'platform', 'libraryType', 'libraryInfo', 'readLength' ,'runIds'), "\n";

# Retrieve SRA information and IDs for each SRX ID, and write down output file. This was originally in script get_sra_id.pl, but it is simplified here (only 1 query, see email Sebastien 15/12/15). Example of a query URL: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?email=bgee@sib.swiss&api_key=a2546089861bb524068974de020c591a1307&save=efetch&db=sra&rettype=xml&term=SRX1152842
SAMPLE:
for my $i ( 0..$#{$tsv{'libraryId'}} ) {
    my $strategy;
    my $libraryType  = '';
    my $libraryInfo  = '';
    my $readLength   = '';
    my $libraryId    = $tsv{'libraryId'}[$i];
    my $experimentId = $tsv{'experimentId'}[$i];
    my $tag          = $tsv{'tags'}[$i] // '';
    my $anatID       = $tsv{'uberonId'}[$i];
    my $stageID      = $tsv{'stageId'}[$i];

    #NOTE Restrict analysis to this specific sample id (library or experiment)
    if ( $sample ne '' ){
        next SAMPLE  if ( $libraryId ne $sample && $experimentId ne $sample );
    }

    # line commented: skipped
    next SAMPLE  if ( $libraryId =~ /^#/ );
    next SAMPLE  if ( $tag =~ /^ScRNA-seq$/i ); #NOTE Discard Single-cell sequencing, not this pipeline part
    print "\t$libraryId\t$experimentId\n";

    # skipped if no species matching the annotated speciesId in the database
    if ( ! exists($species{ $tsv{'speciesId'}[$i] }) ){
        warn "\tWarning: No species/genome assembly defined in the database for [$libraryId] [$experimentId]: [$tsv{'speciesId'}[$i]]. This library was not printed in output file.\n";
        next SAMPLE;
    }
    # skipped if anatID or stageID not in bgee db (e.g. too recent Uberon used by annotators)
    if ( (! exists($organs{ $anatID })  && ! exists($extra{ $anatID }))  || !$anatID ){
        warn "\tWarning: [$anatID] anatId  does not exist in the database for [$libraryId] [$experimentId]. This library was not printed in output file.\n";
        next SAMPLE;
    }
    if ( (! exists($stages{ $stageID }) && ! exists($extra{ $stageID })) || !$stageID ){
        warn "\tWarning: [$stageID] stageId does not exist in the database for [$libraryId] [$experimentId]. This library was not printed in output file.\n";
        next SAMPLE;
    }

    if ( $libraryId !~ /^[SEDC]RX\d+$/ ){
        warn "\tProblem: [$libraryId] [$experimentId] not in SRX format. This library was not printed in output file.\n";
        next SAMPLE;
    }

    # perform query to NCBI
    my $info = get("https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?email=bgee\@sib.swiss&api_key=a2546089861bb524068974de020c591a1307&save=efetch&db=sra&rettype=xml&term=$libraryId");
    my ($organism, $platform, $source, $series) = ('', '', '', '');
    my @SRR;

    # if query answer is empty / undef: connection unsuccessful or wrong ID?
    if ( !$info ){
        warn "\tProblem: No SRA information for [$libraryId] [$experimentId]: connection unsuccessful? Incorrect Id? This library was not printed in output file.\n";
        next SAMPLE;
    }
    if ( $info =~ /<ERROR.+SRA Experiment $libraryId is not public/ ){
        warn "\tWarning: [$libraryId] [$experimentId] is not public\n";
        next SAMPLE;
    }

    # if query was successful, an XML file is returned:
    else {
        # If it is a Wormbase library that was previously manually checked for inclusion, do not make additional checks
        if ( exists $checked_libraries_worm{$libraryId} ){
            #NOTE See https://gitlab.sib.swiss/Bgee/expression-annotations/issues/33
            # if the library was excluded:
            #SRP000401(SRX035162, SRX036882, SRX036967, SRX036969, SRX036970, SRX047787), SRP001010, SRP015688
            if ( $checked_libraries_worm{$libraryId} eq 'TRUE' ){
                warn "\tInfo: [$libraryId] [$experimentId] from Wormbase was excluded following manual curation (see file $RNAseqLibWormExclusion). This library was not printed in output file.\n";
                next SAMPLE
            }
            else {
                # Validated in https://gitlab.sib.swiss/Bgee/expression-annotations/issues/33
                #GSE16552, GSE22410, SRP000401(SRX037197, SRX050630), SRP015688
                warn "\tInfo: [$libraryId] [$experimentId] from Wormbase was manually curated (see file $RNAseqLibWormExclusion), and it was decided to include it. Next warning messages can probably be ignored.\n";
            }
        }
        # Perform a series of sanity checks
        # Check it is a RNA-seq library
        my @valid_lib_strategy = ('GSE30606', 'GSE30617', 'SRP000401', 'E-MTAB-2449', 'SRP003905', 'GSE16552', 'DRP000571', 'SRP000304', 'SRP004363', 'SRP005402', 'SRP033585');
        #NOTE See https://gitlab.sib.swiss/Bgee/expression-annotations/issues/31
        #     See https://gitlab.sib.swiss/Bgee/expression-annotations/issues/82
        # GSE30617      OK for Anne
        # SRP000401     OK for Anne
        # E-MTAB-2449   'validated with lower confidence'
        # SRP003905     '"EST" but that it is Illumina paired-end, so I think it is RNA-seq'
        # GSE16552      'Contradictory info but for now kept'. Also in https://gitlab.sib.swiss/Bgee/expression-annotations/issues/33
        # DRP000571     OK for Anne
        # SRP000304     OK, FL-cDNA is RNA-Seq
        my @invalid_lib_strategies = ('ncRNA-Seq', 'ATAC-seq', 'MAINE-Seq', 'MNase-Seq', 'FAIRE-seq', 'DNase-Hypersensitivity', 'DNase-seq');
        $info =~ /<LIBRARY_STRATEGY>([^<]+)<\/LIBRARY_STRATEGY>/; # [^<] prevents matching to '<' character
        $strategy = $1;
        if ( any { lc($strategy) eq lc($_) } @invalid_lib_strategies ){
            warn "\tProblem: [$libraryId] [$experimentId] is [$strategy], which is not supported by our pipeline for now. Please comment out. This library was not printed in output file.\n";
            next SAMPLE;
        }
        elsif ( $strategy !~ /^RNA-Seq$/i and (all { $experimentId ne $_ } @valid_lib_strategy) ){
            warn "\tProblem: [$libraryId][$experimentId] does not seem to be RNA-seq but [$strategy]. Please check. This library was printed in output file.\n";
        }

        # Particular types of RNA-seq. See https://www.ncbi.nlm.nih.gov/books/NBK49283/ or https://www.ebi.ac.uk/ena/submit/preparing-xmls
        my @valid_selection_methods = ('cDNA', 'Inverse rRNA', 'oligo-dT', 'Oligo-dT', 'PCR', 'PolyA', 'RANDOM', 'RANDOM PCR', 'RT-PCR');
        #NOTE See https://gitlab.sib.swiss/Bgee/expression-annotations/issues/30
        #     See https://gitlab.sib.swiss/Bgee/expression-annotations/issues/82
        my @valid_lib_selection     = ('DRP003809', 'E-MTAB-5895', 'SRP012049', 'SRP018740', 'SRP021223', 'SRP051959', 'SRP055497', 'SRP058036', 'SRP081080', 'SRP082291', 'SRP082342', 'SRP082454', 'SRP092904', 'SRP106023', 'SRP116580', 'SRP045680', 'SRP018725', 'SRP053164', 'SRP056073', 'SRP072263', 'SRP036185', 'SRP058798');
        $info =~ /<LIBRARY_SELECTION>([^<]+)<\/LIBRARY_SELECTION>/; # [^<] prevents matching to '<' character
        my $selection = $1;
        if ( $selection =~ /CAGE/ ){
            warn "\tProblem: [$libraryId][$experimentId] seems to be CAGE-seq, which is not supported by our pipeline for now. Please check. This library was not printed in output file.\n";
            next SAMPLE;
        }
        elsif ( $selection =~ /RACE/ ){
            warn "\tProblem: [$libraryId][$experimentId] seems to be RACE-seq, which is not supported by our pipeline for now. Please check. This library was not printed in output file.\n";
            next SAMPLE;
        }
        elsif ( $selection =~ /size fractionation/ ){
            warn "\tProblem: [$libraryId][$experimentId] seems to be size-fractionated, which biases the expression estimates. Please comment out. This library was not printed in output file.\n";
            next SAMPLE;
        }
        elsif ( (all { $selection !~ /^$_$/ } @valid_selection_methods) && (all { $experimentId ne $_ } @valid_lib_selection) ){
            warn "\tWarning: for [$libraryId][$experimentId], the library selection is indicated as [$selection], which could be incompatible with our pipeline. Please check (this library was printed in output file).\n";
        }

        # species
        $info =~ /<SCIENTIFIC_NAME>([^<]+)<\/SCIENTIFIC_NAME>/;
        $organism = $1;
        if ( $organism ne $species{ $tsv{'speciesId'}[$i] }->{'organism'} ){
            warn "\tProblem: the organism (scientific name) is not matching between the annotation file [", $species{ $tsv{'speciesId'}[$i] }->{'organism'}, "] and the SRA record [$organism], please verify for [$libraryId][$experimentId]. The information from the annotation file is printed in output file.\n";
        }
        # platform
        $info =~ /<PLATFORM><[^<]+><INSTRUMENT_MODEL>([^<]+)<\/INSTRUMENT_MODEL><\/[^<]+><\/PLATFORM>/;
        $platform = $1;
        if (($platform ne $tsv{'platform'}[$i]) and (!exists $checked_libraries{$libraryId})) {
            warn "\tProblem: the platform is not matching between the annotation file [", $tsv{'platform'}[$i], "] and the SRA record [$platform], please verify for [$libraryId][$experimentId] and update $RNAseqLibChecks. The information from the annotation file is printed in output file.\n";
        }
        # Run IDs
        while ( $info =~ /<PRIMARY_ID>([SEDC]RR\d+)<\/PRIMARY_ID>/g ) {
            push @SRR, $1;
        }
        if (scalar @SRR eq 0) {
            warn "\tProblem: no run ID found for $libraryId.\n";
        }

        # Below, instead of checking info from SRA, we get new info from SRA (if possible)
        # First, library type: single-end or paired-end
        # <LIBRARY_LAYOUT><PAIRED NOMINAL_SDEV="280.488123947368" NOMINAL_LENGTH="400"/></LIBRARY_LAYOUT> # PE
        # <LIBRARY_LAYOUT><SINGLE/></LIBRARY_LAYOUT> # SE
        if ( $info =~ /<LIBRARY_LAYOUT><(\w+)(\s.+)?\/><\/LIBRARY_LAYOUT>/ ) {
            $libraryType = $1;
        }
        else {
            warn "\tProblem: no library type specified for $libraryId.\n";
        }
        # Additional info on library (stranded or not for example), not compulsory
        # <EXPERIMENT_ATTRIBUTE><TAG>library_type</TAG><VALUE>cDNAShotgunReadTwoSense</VALUE></EXPERIMENT_ATTRIBUTE><EXPERIMENT_ATTRIBUTE>
        if ( $info =~ /<TAG>library_type<\/TAG><VALUE>([^<]+?)<\/VALUE>/ ) { # the ? is here to get non greedy matching (there are several <VALUE>...<\/VALUE> tags in the XML files)
            $libraryInfo = $1;
        }
        # Read length, not compulsory
        # Beware, this is the sum of both pairs for paired-end! E.g., SRX1152842, lenght=500, corresponding to 2*250bp
        if ( $info =~ /<SPOT_LENGTH>(\d+)<\/SPOT_LENGTH>/ ) {
            $readLength = $1;
        }

        # If read length defined and it seems very short, issue a warning
        if ( $readLength ne '' ){
            #NOTE Currently those short read length libraries look fine (in those experiments), need to check how read alignments are
            #     GSE36026 GSE38998 DRP000415 ERP000787 SRP003276 SRP003822 SRP003823 SRP003826 SRP003829 SRP003831
            if ( (($libraryType eq 'SINGLE') or ($libraryType eq '')) and ($readLength < 36) ){
                warn "\tInfo: Read length is [$readLength] for SE library [$libraryId][$experimentId], which seems low and could indicate that the library is not a classical RNA-seq library. Please check.\n";
            }
            elsif ( ($libraryType eq 'PAIRED') and ($readLength < 72) ){
                warn "\tInfo: Read length is [$readLength / 2] for PE library [$libraryId][$experimentId], which seems low and could indicate that the library is not a classical RNA-seq library. Please check.\n";
            }
        }

        ## Issue warning is the XML entry includes keywords suggesting that the library is not classical RNA-seq
        #NOTE See https://gitlab.sib.swiss/Bgee/expression-annotations/issues/30
        #FIXME SRP070951 is kept till we know what to do with "globin reduction" see https://gitlab.sib.swiss/Bgee/expression-annotations/issues/30
        my @not_traditional = ('DeepSAGE', 'DeepCAGE', 'LongSAGE', 'SuperSAGE', 'CAGE', 'RACE', 'SAGE', 'DpnII', 'DpnIII', 'NlaIII', 'capture', 'CEL-seq', 'DGE-Seq', 'TagSeq', 'globin reduction', 'globin depletion', 'UMI', 'UMIs', 'ATAC-seq', 'MAINE-Seq', 'Mnase-Seq', 'FAIRE-Seq', 'DNase-seq', 'MACE-Seq', 'QuantSeq');
        my @verified        = ('ERP000787', 'ERP001694', 'ERP013973', 'ERP104395',
                               'GSE22410', 'GSE64283',
                               'SRP000401', 'SRP013825', 'SRP021940', 'SRP022567', 'SRP041131', 'SRP076617', 'SRP082284', 'SRP090001',
                               'SRP091779', 'SRP092799', 'SRP098705', 'SRP099849', 'SRP112616', 'SRP123447', 'SRP125768', 'SRP125959');
        if ( all { $experimentId ne $_ } @verified and any { $info =~ /\W$_\W/i } @not_traditional ){
            my @words = grep { $info =~ /\W$_\W/i } @not_traditional;
            warn "\tWarning: [$libraryId][$experimentId] may not be traditional RNA-seq, found word(s) [".join('/', @words)."] subject to caution. Please check.\n";
        }

        # TODO Any other info to get?
    }
    # Printing output file
    print {$OUT} join("\t", $libraryId,
                            $experimentId,
                            $tsv{'speciesId'}[$i],
                            $species{ $tsv{'speciesId'}[$i] }->{'organism'},
                            $species{ $tsv{'speciesId'}[$i] }->{'genomeFilePath'},
                            $species{ $tsv{'speciesId'}[$i] }->{'database'},
                            $tsv{'platform'}[$i],
                            $libraryType,
                            $libraryInfo,
                            $readLength,
                            join(',', sort @SRR),
                     ), "\n";
}
close $OUT;
print "Done\n";
exit 0;

