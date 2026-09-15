import psutil, os
import mysql.connector
import pandas as pd
import numpy as np
from scipy.io import mmread
from scipy import sparse
from scipy.sparse import csr_matrix
from pathlib import Path
import anndata as ad
import argparse
import gc
import logging
from tqdm import tqdm

desc="""Generates .h5ad files from single-cell data"""

def setup_logger(log_dir):
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    log_path = os.path.join(log_dir, "h5ad.log")
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_path),
            # logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

# If some experiments are problematic you can put them on these list to avoid processing them.
ignore_full_length_exp = []
ignore_dropletBased_exp = []

# ---------------------------------------------------------------------------
# scFAIR schema support
#
# The .h5ad files we publish must follow the scFAIR schema , see https://github.com/scFAIR/scFAIR_schema.
# Compliance is checked at https://www.sc-fair.org/stats/compliance_all.
#
# The schema reserves a set of obs/var/uns keys. The Bgee-native columns are renamed
# into those keys rather than duplicated, so a column never appears twice under two
# names. Bgee columns with no scFAIR counterpart are kept untouched, which the schema
# allows. Note that the .tsv companion file follows obs, so its header changes too.
#
# ---------------------------------------------------------------------------
SCFAIR_SCHEMA_VERSION = "7.1.0+scfair1.0"
SCFAIR_SCHEMA_REFERENCE = "https://github.com/scFAIR/scFAIR/edit/main/schema/7.1.0/schema.md"

# cond.sex is an enum, so every possible value is mapped here and an unexpected
# value means the enum changed in the schema and this mapping must be updated.
# The schema requires "unknown" when the sex is unavailable, which is what the
# three non-informative Bgee values amount to.
SEX_TO_PATO = {
    "female": ("PATO:0000383", "female"),
    "male": ("PATO:0000384", "male"),
    "hermaphrodite": ("PATO:0001340", "hermaphrodite"),
    "not annotated": ("unknown", "unknown"),
    "mixed": ("unknown", "unknown"),
    "NA": ("unknown", "unknown"),
}
# Mapping of the Bgee protocol annotations onto EFO, keyed by the
# (rnaSeqTechnologyName, sequencedTranscriptPart) pair: the transcript part is what
# separates 10x 3' from 10x 5', which rnaSeqTechnologyName alone does not record.
# Every term below was checked to be a descendant of "EFO:0010183" for single cell
# library construction, or of "EFO:0002772" for assay by molecule, as the schema
# requires. An unknown pair is a curation gap and raises rather than defaulting to
# a generic term, so that new protocols are annotated deliberately.
ASSAY_TO_EFO = {
    ("Smart-Seq2", "full length"): ("EFO:0008931", "Smart-seq2"),
    ("Adapted Smart-Seq2", "full length"): ("EFO:0008931", "Smart-seq2"),
    ("Smart-Seq", "full length"): ("EFO:0008930", "Smart-seq"),
    ("SMARTer Ultra Low", "full length"): ("EFO:0010184", "Smart-like"),
    ("Fluidigm C1 + SMARTer Ultra Low", "full length"):
        ("EFO:0010058", "Fluidigm C1-based SMARTer library preparation"),
    ("Fluidigm C1 instrument and Nextera XT protocol", "full length"):
        ("EFO:0010058", "Fluidigm C1-based SMARTer library preparation"),
    ("C1 autoprep", "full length"):
        ("EFO:0010058", "Fluidigm C1-based SMARTer library preparation"),
    ("10X Genomics V2", "3prime"): ("EFO:0009899", "10x 3' v2"),
    ("10X Genomics V3", "3prime"): ("EFO:0009922", "10x 3' v3"),
    # A generic library prep kit: the single-cell method itself is not recorded,
    # so the most precise honest term is the generic one.
    ("NEBNext Ultra DNA Library Prep Kit for Illumina", "full length"):
        ("EFO:0008913", "single-cell RNA sequencing"),
}
# rnaSeqLibrary.cellCompartment is an enum('NA', 'nucleus', 'cell'), which maps
# one to one onto the values allowed for obs["suspension_type"].
CELL_COMPARTMENT_TO_SUSPENSION_TYPE = {"cell": "cell", "nucleus": "nucleus", "NA": "na"}
# Bgee stores the Ensembl release it used as a dataSource; the schema accepts a
# closed list of database names.
ENSEMBL_DB_NAMES = {"Ensembl": "Ensembl", "EnsemblMetazoa": "Ensembl Metazoa"}
# Bgee uses the root of the cell type ontology as a placeholder when no cell type
# was annotated. The schema requires "unknown" in that case. A NULL cell type,
# on the other hand, is a data error: it MUST NOT happen for single-cell data.
UNANNOTATED_CELL_TYPE_IDS = {"GO:0005575", "UBERON:0000000"}
HUMAN_SPECIES_ID = 9606


def normalize_chromosome(seq_region_name):
    """Return a chromosome name following the scFAIR requirements: no "chr"
    prefix, and the mitochondrial designator always spelled "MT"."""
    if not seq_region_name:
        return "unknown"
    name = str(seq_region_name).strip()
    if name.lower().startswith("chr"):
        name = name[3:]
    # Ensembl spells the mitochondrion differently per species: "MT" for most,
    # "MtDNA" for C. elegans, "mitochondrion_genome" for Drosophila...
    if name.upper() in ("M", "MT", "MTDNA") or "mitochondr" in name.lower():
        return "MT"
    return name


def get_species_info(cursor):
    """Return, per species ID, the name used for output paths plus everything
    the scFAIR uns fields need (organism, genome assembly, Ensembl release)."""
    cursor.execute("""
    SELECT DISTINCT sp.speciesId, sp.genus, sp.species, sp.genomeVersion,
           ds.dataSourceName, ds.releaseVersion
    FROM species AS sp
    LEFT JOIN dataSource AS ds ON sp.dataSourceId = ds.dataSourceId
    """)
    species_info = {}
    for (species_id, genus, species, genome_version,
         data_source_name, release_version) in cursor.fetchall():
        species_info[species_id] = {
            "name": genus.replace(" ", "_") + "_" + species.replace(" ", "_"),
            "organism": f"{genus} {species}",
            "organism_ontology_term_id": f"NCBITaxon:{species_id}",
            "feature_reference": f"NCBITaxon:{species_id}",
            "ensembl_assembly": genome_version or "unknown",
            "ensembl_database": ENSEMBL_DB_NAMES.get(data_source_name, "Ensembl"),
            "ensembl_release": str(release_version) if release_version else "unknown",
        }
    return species_info


# Gene annotations are the same for every experiment of a species, and the
# transcript join is expensive, so they are queried once per species.
gene_metadata_per_species = {}


def get_gene_metadata(cursor, species_id, logger):
    """Return the scFAIR var annotations of every gene of a species, keyed by
    Ensembl gene ID."""
    if species_id in gene_metadata_per_species:
        return gene_metadata_per_species[species_id]
    cursor.execute("""
    SELECT gene.geneId, gene.geneName, gene.seqRegionName, bioType.geneBioTypeName, gene.geneLength
    FROM gene
    LEFT JOIN geneBioType AS bioType ON gene.geneBioTypeId = bioType.geneBioTypeId
    WHERE gene.speciesId = %s
    """, [species_id])
    gene_metadata_per_species[species_id] = {
        gene_id: {
            "feature_name": gene_name if gene_name else gene_id,
            "feature_type": biotype if biotype else "unknown",
            "feature_chromosome": normalize_chromosome(seq_region_name),
            # gene.geneLength is the median of the lengths of the isoforms of the gene,
            # which is what the schema asks for. It is NULL for the species inserted
            # from a non-Ensembl source; the schema has no "unknown" value for a
            # length, so 0 flags those.
            "feature_length": gene_length if gene_length else 0,
        }
        for gene_id, gene_name, seq_region_name, biotype, gene_length in cursor.fetchall()
    }
    return gene_metadata_per_species[species_id]


def build_scfair_var(gene_ids, gene_metadata, species_info, logger):
    """Build the var DataFrame required by the schema for the given gene IDs."""
    missing = [gene_id for gene_id in gene_ids if gene_id not in gene_metadata]
    if missing:
        logger.warning(f"{len(missing)} feature(s) are absent from the Bgee gene table "
                       f"(e.g. {missing[:5]}); their var annotations default to placeholders.")
    without_length = [gene_id for gene_id in gene_ids
                      if gene_metadata.get(gene_id, {}).get("feature_length", 0) == 0]
    if without_length:
        logger.warning(f"{len(without_length)} of {len(gene_ids)} feature(s) have no geneLength in "
                       f"Bgee (e.g. {without_length[:5]}); their feature_length is set to 0, which "
                       "does not comply with scFAIR. geneLength is filled by the genes pipeline, "
                       "and stays NULL for species inserted from a non-Ensembl source.")
    var = pd.DataFrame({
        "feature_name": [gene_metadata.get(g, {}).get("feature_name", g) for g in gene_ids],
        "feature_type": [gene_metadata.get(g, {}).get("feature_type", "unknown") for g in gene_ids],
        "feature_chromosome": [gene_metadata.get(g, {}).get("feature_chromosome", "unknown")
                               for g in gene_ids],
        "feature_length": [gene_metadata.get(g, {}).get("feature_length", 0) for g in gene_ids],
        # Bgee only distributes real genes, never ERCC spike-ins, and never filters
        # genes out of the matrix it publishes.
        "feature_biotype": "gene",
        "feature_reference": species_info["feature_reference"],
        "feature_is_filtered": False,
    }, index=gene_ids)
    for column in ("feature_type", "feature_chromosome", "feature_biotype", "feature_reference"):
        var[column] = var[column].astype("category")
    return var


def to_scfair_obs(obs, species_info, experiment_id):
    """Rename the Bgee-native obs columns into the names the schema reserves, and
    derive the fields Bgee does not store as such. Expects the unified column names
    anatEntityId / rnaSeqLibraryId (see the callers). Columns that have no scFAIR
    counterpart (the author annotations, the sequencer, the barcode...) are kept as
    they are: the schema allows additional metadata."""
    # Derive the columns that are not a plain rename first, while their source
    # columns are still around.
    if obs["cellTypeId"].isna().any() or (obs["cellTypeId"].astype(str).str.strip() == "").any():
        raise ValueError(f"Experiment {experiment_id} has annotated samples without a cell type. "
                         "Every single-cell annotated sample in Bgee must have one.")
    obs["cell_type_ontology_term_id"] = [
        "unknown" if cell_type_id in UNANNOTATED_CELL_TYPE_IDS else cell_type_id
        for cell_type_id in obs["cellTypeId"]
    ]
    obs["cell_type"] = [
        "unknown" if cell_type_id in UNANNOTATED_CELL_TYPE_IDS else cell_type_name
        for cell_type_id, cell_type_name in zip(obs["cellTypeId"], obs["cellTypeName"])
    ]
    protocols = list(zip(obs["rnaSeqTechnologyName"], obs["sequencedTranscriptPart"]))
    unmapped_protocols = {protocol for protocol in protocols if protocol not in ASSAY_TO_EFO}
    if unmapped_protocols:
        raise ValueError(f"No EFO term mapped for protocol(s) {sorted(unmapped_protocols)}. "
                         "Add them to ASSAY_TO_EFO.")
    obs["assay_ontology_term_id"] = [ASSAY_TO_EFO[protocol][0] for protocol in protocols]
    obs["assay"] = [ASSAY_TO_EFO[protocol][1] for protocol in protocols]
    unexpected_sexes = set(obs["sex"]) - set(SEX_TO_PATO)
    if unexpected_sexes:
        raise ValueError(f"Unexpected value(s) {sorted(unexpected_sexes)} in cond.sex. "
                         "Update SEX_TO_PATO to cover the whole enum.")
    # The schema reserves "sex" for the label of the PATO term, so the Bgee value is
    # replaced in place. It is kept in its own column because the schema collapses
    # 'mixed' and 'not annotated' into the same "unknown".
    obs["bgeeSex"] = obs["sex"]
    obs["sex_ontology_term_id"] = [SEX_TO_PATO[sex][0] for sex in obs["bgeeSex"]]
    obs["sex"] = [SEX_TO_PATO[sex][1] for sex in obs["bgeeSex"]]
    obs["suspension_type"] = [
        CELL_COMPARTMENT_TO_SUSPENSION_TYPE.get(compartment, "na")
        for compartment in obs["cellCompartment"]
    ]

    # Plain renames, then drop the columns the derived ones consumed.
    obs.rename(columns={
        "anatEntityId": "tissue_ontology_term_id",
        "anatEntityName": "tissue",
        "stageId": "development_stage_ontology_term_id",
        "stageName": "development_stage",
        # Free text in Bgee, so the ontologized counterpart cannot be filled.
        "strain": "strain_or_genetic_background",
    }, inplace=True)
    # cellCompartment maps one to one onto suspension_type, and the cell type columns
    # only differ from their scFAIR counterparts by the unannotated placeholder, so
    # nothing is lost by dropping them. rnaSeqTechnologyName and sequencedTranscriptPart
    # are kept: several of them share one EFO term, so assay alone loses the protocol.
    obs.drop(columns=["cellTypeId", "cellTypeName", "cellCompartment"], inplace=True)
    # Bgee does not track individuals; a library is the finest-grained proxy we have.
    obs["donor_id"] = obs["rnaSeqLibraryId"].astype(str)

    # Fields Bgee has no column for.
    obs["experiment_id"] = experiment_id
    # Bgee only annotates samples taken from a tissue, never cell lines or organoids.
    obs["tissue_type"] = "tissue"
    # Bgee only integrates samples from healthy, wild-type-like individuals.
    obs["disease_ontology_term_id"] = "PATO:0000461"
    obs["disease"] = "normal"
    # Bgee holds no self-reported ancestry as such, and the schema mandates "na"
    # outside of Homo sapiens.
    ethnicity = "unknown" if species_info["organism_ontology_term_id"] == \
        f"NCBITaxon:{HUMAN_SPECIES_ID}" else "na"
    obs["self_reported_ethnicity_ontology_term_id"] = ethnicity
    obs["self_reported_ethnicity"] = ethnicity
    obs["is_primary_data"] = True

    categorical_columns = ["physiologicalStatus", "bgeeSex", "rnaSeqTechnologyName",
                           "sequencedTranscriptPart", "assay_ontology_term_id", "assay",
                           "tissue_type", "tissue_ontology_term_id", "tissue",
                           "cell_type_ontology_term_id", "cell_type",
                           "development_stage_ontology_term_id", "development_stage",
                           "sex_ontology_term_id", "sex", "disease_ontology_term_id", "disease",
                           "self_reported_ethnicity_ontology_term_id", "self_reported_ethnicity",
                           "strain_or_genetic_background", "suspension_type"]
    for column in categorical_columns:
        obs[column] = obs[column].astype("category")
    return obs


def set_scfair_uns(adata, species_info, experiment_id, name, doi, assay_description):
    """Fill the dataset-level metadata reserved by the schema."""
    adata.uns["schema_version"] = SCFAIR_SCHEMA_VERSION
    adata.uns["schema_reference"] = SCFAIR_SCHEMA_REFERENCE
    adata.uns["organism_ontology_term_id"] = species_info["organism_ontology_term_id"]
    adata.uns["organism"] = species_info["organism"]
    adata.uns["ensembl_release"] = species_info["ensembl_release"]
    adata.uns["ensembl_database"] = species_info["ensembl_database"]
    adata.uns["ensembl_assembly"] = species_info["ensembl_assembly"]
    # The title has to be unique across datasets, hence the assay and species.
    adata.uns["title"] = (f"{name if name else experiment_id} - {species_info['organism']} "
                          f"({assay_description}, {experiment_id})")
    if doi:
        adata.uns["citation"] = f"https://doi.org/{doi}" if not str(doi).startswith("http") else doi
    # Libraries are the batches Bgee integrates within an experiment.
    adata.uns["batch_condition"] = ["rnaSeqLibraryId"]

def get_args():
    """Parse the arguments """
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--output_dir", help="Location of the directory where to save output file")
    parser.add_argument("--exp_id", help="Specific SRP/ERP/DRP experiment ID to process (optional)", default=None)
    parser.add_argument("--species_id", help="NCBI id of the desired species (optional). If no species ID is provided then H5AD files are created for all Bgee species", default=0, type=int)
    parser.add_argument("--server", help="server address e.g rbioinfo.unil.ch", default=0)
    parser.add_argument("--db", help="mysql database, e.g bgee_15_h5ad", default=0)
    parser.add_argument("--usr", help=" mysql user name, e.g bgee", default=0)
    parser.add_argument("--pwd", help="password of the mysql user, e.g bgee", default=0)
    parser.add_argument("--result_dir", help="directory containing sparse matrices for all libraries", default=0)
    parser.add_argument("--intergenic_prefixes", help="comma-separated list of prefixes used for intergenic regions in gene files (e.g. upstream_,downstream_)", default="upstream_,downstream_")

    args = parser.parse_args()
    # Check if all required arguments are provided
    required_args = ["output_dir", "server", "db", "usr", "pwd", "result_dir"]
    missing_args = [arg for arg in required_args if not getattr(args, arg)]

    if missing_args:
        print("Error: The following arguments are required:")
        for arg in missing_args:
            print(f"--{arg}")
        parser.print_help()
        exit(1)
    return parser.parse_args()

def return_experiment_ids(species_id, exp_id, cursor, logger):
    # define query
    query_all_exps="""
    SELECT lib.rnaSeqExperimentId, exp.rnaSeqExperimentName, exp.rnaSeqExperimentDescription, exp.DOI, cond.speciesId,
    MAX(CASE WHEN annots.multipleLibraryIndividualSample = 0 THEN 1 ELSE 0 END) AS hasFullLength,
    MAX(CASE WHEN annots.multipleLibraryIndividualSample = 1 THEN 1 ELSE 0 END) AS hasDropletBased
    FROM cond
    INNER JOIN rnaSeqLibraryAnnotatedSample AS annots ON annots.conditionId = cond.conditionId
    INNER JOIN rnaSeqLibrary AS lib ON lib.rnaSeqLibraryId = annots.rnaSeqLibraryId
    INNER JOIN rnaSeqExperiment AS exp ON lib.rnaSeqExperimentId = exp.rnaSeqExperimentId
    WHERE lib.rnaSeqTechnologyIsSingleCell = 1
    GROUP BY lib.rnaSeqExperimentId, cond.speciesId
    ORDER BY lib.rnaSeqExperimentId, cond.speciesId;
    """
        # Execute the MySQL query
    cursor.execute(query_all_exps)
    # Fetch the results of the query
    results = cursor.fetchall()
     # subset results only for targeted species:
    filtered_results = results.copy()
    if species_id:
        filtered_results = [result for result in results if result[4] == species_id]
    # subset results only for targeted experiment:
    if exp_id:
        filtered_results = [result for result in filtered_results if result[0] == exp_id]
    logger.info(f"Number of experiments/species to process: {len(filtered_results)}")
    return filtered_results

def exp_to_h5ad_full_length(species_ID, expID, name, description, doi, output, species_info, cursor, logger):
    species_name = species_info["name"]
    full_length_dir = "{output}/{species_name}".format(output=output, species_name=species_name)
    if not os.path.exists(full_length_dir):
        os.makedirs(full_length_dir)
    h5ad_file_path = "{full_length_dir}/{species_name}_{expID}_full_length.h5ad".format(full_length_dir=full_length_dir, species_name=species_name, expID=expID)
    # check if file already exist
    if os.path.isfile(h5ad_file_path):
        logger.info(f"There is already an existing H5ad file for {expID} experiment and species {species_ID}")
    else:
        # Define 1st query (retrieve all metadata for one experiment)
        query_per_library = """
        SELECT DISTINCT annots.rnaSeqLibraryAnnotatedSampleId, cond.anatEntityId,
        cond.stageId, cond.cellTypeId, cond.strain, cond.sex,cond.speciesId, annots.rnaSeqLibraryId, anat.anatEntityName, annots.anatEntityAuthorAnnotation, stage.stageName, annots.stageAuthorAnnotation,
        cellType.anatEntityName as cellTypeName, annots.cellTypeAuthorAnnotation, lib.rnaSeqSequencerName,
        lib.cellCompartment, lib.libraryType, annots.physiologicalStatus,
        lib.rnaSeqTechnologyName, lib.sequencedTranscriptPart
        FROM rnaSeqLibraryAnnotatedSample AS annots
        INNER JOIN rnaSeqLibrary AS lib ON annots.rnaSeqLibraryId = lib.rnaSeqLibraryId
        INNER JOIN cond ON cond.conditionId = annots.conditionId
        INNER JOIN anatEntity AS anat ON cond.anatEntityId = anat.anatEntityId
        INNER JOIN anatEntity AS cellType ON cond.cellTypeId = cellType.anatEntityId
        INNER JOIN stage as stage ON stage.stageId = cond.stageId
        WHERE annots.multipleLibraryIndividualSample = 0 AND lib.rnaSeqTechnologyIsSingleCell = 1
        AND lib.rnaSeqExperimentId =  %s AND cond.speciesId = %s;
        """
        # Execute the MySQL query
        cursor.execute(query_per_library, (expID, species_ID))
        #print('Total Row(s):', cursor.rowcount) # nbre of libraries
        results = cursor.fetchall()
        # Extract metadata from results
        SampleId=[str(result[0]) for result in results]
        anatEntityId=[result[1] for result in results]
        stageId=[result[2] for result in results]
        cellTypeId=[result[3] for result in results]
        strain=[result[4] for result in results]
        sex=[result[5] for result in results]
        speciesId=[result[6] for result in results]
        anatEntityName=[result[8] for result in results]
        anatEntityAuthorAnnotation =[result[9] for result in results]
        stageName=[result[10] for result in results]
        stageAuthorAnnotation=[result[11] for result in results]
        cellTypeName=[result[12] for result in results]
        cellTypeAuthorAnnotation = [result[13] for result in results]
        rnaSeqSequencerName = [result[14] for result in results]
        cellCompartment=[result[15] for result in results]
        libraryType = [result[16] for result in results]
        physiologicalStatus = [result[17] for result in results]
        rnaSeqTechnologyName = [result[18] for result in results]
        sequencedTranscriptPart = [result[19] for result in results]
        libID=[result[7] for result in results]
        query_per_lib = """
        SELECT gene.geneId, result.readsCount
        FROM rnaSeqLibraryAnnotatedSampleGeneResult AS result
        INNER JOIN gene ON result.bgeeGeneId = gene.bgeeGeneId
        WHERE result.rnaSeqLibraryAnnotatedSampleId = %s
        """

        counts_dict = {}
        for libSamp in tqdm(SampleId):
            cursor.execute(query_per_lib, [libSamp])
            results = cursor.fetchall()
            #print(results)
            for gene_id, count in results:
                if libSamp not in counts_dict:
                    counts_dict[libSamp] = {}
                counts_dict[libSamp][gene_id] = count
        # Create the count matrix by iterating over the libSam and gene IDs
        unique_gene_ids = list(set(gene_id for libSamp in counts_dict for gene_id in counts_dict[libSamp]))
        count_matrix = []
        for libSamp in tqdm(SampleId):
            count_matrix.append([counts_dict[libSamp].get(gene_id, 0) for gene_id in unique_gene_ids])
        count_matrix= csr_matrix(count_matrix, dtype=np.float32)
        # Create a dictionary libSamp IDs to metadata values, then transform to df for anndata implementation
        metadata_dict = {SampleId: {"library_id": libID, "anatEntityId": anatEntityId, "anatEntityName": anatEntityName, "anatEntityAuthorAnnotation": anatEntityAuthorAnnotation, "stageId": stageId, "stageName": stageName, "stageAuthorAnnotation": stageAuthorAnnotation, "cellTypeId":cellTypeId, "cellTypeName": cellTypeName, "cellTypeAuthorAnnotation": cellTypeAuthorAnnotation, "strain": strain, "sex":sex, "speciesId":speciesId, "rnaSeqSequencerName":rnaSeqSequencerName, "libraryType": libraryType, "cellCompartment":cellCompartment, "physiologicalStatus": physiologicalStatus, "rnaSeqTechnologyName": rnaSeqTechnologyName, "sequencedTranscriptPart": sequencedTranscriptPart } for SampleId, libID, anatEntityId, anatEntityName, stageId, stageName, cellTypeId, cellTypeName, strain, sex, speciesId, anatEntityAuthorAnnotation, stageAuthorAnnotation, cellTypeAuthorAnnotation,rnaSeqSequencerName, cellCompartment, libraryType, physiologicalStatus, rnaSeqTechnologyName, sequencedTranscriptPart in zip(SampleId, libID, anatEntityId, anatEntityName, stageId, stageName, cellTypeId, cellTypeName, strain, sex, speciesId, anatEntityAuthorAnnotation, stageAuthorAnnotation, cellTypeAuthorAnnotation, rnaSeqSequencerName, cellCompartment, libraryType, physiologicalStatus, rnaSeqTechnologyName, sequencedTranscriptPart)}
        metadata_df = pd.DataFrame.from_dict(metadata_dict, orient='index') # index automatically libSamp_ids
        metadata_df.fillna(value=np.nan, inplace=True)  # replace None with NaN
        # Use the same column names as the droplet-based path, so that a single
        # function can add the scFAIR obs fields for both.
        metadata_df.rename(columns={"library_id": "rnaSeqLibraryId"}, inplace=True)
        metadata_df = to_scfair_obs(metadata_df, species_info, expID)
        # Create a DataFrame with the gene metadata (for anndata implementation)
        gene_metadata = build_scfair_var(unique_gene_ids, get_gene_metadata(cursor, species_ID, logger),
                                         species_info, logger)
        #Create anndata object
        adata = ad.AnnData(X=count_matrix, obs=metadata_df, var=gene_metadata) #as metadata dict same order than libSamp_ids from which the count table have been created it's ok
        set_scfair_uns(adata, species_info, expID, name, doi, "full-length")
        # Only the raw counts are published: rnaSeqLibraryAnnotatedSample.abundanceUnit
        # is an enum('tpm','cpm'), so a single "abundance" layer could not be labelled
        # with one unit without lying about the libraries using the other one.
        adata.uns['matrix_descriptions'] = {'X': 'Raw read counts'}
        adata.obs_names = SampleId # comme for loop ils seront dans le bon ordre
        adata.var_names = unique_gene_ids
        # Same as for droplet-based data: write to temporary paths and rename them only once both
        # files are complete, so that a run killed mid-write does not leave a truncated .h5ad that
        # the "file already exist" check at the top of this function would skip on the next run.
        tsv_file_path = h5ad_file_path.replace(".h5ad", ".tsv")
        tmp_h5ad_file_path = h5ad_file_path + ".tmp"
        tmp_tsv_file_path = tsv_file_path + ".tmp"
        adata.write(tmp_h5ad_file_path)
        adata.obs.to_csv(tmp_tsv_file_path, sep="\t", index=True, header=True)
        os.replace(tmp_h5ad_file_path, h5ad_file_path)
        os.replace(tmp_tsv_file_path, tsv_file_path)

def exp_to_h5ad_dropletBased(species_id, exp_id, name, description, doi, output, species_info, cursor, result_dir, intergenic_prefixes, logger):
    """
    Generate an .h5ad file for a droplet-based single-cell RNA-seq experiment by
    combining data from matrices files. Loads UMI count matrices from one or
    more libraries, concatenates them (cells as rows, genes as columns), and writes
    to an AnnData file.
    """
    species_name = species_info["name"]
    species_dir = os.path.join(output, species_name)
    if not os.path.exists(species_dir):
        os.makedirs(species_dir)
    h5ad_file_path = os.path.join(species_dir, f"{species_name}_{exp_id}_droplet_based.h5ad")
    if os.path.isfile(h5ad_file_path):
        logger.info(f"H5AD file for experiment {exp_id} already exists at {h5ad_file_path}. Skipping.")
        return
    # Query retrieving metadata for all barcodes in the experiment
    #XXX To optimize memory usage, we could run that query per library. We could also retrieve the metadata per annotatedSampleId,
    #    and separatly the mapping between 1. libraryId and annotedSampleId, 2. barcode and annotatedSampleId
    query_per_cell = """
        SELECT DISTINCT indivs.barcode, annots.rnaSeqLibraryAnnotatedSampleId, cond.anatEntityId,
               cond.stageId, cond.cellTypeId, cond.strain, cond.sex, cond.speciesId, annots.rnaSeqLibraryId,
               anat.anatEntityName, annots.anatEntityAuthorAnnotation,
               stage.stageName, annots.stageAuthorAnnotation,
               cellType.anatEntityName AS cellTypeName, annots.cellTypeAuthorAnnotation,
               lib.rnaSeqSequencerName, lib.cellCompartment, lib.libraryType,
               annots.physiologicalStatus, lib.rnaSeqTechnologyName, lib.sequencedTranscriptPart
        FROM rnaSeqLibraryIndividualSample AS indivs
        INNER JOIN rnaSeqLibraryAnnotatedSample AS annots
            ON indivs.rnaSeqLibraryAnnotatedSampleId = annots.rnaSeqLibraryAnnotatedSampleId
        INNER JOIN rnaSeqLibrary AS lib
            ON annots.rnaSeqLibraryId = lib.rnaSeqLibraryId
        INNER JOIN cond
            ON cond.conditionId = annots.conditionId
        INNER JOIN anatEntity AS anat
            ON cond.anatEntityId = anat.anatEntityId
        INNER JOIN anatEntity AS cellType
            ON cond.cellTypeId = cellType.anatEntityId
        INNER JOIN stage AS stage
            ON stage.stageId = cond.stageId
        WHERE annots.multipleLibraryIndividualSample = 1
          AND lib.rnaSeqTechnologyIsSingleCell = 1
          AND lib.rnaSeqExperimentId = %s AND cond.speciesId = %s
    """
    cursor.execute(query_per_cell, (exp_id, species_id))
    metadata_results = cursor.fetchall()

    # Build the dictionary of metadata
    library_barcode_info = {}
    for row in metadata_results:
        library_id = row[8]
        barcode = row[0]
        if library_id not in library_barcode_info:
            library_barcode_info[library_id] = {}
        # barcodes are unique per library but not per experiment. Thats why
        # the first key of the dictionary is the library id
        library_barcode_info[library_id][barcode] = {
            "anatId": row[2],
            "cellTypeId": row[4],
            "stageId": row[3],
            "strain": row[5],
            "sex": row[6],
            "speciesId": row[7],
            "rnaSeqLibraryId": row[8],
            "anatEntityName": row[9],
            "anatEntityAuthorAnnotation": row[10],
            "stageName": row[11],
            "stageAuthorAnnotation": row[12],
            "cellTypeName": row[13],
            "cellTypeAuthorAnnotation": row[14],
            "rnaSeqSequencerName": row[15],
            "cellCompartment": row[16],
            "libraryType": row[17],
            "physiologicalStatus": row[18],
            "rnaSeqTechnologyName": row[19],
            "sequencedTranscriptPart": row[20]
        }
    del metadata_results
    gc.collect()
    # commented out because it is not needed for now but it coult be necessary in the future
    # if we use an approach that write directly on disk
    #XXX We already tried to write the h5ad file in backed mode.
    #    The idea was to first create an empty csr sparse matrix and then fill it with the data
    #    The problem we faced was that it is not possible to  write a sparse matrix in backed mode because it is
    #    not possible to change sparsality of the matrix. It has to be a dense matrix.
    #    The optimum implementation in terms of memory usage would be to create a Zarr dense matrix
    #    and then write into it per library. For a large matrix (4M barcodes and 80K genes) it would
    #    require 1.2TB of disk space. We did not have such amount of disk space available so we decided not
    #    to implement that option but it will be necessary in the future. We did not test the amount of memory
    #    required to load the Zarr matrix and save it as a sparse matrix.
    # query_num_barcodes = """
    #     SELECT COUNT(DISTINCT indivs.rnaSeqLibraryIndividualSampleId)
    #     FROM rnaSeqLibraryIndividualSample AS indivs
    #     INNER JOIN rnaSeqLibraryAnnotatedSample AS annots
    #     ON indivs.rnaSeqLibraryAnnotatedSampleId = annots.rnaSeqLibraryAnnotatedSampleId
    #     INNER JOIN rnaSeqLibrary AS lib
    #     ON annots.rnaSeqLibraryId = lib.rnaSeqLibraryId
    #     INNER JOIN cond
    #     ON annots.conditionId = cond.conditionId
    #     WHERE lib.rnaSeqTechnologyIsSingleCell = 1
    #     AND lib.rnaSeqExperimentId = %s AND cond.speciesId = %s
    # """
    # cursor.execute(query_num_barcodes, (exp_id, species_id))
    # num_barcodes_result = cursor.fetchone()
    # num_annotated_barcodes = num_barcodes_result[0] if num_barcodes_result else 0
    library_ids = sorted(library_barcode_info.keys())

    # Genes are the same and in the same order for all libraries within the same experiment/species,
    # because all libraries are processed with the same kallisto index. We load from the first library
    # and assert consistency across all others.
    gene_file_path = os.path.join(result_dir, library_ids[0], "gene_counts", "gene.genes.txt")
    all_genes = []
    with open(gene_file_path, 'r') as gf:
        all_genes = [line.strip() for line in gf]
    # It is not safe to only assess the number of genes and gene order, because some libraries may have been processed with a different kallisto index. To ensure that all libraries are processed with the same kallisto index, we check that the gene files are identical across all libraries.
    for lib_id in library_ids[1:]:
        other_gene_file = os.path.join(result_dir, lib_id, "gene_counts", "gene.genes.txt")
        with open(other_gene_file, 'r') as gf:
            other_genes = [line.strip() for line in gf]
        if other_genes != all_genes:
            raise ValueError(f"Gene file mismatch between library {library_ids[0]} and {lib_id} for experiment {exp_id}. "
                             "All libraries must be processed with the same kallisto index.")
    # Filter out intergenic regions using prefixes defined in Makefile.common (INTERGENIC_PREFIXES)
    gene_mask = np.array([not any(g.startswith(p) for p in intergenic_prefixes) for g in all_genes])
    genes = [g for g, keep in zip(all_genes, gene_mask) if keep]
    n_vars = len(genes)
    # A malformed --intergenic_prefixes value (e.g. INTERGENIC_PREFIXES containing a space, in which
    # case argparse only receives the first prefix) would keep every feature and silently ship the
    # intergenic regions in the h5ad file. Fail loudly instead.
    if n_vars == len(all_genes):
        raise ValueError(f"No intergenic region matched the prefixes {intergenic_prefixes} in {gene_file_path}. "
                         "Check the value of INTERGENIC_PREFIXES in Makefile.common.")
    logger.debug(f"Loaded {len(all_genes)} features, kept {n_vars} genes after removing intergenic regions.")

    # init variables retrieved per library
    obs_metadata = []
    all_barcode_names = []
    #final_matrix = None
    subset_matrices = []

    # Loop through each library and process the data
    for library_id in library_ids:
        logger.debug(f"Processing library: {library_id}")
        barcodes = []
        barcode_file_path = os.path.join(result_dir, library_id, "gene_counts", "gene.barcodes.txt")
        with open(barcode_file_path, 'r') as gf:
            barcodes = [line.strip() for line in gf]

        annotated_barcodes = list(library_barcode_info[library_id].keys())

        # load the sparse matrix with scipy
        matrix_file_path = os.path.join(result_dir, library_id, "gene_counts", "gene.mtx")
        matrix = mmread(matrix_file_path).tolil()

        # Create a mask for the rows to keep
        mask = np.isin(barcodes, annotated_barcodes)

        # Get the indices of the barcodes to keep using the mask
        indices_to_keep = np.where(mask)[0]

        # Subset the sparse matrix and the barcodes using the indices to keep, and filter intergenic columns
        subset_sparse_matrix = sparse.csr_matrix(matrix[indices_to_keep, :][:, gene_mask])
        subset_barcodes = [barcodes[i] for i in indices_to_keep]
        subset_matrices.append(subset_sparse_matrix)
        # log the max score of the subset_sparse_matrix
        logger.debug(f"Max score of subset_sparse_matrix: {np.max(subset_sparse_matrix.data)}")
        for barcode in subset_barcodes:
            all_barcode_names.append(f"{barcode}_{library_id}")
            obs_metadata.append({
                "barcode": barcode,
                "anatEntityId": library_barcode_info[library_id][barcode]["anatId"],
                "cellTypeId": library_barcode_info[library_id][barcode]["cellTypeId"],
                "stageId": library_barcode_info[library_id][barcode]["stageId"],
                "strain": library_barcode_info[library_id][barcode]["strain"],
                "sex": library_barcode_info[library_id][barcode]["sex"],
                "speciesId": library_barcode_info[library_id][barcode]["speciesId"],
                "rnaSeqLibraryId": library_id,
                "anatEntityName": library_barcode_info[library_id][barcode]["anatEntityName"],
                "anatEntityAuthorAnnotation": library_barcode_info[library_id][barcode]["anatEntityAuthorAnnotation"],
                "stageName": library_barcode_info[library_id][barcode]["stageName"],
                "stageAuthorAnnotation": library_barcode_info[library_id][barcode]["stageAuthorAnnotation"],
                "cellTypeName": library_barcode_info[library_id][barcode]["cellTypeName"],
                "cellTypeAuthorAnnotation": library_barcode_info[library_id][barcode]["cellTypeAuthorAnnotation"],
                "rnaSeqSequencerName": library_barcode_info[library_id][barcode]["rnaSeqSequencerName"],
                "cellCompartment": library_barcode_info[library_id][barcode]["cellCompartment"],
                "libraryType": library_barcode_info[library_id][barcode]["libraryType"],
                "physiologicalStatus": library_barcode_info[library_id][barcode]["physiologicalStatus"],
                "rnaSeqTechnologyName": library_barcode_info[library_id][barcode]["rnaSeqTechnologyName"],
                "sequencedTranscriptPart": library_barcode_info[library_id][barcode]["sequencedTranscriptPart"]
            })
        # if final_matrix is None:
        #     final_matrix = subset_sparse_matrix
        # else:
        #     final_matrix = sparse.vstack([final_matrix, subset_sparse_matrix], format="csr", dtype=np.float32)
        del subset_sparse_matrix, matrix, subset_barcodes
        gc.collect()
    logger.debug(f"Start to concatenante all sparse matrices")
    # Concatenate all sparse matrices
    logger.debug(f"Memory usage before creation of final matrix: {psutil.Process(os.getpid()).memory_info().rss / (1024**3):.2f} GB")

    final_matrix = sparse.vstack(subset_matrices, format="csr", dtype=np.float32)
    logger.debug(f"Memory usage when final matrix has been created: {psutil.Process(os.getpid()).memory_info().rss / (1024**3):.2f} GB")

    del subset_matrices
    gc.collect()
    logger.debug(f"Finished to concatenantte all sparse matrices")
    # Create the AnnData object
    logger.debug(f"Start to generate h5ad object with {len(all_barcode_names)} observations and {n_vars} variables.")
    # log the memory usage in GB
    logger.debug(f"Memory usage before creating AnnData object: {psutil.Process(os.getpid()).memory_info().rss / (1024**3):.2f} GB")
    adata = ad.AnnData(X=final_matrix)
    # now that we created the AnnData object, we can delete the final_matrix to save memory
    del final_matrix
    gc.collect()
    obs = pd.DataFrame(obs_metadata, index=all_barcode_names)
    adata.obs = to_scfair_obs(obs, species_info, exp_id)
    adata.obs_names = all_barcode_names
    adata.var = build_scfair_var(genes, get_gene_metadata(cursor, species_id, logger), species_info, logger)
    adata.var_names = genes
    set_scfair_uns(adata, species_info, exp_id, name, doi, "droplet-based")

    # Write the AnnData file on disk. We write to temporary paths and rename them only once both
    # files are complete: these matrices are large enough to hit memory limits, and a run killed
    # mid-write would otherwise leave a truncated .h5ad at the final path, which the "file already
    # exists" check at the top of this function would silently skip on the next run.
    tsv_file_path = h5ad_file_path.replace(".h5ad", ".tsv")
    tmp_h5ad_file_path = h5ad_file_path + ".tmp"
    tmp_tsv_file_path = tsv_file_path + ".tmp"
    logger.debug(f"Start to write h5ad file to {h5ad_file_path}")
    adata.write(tmp_h5ad_file_path, compression="gzip")
    adata.obs.to_csv(tmp_tsv_file_path, sep="\t", index=True, header=True)
    os.replace(tmp_h5ad_file_path, h5ad_file_path)
    os.replace(tmp_tsv_file_path, tsv_file_path)
    del adata
    gc.collect()
    logger.info(f"AnnData .h5ad file saved to {h5ad_file_path}")
    # log memory usage in GB
    logger.debug(f"Memory usage after writing AnnData h5ad file: {psutil.Process(os.getpid()).memory_info().rss / (1024**3):.2f} GB")

def main():
    args = get_args()
    # Ensure the output directory exists
    output_path = Path(args.output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    logger = setup_logger(str(output_path))

    # Connect to the MySQL database
    cnx = mysql.connector.connect(
    host=args.server,
    user=args.usr,
    password=args.pwd,
    database=args.db)
    # Create a cursor object to interact with the database
    cursor = cnx.cursor()
    # Get the species ID from the command line argument
    species_id = args.species_id
    # Get the species name and the scFAIR dataset-level metadata from the database
    species_info_by_id = get_species_info(cursor)
    # Get experiment info for the specified species and experiment
    experiments = return_experiment_ids(species_id, args.exp_id, cursor, logger)
    intergenic_prefixes = [p.strip() for p in args.intergenic_prefixes.split(",") if p.strip()]

    for exp_id, name, description, doi, species_id, has_full_length, has_droplet in experiments:
        # Process experiments with full-length if has_full_length is 1
        if has_full_length == 1 and exp_id not in ignore_full_length_exp:
            full_length_output_path = output_path / "full_length"
            full_length_output_path.mkdir(parents=True, exist_ok=True)
            logger.info(f"Processing full-length for experiment ID: {exp_id} and species ID: {species_id}")
            exp_to_h5ad_full_length(species_id, exp_id, name, description, doi, full_length_output_path,
                                    species_info_by_id[species_id], cursor, logger)
        # Process experiments with droplet-based if has_droplet is 1
        if has_droplet == 1 and exp_id not in ignore_dropletBased_exp:
            droplet_output_path = output_path / "droplet"
            droplet_output_path.mkdir(parents=True, exist_ok=True)
            logger.info(f"Processing droplet-based for experiment ID: {exp_id} and species ID: {species_id}")
            exp_to_h5ad_dropletBased(species_id, exp_id, name, description, doi, droplet_output_path,
                                    species_info_by_id[species_id], cursor, args.result_dir, intergenic_prefixes, logger)
    # Close the cursor and connection
    cursor.close()
    cnx.close()

if __name__ == "__main__":
    main()
