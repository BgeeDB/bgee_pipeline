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
ignore_dropletBased_exp = ["SRP467631"]

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

def get_species_names(cursor):
    # define query
    query_all_species="""
    SELECT distinct speciesId, genus, species FROM species
    """
        # Execute the MySQL query
    cursor.execute(query_all_species)
    # Fetch the results of the query
    results = cursor.fetchall()
    speciesId_to_name = {}
    for speciesId, genus, species in results:
        speciesId_to_name[speciesId] = genus.replace(" ", "_") + "_" + species.replace(" ", "_")
    return speciesId_to_name

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

def exp_to_h5ad_full_length(species_ID, expID, name, description, doi, output, species_name, cursor):
    full_length_dir = "{output}/full_length/{species_name}".format(output=output, species_name=species_name)
    if not os.path.exists(full_length_dir):
        os.makedirs(full_length_dir)
    h5ad_file_path = "{full_length_dir}/{species_name}_{expID}_full_length.h5ad".format(full_length_dir=full_length_dir, species_name=species_name, expID=expID)
    # check if file already exist
    if os.path.isfile(h5ad_file_path):
        print("There is already an existing H5ad file for {exp_id} experiment and species {species_id}".format(exp_id=expID, species_id=species_ID))
    else:
        # Define 1st query (retrieve all metadata for one experiment)
        query_per_library = """
        SELECT DISTINCT annots.rnaSeqLibraryAnnotatedSampleId, cond.anatEntityId,
        cond.stageId, cond.cellTypeId, cond.strain, cond.sex,cond.speciesId, annots.rnaSeqLibraryId, anat.anatEntityName, annots.anatEntityAuthorAnnotation, stage.stageName, annots.stageAuthorAnnotation,
        cellType.anatEntityName as cellTypeName, annots.cellTypeAuthorAnnotation, lib.rnaSeqSequencerName,
        lib.cellCompartment, lib.libraryType
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
        libID=[result[7] for result in results]
        query_per_lib = """
        SELECT gene.geneId, result.abundanceUnit, result.abundance, result.readsCount, result.UMIsCount
        FROM rnaSeqLibraryAnnotatedSampleGeneResult AS result
        INNER JOIN gene ON result.bgeeGeneId = gene.bgeeGeneId
        WHERE result.rnaSeqLibraryAnnotatedSampleId = %s
        """

        counts_dict = {}
        for libSamp in tqdm(SampleId):
            cursor.execute(query_per_lib, [libSamp])
            results = cursor.fetchall()
            #print(results)
            for j, result in enumerate(results):
                gene_id = result[0]
                count = result[3]
                tpm = result[2]
                if libSamp not in counts_dict:
                    counts_dict[libSamp] = {}
                if gene_id not in counts_dict[libSamp]:
                    counts_dict[libSamp][gene_id] = {}
                counts_dict[libSamp][gene_id]["count"] = count
                counts_dict[libSamp][gene_id]["tpm"] = tpm
        # Create the count matrix by iterating over the libSam and gene IDs
        unique_gene_ids = list(set(gene_id for libSamp in counts_dict for gene_id in counts_dict[libSamp]))
        count_matrix = []
        tpm_matrix = []
        for libSamp in tqdm(SampleId):
            libSamp_counts = []
            libSamp_tpms = []
            for gene_id in unique_gene_ids:
                if gene_id in counts_dict[libSamp]:
                    libSamp_counts.append(counts_dict[libSamp][gene_id]["count"])
                    libSamp_tpms.append(counts_dict[libSamp][gene_id]["tpm"])
                else:
                    libSamp_counts.append(0)
                    libSamp_tpms.append(0)
            count_matrix.append(libSamp_counts)
            tpm_matrix.append(libSamp_tpms)
        count_matrix= csr_matrix(count_matrix, dtype=np.float32)
        tpm_matrix= csr_matrix(tpm_matrix, dtype=np.float32)
        # Create a dictionary libSamp IDs to metadata values, then transform to df for anndata implementation
        metadata_dict = {SampleId: {"library_id": libID, "anatEntityId": anatEntityId, "anatEntityName": anatEntityName, "anatEntityAuthorAnnotation": anatEntityAuthorAnnotation, "stageId": stageId, "stageName": stageName, "stageAuthorAnnotation": stageAuthorAnnotation, "cellTypeId":cellTypeId, "cellTypeName": cellTypeName, "cellTypeAuthorAnnotation": cellTypeAuthorAnnotation, "strain": strain, "sex":sex, "speciesId":speciesId, "rnaSeqSequencerName":rnaSeqSequencerName, "libraryType": libraryType, "cellCompartment":cellCompartment } for SampleId, libID, anatEntityId, anatEntityName, stageId, stageName, cellTypeId, cellTypeName, strain, sex, speciesId, anatEntityAuthorAnnotation, stageAuthorAnnotation, cellTypeAuthorAnnotation,rnaSeqSequencerName, cellCompartment, libraryType in zip(SampleId, libID, anatEntityId, anatEntityName, stageId, stageName, cellTypeId, cellTypeName, strain, sex, speciesId, anatEntityAuthorAnnotation, stageAuthorAnnotation, cellTypeAuthorAnnotation, rnaSeqSequencerName, cellCompartment, libraryType)}
        metadata_df = pd.DataFrame.from_dict(metadata_dict, orient='index') # index automatically libSamp_ids
        metadata_df.insert(1, 'experiment_id', expID)  # add experiment id column , str(exp_ID[0])
        metadata_df.fillna(value=np.nan, inplace=True)  # replace None with NaN
        # Create a DataFrame with the gene metadata (for anndata implementation)
        gene_metadata = pd.DataFrame({"gene_id": unique_gene_ids})
        #Create anndata object
        adata = ad.AnnData(X=count_matrix, obs=metadata_df, var=gene_metadata) #as metadata dict same order than libSamp_ids from which the count table have been created it's ok
        adata.layers["abundance"] = tpm_matrix
        # Document what each matrix contains
        adata.uns['matrix_descriptions'] = {
            'X': 'Raw read counts',
            'layers': {
                'abundance': 'Transcripts Per Million',
            }
        }
        adata.obs_names = SampleId # comme for loop ils seront dans le bon ordre
        adata.var_names = unique_gene_ids
        adata.write(h5ad_file_path)
        adata.obs.to_csv(h5ad_file_path.replace(".h5ad", ".tsv"), sep="\t", index=True, header=True)

def exp_to_h5ad_dropletBased(species_id, exp_id, name, description, doi, output, species_name, cursor, result_dir, logger):
    """
    Generate an .h5ad file for a droplet-based single-cell RNA-seq experiment by
    combining data from matrices files. Loads UMI count matrices from one or
    more libraries, concatenates them (cells as rows, genes as columns), and writes
    to an AnnData file.
    """
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
               lib.rnaSeqSequencerName, lib.cellCompartment, lib.libraryType
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
            "libraryType": row[17]
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

    # Genes are the same and in the same order for all libraries. We can load them from the first library's gene file.
    gene_file_path = os.path.join(result_dir, library_ids[0], "gene_counts", "gene.genes.txt")
    genes = []
    with open(gene_file_path, 'r') as gf:
        genes = [line.strip() for line in gf]
    n_vars = len(genes)

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

        # Subset the sparse matrix and the barcodes using the indices to keep
        subset_sparse_matrix = sparse.csr_matrix(matrix[indices_to_keep, :])
        subset_barcodes = [barcodes[i] for i in indices_to_keep]
        subset_matrices.append(subset_sparse_matrix)
        # log the max score of the subset_sparse_matrix
        logger.debug(f"Max score of subset_sparse_matrix: {np.max(subset_sparse_matrix.data)}")
        for barcode in subset_barcodes:
            all_barcode_names.append(f"{barcode}_{library_id}")
            obs_metadata.append({
                "barcode": barcode,
                "anatId": library_barcode_info[library_id][barcode]["anatId"],
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
                "libraryType": library_barcode_info[library_id][barcode]["libraryType"]
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
    adata = ad.AnnData(X=final_matrix, dtype=np.float32)
    # now that we created the AnnData object, we can delete the final_matrix to save memory
    del final_matrix
    gc.collect()
    adata.obs = pd.DataFrame(obs_metadata)
    adata.obs_names = all_barcode_names
    gene_metadata = pd.DataFrame({"gene_id": genes})
    adata.var = gene_metadata
    adata.var_names = genes

    # Write the AnnData file on disk
    logger.debug(f"Start to write h5ad file to {h5ad_file_path}")
    adata.write(h5ad_file_path, compression="gzip")
    adata.obs.to_csv(h5ad_file_path.replace(".h5ad",".tsv"), sep="\t", index=True, header=True)
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
    # Get the species name from the database
    species_id_to_name = get_species_names(cursor)
    # Get experiment info for the specified species and experiment
    experiments = return_experiment_ids(species_id, args.exp_id, cursor, logger)

    for exp_id, name, description, doi, species_id, has_full_length, has_droplet in experiments:
        # Process experiments with full-length if has_full_length is 1
        if has_full_length == 1 and exp_id not in ignore_full_length_exp:
            full_length_output_path = output_path / "full_length"
            full_length_output_path.mkdir(parents=True, exist_ok=True)
            print(f"Processing full-length for experiment ID: {exp_id} and species ID: {species_id}")
            exp_to_h5ad_full_length(species_id, exp_id, name, description, doi, full_length_output_path, species_id_to_name[species_id], cursor)
        # Process experiments with droplet-based if has_droplet is 1
        if has_droplet == 1 and exp_id not in ignore_dropletBased_exp:
            droplet_output_path = output_path / "droplet"
            droplet_output_path.mkdir(parents=True, exist_ok=True)
            logger.info(f"Processing droplet-based for experiment ID: {exp_id} and species ID: {species_id}")
            #exp_to_h5ad_dropletBased(species_id, exp_id, name, description, doi, droplet_output_path, species_id_to_name[species_id],
            #                         cursor, args.result_dir, logger)
    # Close the cursor and connection
    cursor.close()
    cnx.close()

if __name__ == "__main__":
    main()
