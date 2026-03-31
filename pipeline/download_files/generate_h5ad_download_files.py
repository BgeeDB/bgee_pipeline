#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 17:26:28 2023

@author: sagane
"""

desc="""Generates .h5ad files from single-cell data"""

import pandas as pd
import mysql.connector
import numpy as np
from scipy.sparse import csr_matrix
import anndata as ad
from tqdm import tqdm
import argparse
import os


# Generates .h5ad files from single-cell experiments, supporting both droplet-based and full-length sequencing data.
#Two modes: can either process a specific experiment from a specific species or all the single cell experiments for a given species:
# A third mode has been added. It allows to generate the h5ad files for one exp_id and one species_id. It is useful to generate the h5ad file for multi species experiments if we want to generate the h5ad file only for one of them. In that case the user has to provide the experiment_id and the species_id as argument of the script.
# 1. Processing a Specific Experiment by providing it's ID (--species_id is still mandatory here, to allow the processing of the libraries belonging only to this species)
# python3 cells_to_h5ad_v4.pyy /path/to/output/directory --exp_id [experiment_ID] --species_id [species_id] --server [db_server] --db [database_name] --usr [username] --pwd [password]

# 2. Processing all Experiments by Species ID
#python3 cells_to_h5ad_v4.py /path/to/output/directory --species_id [species_ID] --server [db_server] --db [database_name] --usr [username] --pwd [password]


#server="dbbioinfo.unil.ch"
#db="bgee_v15_2"
#pwd="bgee"
#usr="bgee"
#species_ID=10090
#output="/Users/sagane/Desktop/test"
#output_dir="/Users/sagane/Desktop/test"

# If some experiments are problematic you can put them on these list to avoid processing them.
ignore_full_length_exp = ["SRP201320"]
ignore_dropletBased_exp = []


def get_args():
    """Parse the arguments """
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("output_dir", help="Location of the directory where to save output file")
    parser.add_argument("--exp_id", help="Specific SRP/ERP/DRP experiment ID to process (optional)", default=None)
    parser.add_argument("--species_id", help="NCBI id of the desired species (optional). If no species ID is provided then H5AD files are created for all Bgee species", default=0, type=int)
    parser.add_argument("--server", help="server address e.g rbioinfo.unil.ch", default=0)
    parser.add_argument("--db", help="mysql database, e.g bgee_15_h5ad", default=0)
    parser.add_argument("--usr", help=" mysql user name, e.g bgee", default=0)
    parser.add_argument("--pwd", help="password of the mysql user, e.g bgee", default=0)
    args = parser.parse_args()
    # Check if all required arguments are provided
    required_args = ["output_dir", "server", "db", "usr", "pwd"]
    missing_args = [arg for arg in required_args if not getattr(args, arg)]

    if missing_args:
        print("Error: The following arguments are required:")
        for arg in missing_args:
            print(f"--{arg}")
        parser.print_help()
        exit(1)
    return args

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

def return_expIDs(species_ID, exp_ID, cursor):
    # define query
    query_all_exps="""
    SELECT lib.rnaSeqExperimentId, exp.rnaSeqExperimentName, exp.rnaSeqExperimentDescription, exp.DOI, cond.speciesId,
      IF(EXISTS(SELECT 1
        FROM cond AS cond2
        INNER JOIN rnaSeqLibraryAnnotatedSample AS annots2 ON annots2.conditionId = cond2.conditionId
        INNER JOIN rnaSeqLibrary AS lib2 ON lib2.rnaSeqLibraryId = annots2.rnaSeqLibraryId
        WHERE lib2.rnaSeqLibraryId = lib.rnaSeqLibraryId AND cond2.speciesId = cond.speciesId
        AND annots2.multipleLibraryIndividualSample = 0 LIMIT 1)
      , 1, 0) AS hasFullLength,
      IF(EXISTS(SELECT 1
        FROM cond AS cond2
        INNER JOIN rnaSeqLibraryAnnotatedSample AS annots2 ON annots2.conditionId = cond2.conditionId
        INNER JOIN rnaSeqLibrary AS lib2 ON lib2.rnaSeqLibraryId = annots2.rnaSeqLibraryId
        WHERE lib2.rnaSeqLibraryId = lib.rnaSeqLibraryId AND cond2.speciesId = cond.speciesId
        AND annots2.multipleLibraryIndividualSample = 1 LIMIT 1)
      , 1, 0) AS hasDropletBased
    FROM cond
    INNER JOIN rnaSeqLibraryAnnotatedSample AS annots ON annots.conditionId = cond.conditionId
    INNER JOIN rnaSeqLibrary AS lib ON lib.rnaSeqLibraryId = annots.rnaSeqLibraryId
    INNER JOIN rnaSeqExperiment as exp on lib.rnaSeqExperimentId = exp.rnaSeqExperimentId
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
    if species_ID:
        filtered_results = [result for result in results if result[4] == species_ID]
    # subset results only for targeted experiment:
    if exp_ID:
        filtered_results = [result for result in filtered_results if result[0] == exp_ID]
    print("Number of experiments/species to process: ", len(filtered_results))
    return filtered_results

def exp_to_h5ad_full_length(species_ID, expID, name, description, doi, output, species_id_to_name, cursor):
    species_name = species_id_to_name[species_ID]
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
        #test = adata[adata.obs["stageId"] == "FBdv:00007026"]
        #test.X





def exp_to_h5ad_dropletBased(species_ID, expID, name, description, doi, output, species_id_to_name, cursor):
    species_name = species_id_to_name[species_ID]
    droplet_based_dir = "{output}/droplet_based/{species_name}".format(output=output, species_name=species_name)
    if not os.path.exists(droplet_based_dir):
        os.makedirs(droplet_based_dir)
    h5ad_file_path = "{droplet_based_dir}/{species_name}_{expID}_droplet_based.h5ad".format(droplet_based_dir=droplet_based_dir, species_name=species_name, expID=expID)
    # 1. check is file already exist
    if os.path.isfile(h5ad_file_path):
        print("There is already an existing H5ad file for {exp_id} experiment and species {species_id}".format(exp_id=expID, species_id=species_ID))
    else:
        # Define 1st query (retrieve all metadata for one experiment)
        query_per_cell = """
        SELECT DISTINCT indivs.barcode, annots.rnaSeqLibraryAnnotatedSampleId, cond.anatEntityId,
        cond.stageId, cond.cellTypeId, cond.strain, cond.sex,cond.speciesId, annots.rnaSeqLibraryId, anat.anatEntityName, annots.anatEntityAuthorAnnotation, stage.stageName, annots.stageAuthorAnnotation,
        cellType.anatEntityName as cellTypeName, annots.cellTypeAuthorAnnotation, lib.rnaSeqSequencerName,
        lib.cellCompartment, lib.libraryType
        FROM rnaSeqLibraryIndividualSample AS indivs
        INNER JOIN rnaSeqLibraryAnnotatedSample AS annots ON indivs.rnaSeqLibraryAnnotatedSampleId = annots.rnaSeqLibraryAnnotatedSampleId
        INNER JOIN rnaSeqLibrary AS lib ON annots.rnaSeqLibraryId = lib.rnaSeqLibraryId
        INNER JOIN cond ON cond.conditionId = annots.conditionId
        INNER JOIN anatEntity AS anat ON cond.anatEntityId = anat.anatEntityId
        INNER JOIN anatEntity AS cellType ON cond.cellTypeId = cellType.anatEntityId
        INNER JOIN stage as stage ON stage.stageId = cond.stageId
        WHERE annots.multipleLibraryIndividualSample = 1 AND lib.rnaSeqTechnologyIsSingleCell = 1
        AND lib.rnaSeqExperimentId = %s AND cond.speciesId = %s
        """
        # Execute the MySQL query
        cursor.execute(query_per_cell, (expID, species_ID))
        #print('Total Row(s):', cursor.rowcount)
        # Fetch the results of the query
        results = cursor.fetchall()
        # Extract data from results
        barcodes = [result[0] for result in results]
        SampleId=[str(result[1]) for result in results]
        anatEntityId=[result[2] for result in results]
        stageId=[result[3] for result in results]
        cellTypeId=[result[4] for result in results]
        strain=[result[5] for result in results]
        sex=[result[6] for result in results]
        speciesId=[result[7] for result in results]
        anatEntityName=[result[9] for result in results]
        anatEntityAuthorAnnotation =[result[10] for result in results]
        stageName=[result[11] for result in results]
        stageAuthorAnnotation=[result[12] for result in results]
        cellTypeName=[result[13] for result in results]
        cellTypeAuthorAnnotation = [result[14] for result in results]
        rnaSeqSequencerName = [result[15] for result in results]
        cellCompartment=[result[16] for result in results]
        libraryType = [result[17] for result in results]
        libID=[result[8] for result in results]
        #create unique cell ids that will also allow to feed next query
        unique_cell_ids = list(map('_'.join, zip(barcodes, SampleId))) #ls3 = [x+'-'+y for x in ls1 for y in ls2]
        user_unique_cell_ids= list(map('_'.join, zip(barcodes, libID))) #create unique cell ids for users: replacing our internal 'rnaSeqLibraryAnnotatedSampleId' with SRA library ID
        ####### 2nd Query #######
        query_per_cbarcode = """
        SELECT gene.geneId, cell.abundanceUnit, cell.abundance, cell.readsCount, cell.UMIsCount
        FROM rnaSeqLibraryIndividualSampleGeneResult AS cell
        INNER JOIN gene ON cell.bgeeGeneId = gene.bgeeGeneId
        INNER JOIN rnaSeqLibraryIndividualSample AS indivs ON indivs.rnaSeqLibraryIndividualSampleId = cell.rnaSeqLibraryIndividualSampleId
        WHERE indivs.rnaSeqLibraryAnnotatedSampleId = %s
        AND indivs.barcode = %s;
        """
        # Initialize lists to store CSR matrix data
        rows = []
        cols = []
        data_count = []
        data_cpm = []

        # Create a mapping from gene_id to column index
        unique_gene_ids = [] # keep track of order for var_names
        gene_to_col = {}

        for cell_index, cell_id in enumerate(tqdm(unique_cell_ids)):  # Enumerate for row index
            bar, samp = cell_id.split("_")
            cursor.execute(query_per_cbarcode, (samp, bar))
            results = cursor.fetchall()
            for result in results:
                gene_id = result[0]
                count = result[4]
                cpm = result[2]

                if gene_id not in gene_to_col:
                    gene_to_col[gene_id] = len(unique_gene_ids)
                    unique_gene_ids.append(gene_id) # add in order of appearance

                col_index = gene_to_col[gene_id]
                rows.append(cell_index)
                cols.append(col_index)
                data_count.append(count)
                data_cpm.append(cpm)

        # Create the CSR matrices
        count_matrix = csr_matrix((data_count, (rows, cols)), shape=(len(unique_cell_ids), len(unique_gene_ids)), dtype=np.float32)
        cpm_matrix = csr_matrix((data_cpm, (rows, cols)), shape=(len(unique_cell_ids), len(unique_gene_ids)), dtype=np.float32)

        # Create a dictionary that maps cell IDs to metadata values, then transform to df for anndata implementation
        metadata_dict = {unique_cell_ids: {"barcodes":barcodes, "library_id": libID, "anatEntityId": anatEntityId, "anatEntityName": anatEntityName, "anatEntityAuthorAnnotation": anatEntityAuthorAnnotation, "stageId": stageId, "stageName": stageName, "stageAuthorAnnotation": stageAuthorAnnotation, "cellTypeId":cellTypeId, "cellTypeName": cellTypeName, "cellTypeAuthorAnnotation": cellTypeAuthorAnnotation, "strain": strain, "sex":sex, "speciesId":speciesId, "rnaSeqSequencerName":rnaSeqSequencerName, "libraryType": libraryType, "cellCompartment":cellCompartment } for unique_cell_ids, libID, barcodes, anatEntityId, anatEntityName, stageId, stageName, cellTypeId, cellTypeName, strain, sex, speciesId, anatEntityAuthorAnnotation, stageAuthorAnnotation, cellTypeAuthorAnnotation,rnaSeqSequencerName, cellCompartment, libraryType in zip(unique_cell_ids, libID, barcodes, anatEntityId, anatEntityName, stageId, stageName, cellTypeId, cellTypeName, strain, sex, speciesId, anatEntityAuthorAnnotation, stageAuthorAnnotation, cellTypeAuthorAnnotation, rnaSeqSequencerName, cellCompartment, libraryType)}
        metadata_df = pd.DataFrame.from_dict(metadata_dict, orient='index') # index automatically unique_cell_ids
        metadata_df.insert(1, 'experiment_id', expID)  # add experiment id column , str(exp_ID[0])
        # Create a DataFrame with the gene metadata (for anndata implementation)
        gene_metadata = pd.DataFrame({"gene_id": unique_gene_ids})
        #Create anndata object
        adata = ad.AnnData(X=count_matrix, obs=metadata_df, var=gene_metadata) #as metadata dict same order than unique_cell_ids from which the count table have been created it's ok
        adata.layers["abundance"] = cpm_matrix
        # Document what each matrix contains
        adata.uns['matrix_descriptions'] = {
            'X': 'raw UMIs count',
            'layers': {
                'abundance': 'Counts Per Million',
            }
        }
        adata.obs_names = user_unique_cell_ids # comme for loop ils seront dans le bon ordre
        adata.var_names = unique_gene_ids
        adata.write(h5ad_file_path)
        adata.obs.to_csv(h5ad_file_path.replace(".h5ad",".tsv"), sep="\t", index=True, header=True)
        #test = adata[adata.obs["stageId"] == "FBdv:00007026"]
        #test.X
        #adata[["AAAGATGGTTTCCACC_SRX5979983", "ACACCCTCAGGGTATG_SRX5979983"], ["FBgn0000064", "FBgn0000017"]].X.toarray().tolist()

def main():
    args = get_args()
    # Ensure the output directory exists
    output_path = args.output_dir
    if not os.path.exists(output_path):
        os.makedirs(output_path)
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
    experiments = return_expIDs(species_id, args.exp_id, cursor)

    for expID, name, description, doi, species_id, has_full_length, has_droplet in experiments:
        # Process experiments with full-length if has_full_length is 1
        if has_full_length == 1 and expID not in ignore_full_length_exp:
            print(f"Processing full-length for experiment ID: {expID} and species ID: {species_id}")
            exp_to_h5ad_full_length(species_id, expID, name, description, doi, output_path, species_id_to_name, cursor)
        # Process experiments with droplet-based if has_droplet is 1
        if has_droplet == 1 and expID not in ignore_dropletBased_exp:
            print(f"Processing droplet-based for experiment ID: {expID} and species ID: {species_id}")
            # exp_to_h5ad_dropletBased(species_id, expID, name, description, doi, output_path, species_id_to_name, cursor)
    # Close the cursor and connection
    cursor.close()
    cnx.close()

if __name__ == "__main__":
    main()




# To DO:
# For now I store only one count matrix (= read count), do we want other parallell count matrix with expression level and expression p-values etc?

#- if ERP002 covers one species but two datatype (full-length and droplet based)we also want two files ? ERP002_mouse_full_length and ERP002_mouse_dropletBased
# Two ways of running script, either providing species id, will in this case generate h5ad file for every experiments having this species (only libraries from this specific species, if two differents species, the librairies of the other species will not be included in the h5ad file).
# Or by providing a specific experiment ids (always for a specific species) that will generate a h5ad file containing all the librarires of that specific species (need to check with Fred if my SQL query is correct)
#bien verifier que ca marche quand plusieurs especes !!
