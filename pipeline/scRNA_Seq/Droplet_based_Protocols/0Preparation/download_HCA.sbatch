#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --account=mrobinso_bgee

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=3-00:00:00

#SBATCH --output=/users/smoretti/bgee_pipeline.git/pipeline/scRNA_Seq/Droplet_based_Protocols/download_HCA.out
#SBATCH --error=/users/smoretti/bgee_pipeline.git/pipeline/scRNA_Seq/Droplet_based_Protocols/download_HCA.err
#SBATCH --export=NONE
#SBATCH --job-name=HCA
#SBATCH --mail-user sara.fonsecacosta@unil.ch

## Note: The download is done to /tmp in axiom, since hard links are not allowed in /scratch
## Note: To download HCA data we need to use the manifest file for each experiment (file generated in HCA web-page)
## This script was done with some recommendations from Etienne Orliac (IT service - unil)

source  /dcsrsoft/spack/bin/setup_dcsrsoft
module load gcc python py-virtualenv
python -V

## Paths
export manifest_file=../../../source_files/scRNA_Seq/Manifest_file.tsv
export tmp_folder_Download_data=/tmp/DOWNLOAD_HCA_DATA.$RANDOM
export final_destination=/work/FAC/FBM/DEE/mrobinso/bgee/downloads/scRNA_Seq_All/scRNASeq_libraries_Droplet_10X/

virtualenv hca
source hca/bin/activate
pip install hca
hca --version


## Create a temporary folder in /tmp area
if [ -d $tmp_folder_Download_data ]; then
    rm -r $tmp_folder_Download_data || exit 1
fi
mkdir -v $tmp_folder_Download_data


# Download & process the data
hca dss download-manifest --manifest $manifest_file --download-dir $tmp_folder_Download_data --replica 'aws' --layout bundle --no-metadata


## Read manifest file line by line to identify which directories to move
declare -A dirs_to_move
while read -r el1 el2 the_rest; do
    dirs_to_move[$el1.$el2]=1
done <<< "$(sed 1d $manifest_file)"
echo; echo "List of unique directories = ${!dirs_to_move[@]}"; echo


for dir_name in "${!dirs_to_move[@]}"; do

    echo; echo "@@@ Dealing with unique directory: $dir_name"

    ## delete if already present in final destination folder
    if [[ -d $final_destination/$dir_name ]]; then
	echo "Warning: deleting $final_destination/$dir_name"
	rm -rv $final_destination/$dir_name
    fi
    ## copy the new data
    mv -v $tmp_folder_Download_data/$dir_name $final_destination

done

## delete temporary folder
rm -r $tmp_folder_Download_data

