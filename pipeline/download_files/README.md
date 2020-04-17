# Generation of download files

**General information:**

1. [Rank information](#rank-information)

**Developer guidelines**

1. [Details](#details)
2. [Data generation](#data-generation)
3. [Data verification](#data-verification)
4. [Error handling](#error-handling)
5. [Other notable Makefile targets](#other-notable-makefile-targets)
6. [Notes about Jupyter/iPtython notebooks](#notes-about-jupyter-ipython-notebooks)

# General information

## Rank information

### RNA-Seq data

[See how RNA-Seq ranks are computed](../post_processing/#rna-seq-data)

#### Notes

1. As of Bgee 14.1, ranks are computed in conditions with raw data, considering the anatomical entity
and the developmental stage. It means that all data from different sexes and strains, but the same
anatomical entity and developmental stage, are averaged.
2. The anatomical entities and developmental stages considered are the ones remapped for producing
calls of presence/absence of expression, not the ones used in the raw data annotations. For instance,
if a library is annotated to the developmental stage HsapDv:0000134 "40-year-old human stage",
it is remapped to the more general stage HsapDv:0000090 "25-44 year-old human stage".
Both the raw annotations and the remapped conditions are provided in the download files.
3. The rank assigned to an anatomical entity, not considering the developmental stage,
is the minimum over the ranks computed for the different developmental stages.
Basically, it means that we compute ranks only for the various remapped
anatomical entities-developmental stages, we do not compute ranks in anatomical entities by averaging
all data in it.
4. To compute ranks, we average the ranks of genes in the libraries in a same remapped condition,
weighted by the count of distinct ranks in the corresponding library. We then normalize this weighted mean
by using the max rank in the remapped condition, as compared to the max of the max ranks
over all conditions.
5. As of Bgee 14.1, ranks are not propagated using the graph of conditions, as it is the case for calls
of presence/absence of expression. It means that there is rank information only in conditions with
data annotated. This might change in a future release (for instance, by averaging the ranks
in all sub-conditions with data).

#### Compute expression ranks

The necessary information is present in our processed expression value download files.
Below, `species_name` refers to the scientific name of the species you are interested in.

** 1. Retrieve information from download files **

1. In the folder `species_name_RNA-Seq_experiments_libraries`, open the file `species_name_RNA-Seq_libraries.tsv`. The columns that you will need to use are: `Experiment ID`,
`Library ID`, `Expression mapped anatomical entity ID`, `Expression mapped stage ID`,
`Distinct rank count`, `Max rank in the expression mapped condition`.
2. Identify all the libraries that have the same values for the columns
`Expression mapped anatomical entity ID` and `Expression mapped stage ID`.
3. You can retrieve the raw gene ranks from each library from the files in the folder
`species_name_RNA-Seq_read_counts_TPM_FPKM`. All data are grouped by experiment, so that all libraries
from a same experiment are present in a same file. This means that if the libraries in a same
remapped condition, identified at the previous step, belong to different experiments, you will need
to open these different files. The columns that you will need to use in these files are:
`Experiment ID`, `Library ID`, `Gene ID`, `Rank`.

** 2. compute weighted mean rank of genes for each remapped condition **

1. A remapped condition considers the mapped anatomical entity and the mapped developmental stage.
From the file `species_name_RNA-Seq_libraries.tsv`, identify all the libraries that have
the same values for the columns `Expression mapped anatomical entity ID` and `Expression mapped stage ID`.
2. Retrieve all rank values for a given gene in a same mapped condition from the files in the folder
`species_name_RNA-Seq_read_counts_TPM_FPKM`. You can use the `Experiment ID` of the library
to target the appropriate file (file name format
`species_name_RNA-Seq_read_counts_TPM_FPKM_EXPERIMENTID.tsv.`). In the file, you can use the column
`Library ID` to retrieve the data related to the libraries in the mapped condition. Note:
do not use the `Anatomical entity ID` and `Stage ID` provided in this file, they are not remapped.
3. Compute the weighted mean rank for each gene over all the libraries in the same mapped condition:
`SUM('Rank' * 'Distinct rank count')/SUM('Distinct rank count')`

** 3. normalize the weighted mean ranks **

1. From the column `Max rank in the expression mapped condition` in the file
`species_name_RNA-Seq_libraries.tsv`, retrieve the maximum value over all rows. This will be referred to
as `Max of max ranks` below.
2. For each mapped condition and gene, the weighted mean rank is normalized:
`'Weighted mean rank' * (1 + 'Max of max ranks' / 'Max rank in the expression mapped condition')/2`

** 4. Rank for anatomical entities **

If you want to retrieve the rank for a gene in an anatomical entity, as used in Bgee,
you need to retrieve the minimum normalized weighted mean rank of this gene in this anatomical entity
over all developmental stages (as computed in the previous steps).

#### Compute expression scores from ranks

We also transform these ranks into expression scores, more easily understandable by users:
higher gene expression translates into lower rank but higher expression score, from 0 to 100.

The expression score is computed as
`('Max of max ranks' + 1 - 'Normalized weighted mean rank') * (100/'Max of max ranks')`

# Developer guidelines

1. [Details](#details)
2. [Data generation](#data-generation)
3. [Data verification](#data-verification)
4. [Error handling](#error-handling)
5. [Other notable Makefile targets](#other-notable-makefile-targets)
6. [Notes about Jupyter/iPtython notebooks](#notes-about-jupyter-ipython-notebooks)

## Details

**Requirements**: having successfully run all previous steps, this is a post-processing steps
(TODO: list precisely the final steps required)

**Goal**: generate download files, notebooks, and XRefs file to link from UniProtKB

* `make ../../download_files/generate_processed_files`: generates the processed data download files 
* `make ../../download_files/xref_uniprot_info`: generates file containing XRefs comsumed by UniProt 
to link to Bgee. 

* **IMPORTANT**: You need to update the files `release.tsv` at the root of the Bgee FTP. 
This file is necessary to the BgeeDB R package, and needs to be updated at each release.

## Data generation

* Verify that the IDs of the data sources UniProtKB/TrEMBL and UniProtKB/Swiss-Prot used in the `Makefile` are correct. Check the variables `TREMBL_ID` and `SWISSPROT_ID` in `Makefile`, and compare them to the `dataSourceId` used in [pipeline/db_creation/insert_data_sources.sql](../db_creation/insert_data_sources.sql) for these two data sources.

* If it is the first time you execute this step in this pipeline run:
  `make clean`

* Run Makefile:
  `make`

* After the `Data verification` step (see below):
1. Upload the XRefs for UniProtKB to our [FTP server](ftp://ftp.bgee.org/), to the location `current/XRefBgee.txt` (so that it is accessible at the URL [ftp://ftp.bgee.org/current/XRefBgee.txt](ftp://ftp.bgee.org/current/XRefBgee.txt)). As of Bgee 13, this is done on the machine `ftpbgee`, from the path `/u01/ftp`.
2. Similarly, upload the files containing rank information to FTP.
    1. Copy all zip files in `generated_files/download_files/ranks/anat_entity/` and  `generated_files/download_files/ranks/condition/` to the FTP, so that it is accessible at the URL [ftp://ftp.bgee.org/current/download/ranks/anat_entity/](ftp://ftp.bgee.org/current/download/ranks/anat_entity/) and [ftp://ftp.bgee.org/current/download/ranks/condition/](ftp://ftp.bgee.org/current/download/ranks/condition/). As of Bgee 13, this is done on the machine `ftpbgee`, from the path `/u01/ftp`.
    2. Update the file [generated_files/download_files/ranks/README.md](../../generated_files/download_files/ranks/README.md), notably to provide max rank information for current release, and copy it to FTP to be accessible at [ftp://ftp.bgee.org/current/download/ranks/README.md](ftp://ftp.bgee.org/current/download/ranks/README.md).
3. Upload to FTP the expression call and processed data download files (TO BE WRITTEN)
4. Don't forget to update the file `release.tsv` at the root of the [Bgee FTP](ftp://ftp.bgee.org/).

## Data verification

### For XRefs from UniProtKB:
* Check the file `generated_files/download_files/xref_uniprot_info`
* Check the file `generated_files/download_files/XRefBgee.txt`
* After uploading the new file `XRefBgee.txt` to our FTP at the location `current/XRefBgee.txt`, check that it is correctly accessible from [ftp://ftp.bgee.org/current/XRefBgee.txt](ftp://ftp.bgee.org/current/XRefBgee.txt)

### For download files with rank info:
* Check the file `generated_files/download_files/rank_files_info`. Be notably careful about file sizes.

## Error handling

## Other notable Makefile targets
* generate the XRefs for UniProtKB:
  `make ../../generated_files/download_files/xref_uniprot_info`
* generate the files with rank info:
  `make ../../generated_files/download_files/rank_files_info`


## Notes about Jupyter/iPython notebooks

TODO: to be integrated into pipeline

We use Jupyter/iPython to generate nice-looking figures about the content of the download files (see '.ipynb' files in this directory). Installation instruction to use notebooks:

1. Install with pip command the necessary python modules

    ```
    $ pip install jupyter pandas ipython matplotlib numpy seaborn
    ```
    
2. With the terminal, go to the directory where jupyter notebook file is
3. Run Jupyter Notebook

    ```
    $ jupyter notebook
    ```
    This commnad will open a webpage where jupyter notebook file will appear

4. Double click over jupyter notebook file

    This will open a new webpage with the jupyter notebook.
    Note that you have to change in the jupyter notebook the directory where bgee files are and this files must be uncompressed


Versions of python modules we have been working with:

    ipython==4.0.0
    jupyter==1.0.0
    matplotlib==1.3.1
    numpy==1.8.0rc1
    pandas==0.16.2
    seaborn==0.6.0
