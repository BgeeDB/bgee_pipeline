**Requirements**: having successfully run all previous steps, this is a post-processing steps
(TODO: list precisely the final steps required)

**Goal**: generate download files, notebooks, and XRefs file to link from UniProtKB

## Details

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


## Notes about Jupyter/iPtython notebooks

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
