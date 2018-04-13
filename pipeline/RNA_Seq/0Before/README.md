* `get_GTF_files.pl`
  Downloads the GTF annotation files for all species in Bgee

* `get_genome_files.pl`
  Downloads the genome fasta files for all species in Bgee

* `prepare_GTF.R`:
  * This scripts reads the gene annotation file (GTF format) and outputs processed file gtf_all
  * Exons from all transcripts are included in the output (coding genes or not)
  * A set of intergenic region is included in the output. For all intergenic regions, 500 bps are removed on eahc size. Regions smaller than 2,000 bps are removed. Those larger than 20,000 bps are shortened, and 10,000 bps are kept on each side of the center of the region.

* `check_rna_seq_curation.pl`
  * Script used to detect errors in the annotation file
  * Checks the presence of organs in the database, and if these organs exists at the developmental stages (database `organ` and `stage` tables need to be filled prior to running the script)
  * The script can be run `before` running the pipeline, but also `after`. The `after` option is not useful anymore and was removed from Makefile (it checked if processed data exist for each RNA-seq library, but wihtout considering the good files for the set of libraries to test). Better testing was put at beginning of `3Insertion/insert_rna_seq.pl` script.
  * TODO there are some improvements to make to this script, for example to check if leading and trailing spaces are present in the annotation fields. See comment in Makefile
  * TODO there are some case sensitivity issues with strain annotation that could be noticed by this script. If two conditions differ only by case sensitivity of strain, the `insert_get_conditions` function will try to reinsert it but SQL, which is not case sensitive) will bug because the condition already exists...


* `create_rna_seq_sample_info.pl`
  * This scripts parses the annotation file, and retrieves information on the species using the Bgee database, and on the SRA records (species, platform, SRR IDs) usign NCBI e-utils
  * The script is launched on our annotation merged wiht Wormbase annotation. An external file is checked for libraries which were manually checked (with decision to include them or not).
  * The script `get_sra_id.pl` from the old pipeline was merged into this script for simplification, and to create less intermediate information files.
  * The script issues warnings if the organism and platform information in the annotation do not match the information on the SRA record (before, TRUE or FALSE flags columns were written in output file). 
  * **Thus the output of this script needs to be verified!**. 
   Some warnings are minor and can be ignored, for example: `Problem: the organism (scientific name) is not matching between the annotation file [Xenopus tropicalis] and the SRA record [Xenopus (Silurana) tropicalis], please verify. The information from the annotation file is printed in output file`.
   Some samples are annotated but from species not yet in Bgee, for example not GSM1054987 is from Acomys cahirinus (rodent) Rat
   Some inconsistencies between annotation and SRA record shall be retained over releases if SRA is wrong. For platform, we remember those cases from release to release, so that we don't have to verify everytime. See issue #98 on annotation GitHub.
   Some warnings are issued for miRNA-seq, ncRNA-seq, and other non-conventional RNA-seq libraries (CAGE-seq, RACE-seq, SAGE, DeepSage, etc). Of course some will be missed because meta-data in SRA are often incomplete... 
   Some warnings are issued if read length too short (<36nt): if the library is recent enough, this is suspicious and could denote a non-conventional RNA-seq library

* `get_SRA.pl`
  * Downloads the .sra files for all libraries in `rna_seq_sample_info.txt`
  * Extracts the .fastq files from the .sra files
  * Compress these files
  * For GTEx data, encrypt the fastq files.
  * `get_SRA-reverseOrder.pl*` is launched on vital-IT. It parses the `rna_seq_sample_info.txt` file in the reverse order, using only odd or even number rows, allowing to launch downloads in parallel on two machines (in addition to bigbgee)
  

