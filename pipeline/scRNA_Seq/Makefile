PIPELINEROOT := ../
DIR_NAME := scRNA_Seq/
include $(PIPELINEROOT)Makefile.common

## Run the full length and target-based protocols
run_all_pipelines:
  make bgee_pipeline/pipeline/scRNA_Seq/Full_Length_Protocols/Makefile
  make bgee_pipeline/pipeline/scRNA_Seq/Droplet_based_Protocols/Makefile