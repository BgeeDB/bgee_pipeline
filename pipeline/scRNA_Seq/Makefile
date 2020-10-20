PIPELINEROOT := ../
DIR_NAME := scRNA_Seq/
include $(PIPELINEROOT)Makefile.common

## Run the full length and target-based protocols
run_all_pipelines:
	make -f Full_Length_Protocols/Makefile
	make -f Droplet_based_Protocols/Makefile

