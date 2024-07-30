include makefile.properties
SHELL = /bin/bash
export CONFIG_TEMPLATE

.PHONY: all clean

main: config_file  running_app 

config_file:
	@echo "$$CONFIG_TEMPLATE" > $(GENERATE_HOMOLOGS_APP_DIR)/template.properties
	@touch $@

install_pipenv:
	@curl https://raw.githubusercontent.com/pypa/pipenv/master/get-pipenv.py | $(PYTHON) 1> $@.tmp 2> $@.err
	@mv $@.tmp $@

running_app:
	@cd $(GENERATE_HOMOLOGS_APP_DIR)
	@$(PIPENV) install 1>$@.tmp 2>$@.err
	@$(PIPENV) run $(PYTHON) $(GENERATE_HOMOLOGS_APP_FILE) -c $(TEMPORARY_CONFIG_FILES_DIRECTORY) -n $(N_PROCESSES) 1>> $@.tmp 2>> $@.err
	@mv $@.tmp $@



