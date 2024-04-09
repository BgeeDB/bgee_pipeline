include properties
SHELL = /bin/bash
export PARAMETERS

.PHONY: all clean

##Run install_pipenv, if pipenv is not installed.
oma_bgee_xref: config_file
	@$(PIPENV) install 1>$@.tmp 2>$@.err
	@$(PIPENV) run $(PYTHON) $(APP_FILE) 1>> $@.tmp 2>> $@.err
	@mv $@.tmp $@

config_file:
	@echo "$$PARAMETERS" > $(APP_DIR)/parameters.py
	@touch $@

install_pipenv:
	@curl https://raw.githubusercontent.com/pypa/pipenv/master/get-pipenv.py | $(PYTHON) 1> $@.tmp 2> $@.err
	@mv $@.tmp $@