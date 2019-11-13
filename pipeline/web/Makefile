PIPELINEROOT := ../
DIR_NAME := web/
include $(PIPELINEROOT)Makefile.common

all: $(VERIFICATIONFILE)

sitemap.xml:
	@perl generate_site_map.pl -bgee=$(BGEECMD) >out 2>err
	@if [[ ! -s out ]]; then $(RM) out; fi
	@if [[ ! -s err ]]; then $(RM) err; fi
	@echo
	@echo 'Copy "sitemap*.xml" files to the htdocs/ folder, i.e. in "/var/www/html/bgee.unil.ch"'
	@echo 'The file htdocs/robots.txt must contain a line linking to the index sitemap file "sitemap.xml":'
	@echo '        Sitemap: https://bgee.org/sitemap.xml'
	@echo
	@echo 'Think to validate and submit your "sitemap.xml" in Google Dashboard anytime you modify it!'
	@echo


$(VERIFICATIONFILE): sitemap.xml
	@touch $@

clean:
	-@$(RM) -R $(VERIFICATIONFILE) sitemap* out err
