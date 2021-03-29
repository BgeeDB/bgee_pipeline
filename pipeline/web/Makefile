PIPELINEROOT := ../
DIR_NAME := web/
include $(PIPELINEROOT)Makefile.common

all: $(VERIFICATIONFILE)

sitemap.xml:
	@perl generate_site_map.pl -bgee=$(BGEECMD) >sitemap.out 2>sitemap.err
	@if [[ ! -s sitemap.out ]]; then $(RM) sitemap.out; fi
	@if [[ ! -s sitemap.err ]]; then $(RM) sitemap.err; fi
	@echo
	@echo 'Copy "sitemap*.xml" files to the htdocs/ folder, i.e. in "/var/www/html/bgee.unil.ch"'
	@echo 'The file htdocs/robots.txt must contain a line linking to the index sitemap file "sitemap.xml":'
	@echo '        Sitemap: https://bgee.org/sitemap.xml'
	@echo
	@echo 'Think to validate and submit your "sitemap.xml" in Google Dashboard anytime you modify it!'
	@echo


#TODO update sphinx.conf file and regenerate sphinx indexes

$(VERIFICATIONFILE): sitemap.xml
	@touch $@

clean:
	-@$(RM) -R $(VERIFICATIONFILE) sitemap*
