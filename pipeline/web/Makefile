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
	@echo '        Sitemap: https://www.bgee.org/sitemap.xml'
	@echo
	@echo 'Think to validate and submit your "sitemap.xml" in Google Dashboard anytime you modify it!'
	@echo


#TODO update sphinx.conf file and regenerate sphinx indexes


schema.org_dump:
	@perl generate_schema.org_dump.pl  -bgee=$(BGEECMD)  >$(DBNAME).schema_org.jsonld 2>schema.org_dump.err
	@if [[ ! -s schema.org_dump.err ]]; then $(RM) schema.org_dump.err; fi
	@gzip -9 $(DBNAME).schema_org.jsonld
	@md5sum $(DBNAME).schema_org.jsonld.gz >$(DBNAME).schema_org.jsonld.gz.md5
	@touch $@


$(VERIFICATIONFILE): sitemap.xml schema.org_dump
	@touch $@

clean:
	-@$(RM) -R $(VERIFICATIONFILE) sitemap* schema.org_dump* $(DBNAME).schema_org.jsonld*
