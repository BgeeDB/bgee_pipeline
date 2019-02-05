PIPELINEROOT := ../
DIR_NAME := SPARQL-RDF/
include $(PIPELINEROOT)Makefile.common
RELEASE := 14.0
FTP_RELEASE := `echo $(RELEASE)| sed -e 's@\.@_@'`

all: $(VERIFICATIONFILE)

get_bgeeLite_dump:
	# Download bgee lite dump
	@$(WGET) "ftp://ftp.bgee.org/bgee_v${FTP_RELEASE}/sql_lite_dump.tar.gz" && $(MV) sql_lite_dump.tar.gz /tmp/ || rm -f sql_lite_dump.tar.gz
#Do not store sql_lite_dump.tar.gz in source_files/ because big and in our FTP
	@tar xvfz /tmp/sql_lite_dump.tar.gz && rm -f /tmp/sql_lite_dump.tar.gz
	@touch $@


$(VERIFICATIONFILE): sitemap.xml
	@touch $@

clean:
	-@$(RM) -R $(VERIFICATIONFILE) get_bgeeLite_dump
