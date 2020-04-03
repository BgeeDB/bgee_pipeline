**Requirements**: having successfully run dblite_creation

## Prerequirements 
[Install Virtuoso 7.2](http://vos.openlinksw.com/owiki/wiki/VOS#How%20Do%20I%20Install%20Virtuoso%3F)

Edit [virtuoso.ini](http://docs.openlinksw.com/virtuoso/dbadm/) configuration file to allow for reading/writing files in the directory that contains the RDF-based files.

virtuoso.ini line to edit: **DirsAllowed = ., ../vad, /new/path/bgee/ttl**

Download and extract [Ontop-cli 3.0.1](https://sourceforge.net/projects/ontop4obda/files/ontop-3.0.1/ontop-cli-3.0.1.zip/download) command line tool. 

**Goals**:
* Set up the Bgee SPARQL endpoint



## Running the easyBgee to RDF converter

Set up the Bgee SPARQL endpoint with easyBgee_v*RELEASE* database and mapping files for Ontop tool.

To run the bgee_rdf.sh bash script:
 easyBgee to RDF triples and load into a Virtuoso data store
                [ {-x | --ontop} <ontop directory> ]
                [ {-m | --mapping-file} <mapping file> ]
                [ {-o | --output-dir-path} <Turtle file output> ]
                [ {-p | --ontop-property-file} <property file>  ]
                [ {-t | --ontology-file} <ontology file> ]
                [ {-i | --isql} <iSQL tool from Virtuoso> ]
                [ {-s | --host-port} <host and port to connect to Virtuoso, e.g. locahost:1111> ]
                [ {-w | --virtuoso-pwd} <Virtuoso password> ]
                [ {-u | --virtuoso-user} <Virtuoso username> ]
                [ {-v | --bgee-version} <the bgee version, e.g. '14_1'> ]
                [ {-h | --help} <usage help> ]

   Example:              
    ./bgee_rdf.sh -m ./conf/genex_adapt.obda -o **/new/path/bgee/ttl** -p  ./conf/genex_adapt.properties -x ./ontop-cli-3 -t ./conf/genex_adapt.owl -i ~/not_save/virtuoso/bin/isql -s 1111 -w dba -u dba -v '14_1'



## Data generation

* If it is the first time you execute this step in this pipeline run:
  `make clean`

* Run Makefile:
  `make`

## Data verification

## Error handling

## Other notable Makefile targets

