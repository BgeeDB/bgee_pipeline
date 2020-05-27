**Requirements**: having successfully run dblite_creation

## Prerequirements 
[Install Virtuoso 7.2](http://vos.openlinksw.com/owiki/wiki/VOS#How%20Do%20I%20Install%20Virtuoso%3F)

Edit [virtuoso.ini](http://docs.openlinksw.com/virtuoso/dbadm/) configuration file to allow for reading/writing files in the directory that contains the RDF-based files.

virtuoso.ini line to edit: **DirsAllowed = ., ../vad, /new/path/bgee/ttl**

**Goals**:
* Set up the Bgee SPARQL endpoint


## Data generation

* If it is the first time you execute this step in this pipeline run:
  `make clean`

* Run Makefile:
  `make`


## Running the easyBgee to RDF converter

Generate triples to integrate in the triple store (virtuoso) with easyBgee_v*RELEASE* database and mapping files for Ontop tool.

To run the bgee_rdf.sh bash script:
 easyBgee to RDF triples
 
                [ {-x | --ontop} <ontop directory> ]
                
                [ {-m | --mapping-file} <mapping file> ]
                
                [ {-o | --output-dir-path} <Turtle file output> ]
                
                [ {-p | --ontop-property-file} <property file>  ]
                
                [ {-t | --ontology-file} <ontology file> ]
                                
                [ {-v | --bgee-version} <the bgee version, e.g. '14_1'> ]
                
                [ {-h | --help} <usage help> ]

   Example:              
    ./bgee_rdf.sh -m ./conf/genex_adapt.obda -o **/new/path/bgee/ttl** -p  ./conf/genex_adapt.properties -x ./ontop-cli-3 -t ./conf/genex_adapt.owl -v '14_1'



## Data verification

## Error handling

## Other notable Makefile targets

