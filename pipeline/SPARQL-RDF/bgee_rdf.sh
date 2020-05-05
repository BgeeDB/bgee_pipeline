#!/bin/bash
usage()
{
    echo "usage: easyBgee to RDF triples and load into a Virtuoso data store
                [ {-x | --ontop} <ontop directory> ]
                [ {-m | --mapping-file} <mapping file> ]
                [ {-o | --output-dir-path} <Turtle file output> ]
                [ {-p | --ontop-property-file} <property file>  ]
                [ {-t | --ontology-file} <ontology file> ]
                [ {-v | --bgee-version} <the bgee version, e.g. '14_1'> ]
                [ {-j | --java-ontop-args } <Java Ontop memory arguments > (optional) ]
                [ {-h | --help} <usage help> ]

   Example:              
    ./bgee_rdf.sh -m ./conf/genex_adapt.obda -o ./ttl -p  ./conf/genex_adapt.properties -x ./ontop-cli-3 -t ./conf/genex_adapt.owl -v '14_1' -j '-Xmx128G'
    "

}
#Variables to be assigned, currently the default is no value.
ontop_dir_path=
obda_file_path=
output_dir_path=
ontop_property_file=
ontology_file=
bgee_version=
java_ontop_args=

if [ "$1" == "" ]
then usage
     exit 1
fi

while [ "$1" != "" ]; do
    case $1 in
        -x | --ontop )          shift
                                ontop_dir_path=${1%/}
                                ;;
        -m | --mapping-file )   shift
                                obda_file_path=${1%/}
                                ;;
        -o | --output-dir-path ) shift
                                 output_dir_path=${1%/}
                                ;;
        -p | --ontop-property-file ) shift
                                     ontop_property_file=${1%/}
                                    ;;  
        -t | --ontology-file )  shift
                                ontology_file=${1%/}
                                ;;
        -v | --bgee-version )   shift
                                bgee_version=${1%/}
                                ;;     
        -j | --java-ontop-args )   shift
                                java_ontop_args=${1%/}
                                ;;                                                                                                          
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

echo "Starting RDF creation and loading into Virtuoso data store..."
echo "Creating RDF data and serializing them as TURTLE files..."
#ontop edit java max memory
if [ "$java_ontop_args" != "" ]
then 
    java_ontop_args=`echo $java_ontop_args | tr -d \"`
    sed -i  's,ONTOP_JAVA_ARGS=.*.,ONTOP_JAVA_ARGS="'$java_ontop_args'",g' $ontop_dir_path/ontop
fi

#Run ontop tool
$ontop_dir_path/ontop materialize --separate-files -m $obda_file_path -f turtle -o $output_dir_path -p $ontop_property_file -t $ontology_file  
#--disable-reasoning 

# for MacOS it has the empty string argument ''
#sed -i '' 's," \.,"^^xsd:double.,g' $output_dir_path/http___purl_org_genex_hasExpressionLevelP_*.ttl 
sed -i 's," \.,"^^xsd:double.,g' $output_dir_path/http___purl_org_genex_hasExpressionLevelP_*.ttl 

echo "The RDF data were created and saved as TURTLE files..."


echo "Generating iSQL script to load RDF data into a virtuoso data store..."
dir_temp=`pwd`

cd $output_dir_path
ls *.ttl > "./ttl_files.txt"
input="./ttl_files.txt"
output="./add_bgee_virtuoso.sql"

while IFS= read -r line
do
  file=`pwd`"/$line"   
  echo "DB.DBA.TTLP_MT(file_to_string_output('$file'), '', 'https://bgee.org/rdf_v$bgee_version');" >>  "$output"
done < "$input"
echo "SPARQL LOAD <http://purl.org/genex#> INTO <https://bgee.org/rdf_v$bgee_version>;" >>  "$output"
echo "SPARQL LOAD <http://purl.org/lscr#> INTO <https://bgee.org/rdf_v$bgee_version>;" >>  "$output"
#Remove empty string statements, if any
echo "SPARQL with <https://bgee.org/rdf_v$bgee_version> delete {?z ?d '' } where {?z ?d ''. };" >>  "$output"
rm -f ./ttl_files.txt
cd $dir_temp

echo "Executing iSQL script to load RDF data into a virtuoso data store..."
$isql_file_path  $virtuoso_host $virtuoso_user $virtuoso_pwd "EXEC=LOAD '$output_dir_path'/add_bgee_virtuoso.sql"

echo "Finished. All files were loaded into a virtuoso data store."
