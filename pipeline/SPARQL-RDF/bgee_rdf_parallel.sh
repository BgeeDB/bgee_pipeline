#!/bin/bash
usage()
{
    echo "usage: easyBgee to RDF triples with several Ontop processes working in parallel.
                The mapping is split by slice of genes (see slice_mapping.pl)
                [ {-x | --ontop} <ontop directory> ]
                [ {-m | --mapping-file} <mapping file> ]
                [ {-o | --output-dir-path} <Turtle file output> ]
                [ {-p | --ontop-property-file} <property file>  ]
                [ {-t | --ontology-file} <ontology file> ]
                [ {-c | --xml-catalog} <XML catalog of the imported ontologies, to read them from local files> (optional) ]
                [ {-g | --gene-bounds} <first bgeeGeneId of each slice, comma-separated and ascending> ]
                [ {-d | --delay} <seconds between the start of two Ontop processes> (optional, default 30) ]
                [ {-j | --java-ontop-args } <Java Ontop memory arguments > (optional) ]
                [ {-l | --test-genes} <\"low high\" bgeeGeneId> (optional, test mode) ]
                [ {-n | --slices} <number of slices in test mode> (optional, default 3) ]
                [ {-h | --help} <usage help> ]

   One Ontop process reads each slice of genes, plus one for the other tables. The rule generate_triples
   of the Makefile computes bounds giving slices with about the same number of rows of globalExpression.

   Example:
    ./bgee_rdf_parallel.sh -m ./conf/genex_adapt.obda -o ./ttl -p ./conf/genex_adapt.properties -x ./ontop-cli-3 -t ./conf/genex_adapt.owl -g '257972,410322,523117' -j '-Xmx16G'

   Test mode: only the genes with low <= bgeeGeneId < high are materialized, e.g. a hundred genes.
   They are materialized in one go and in -n slices, and the triples obtained are compared.
   Nothing is written to the output directory but the working files.
    ./bgee_rdf_parallel.sh -m ... -o ./ttl_test -l '300000 300100' -n 3
    "
}
#Variables to be assigned, currently the default is no value.
ontop_dir_path=
obda_file_path=
output_dir_path=
ontop_property_file=
ontology_file=
xml_catalog=
gene_bounds=
slice_count=3
delay=30
java_ontop_args=
test_genes=

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
        -c | --xml-catalog )    shift
                                xml_catalog=${1%/}
                                ;;
        -g | --gene-bounds )    shift
                                gene_bounds=$1
                                ;;
        -n | --slices )         shift
                                slice_count=$1
                                ;;
        -d | --delay )          shift
                                delay=$1
                                ;;
        -j | --java-ontop-args ) shift
                                java_ontop_args=${1%/}
                                ;;
        -l | --test-genes )     shift
                                test_genes=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

for required in ontop_dir_path obda_file_path output_dir_path ontop_property_file ontology_file; do
    if [ -z "${!required}" ]; then
        echo "Missing argument: $required"
        usage
        exit 1
    fi
done
if [ -n "$test_genes" ]; then
    read -r test_low test_high <<< "$test_genes"
    if ! [[ "$test_low" =~ ^[0-9]+$ && "$test_high" =~ ^[0-9]+$ && "$slice_count" =~ ^[1-9][0-9]*$ ]] ||
            [ "$test_low" -ge "$test_high" ]; then
        echo "Invalid test mode, expected \"low high\" bgeeGeneId with low < high, and a number of slices: $test_genes, $slice_count"
        exit 1
    fi
    # Slices of the same number of genes: the window holds only a few genes
    step=$(( (test_high - test_low + slice_count - 1) / slice_count ))
    gene_bounds=$test_low
    for (( first = test_low + step; first < test_high; first += step )); do
        gene_bounds+=",$first"
    done
elif [ -z "$gene_bounds" ]; then
    echo "Missing argument: gene_bounds"
    usage
    exit 1
fi
if [ -n "$xml_catalog" ] && [ ! -f "$xml_catalog" ]; then
    echo "XML catalog not found: $xml_catalog"
    exit 1
fi

here=$(cd "$(dirname "$0")" && pwd)
work_dir="$output_dir_path/parts"

if [ -z "$test_genes" ] && compgen -G "$output_dir_path/*.ttl" > /dev/null; then
    echo "$output_dir_path already contains TTL files, remove them or use another directory"
    exit 1
fi
if [ -e "$work_dir" ]; then
    echo "$work_dir already exists, remove it (it holds the logs of a previous run)"
    exit 1
fi
mkdir -p "$work_dir" || exit 1

#ontop edit java max memory, once for all the processes started afterwards
if [ "$java_ontop_args" != "" ]
then
    java_ontop_args=`echo $java_ontop_args | tr -d \"`
    sed -i  's,ONTOP_JAVA_ARGS=.*.,ONTOP_JAVA_ARGS="'$java_ontop_args'",g' $ontop_dir_path/ontop
fi

echo "Splitting the mapping file..."
slicer_args=(-mapping "$obda_file_path" -prefix "$work_dir/mapping" -bounds "$gene_bounds")
if [ -n "$test_genes" ]; then
    slicer_args+=(-test_window "$test_low $test_high")
fi
# The working directory was created just above and holds nothing else yet
mapping_files=$(perl "$here/slice_mapping.pl" "${slicer_args[@]}") || { echo "The mapping file could not be split"; rm -rf "$work_dir"; exit 1; }

# The ontologies imported are read from the files of the XML catalog rather than downloaded
catalog_args=()
if [ -n "$xml_catalog" ]; then
    catalog_args=(-x "$xml_catalog")
fi
echo "Creating RDF data and serializing them as TURTLE files, in parallel..."
pids=()
names=()
for mapping_file in $mapping_files; do
    name=${mapping_file##*/}
    name=${name#mapping.}
    name=${name%.obda}
    if [ ${#pids[@]} -gt 0 ] && [ "$delay" -gt 0 ]; then
        sleep "$delay"
    fi
    mkdir -p "$work_dir/$name"
    echo "Starting Ontop for the part $name"
    "$ontop_dir_path/ontop" materialize --separate-files -m "$mapping_file" -f turtle -o "$work_dir/$name" -p "$ontop_property_file" -t "$ontology_file" "${catalog_args[@]}" > "$work_dir/$name.log" 2>&1 &
    pids+=($!)
    names+=("$name")
done

failed=0
for i in "${!pids[@]}"; do
    if wait "${pids[$i]}"; then
        echo "Part ${names[$i]} done"
    else
        echo "Part ${names[$i]} FAILED, see $work_dir/${names[$i]}.log"
        failed=1
    fi
done
if [ $failed -ne 0 ]; then
    echo "At least one part failed: the files written are left in $work_dir"
    exit 1
fi

if [ -n "$test_genes" ]; then
    # Prefix declarations are not triples, and the duplicates are not significant in RDF
    for name in "${names[@]}"; do
        if [ "$name" == "ref" ]; then
            cat "$work_dir/ref/"*.ttl | grep -v '^@' | LC_ALL=C sort -u > "$work_dir/ref.sorted"
        else
            cat "$work_dir/$name/"*.ttl
        fi
    done | grep -v '^@' | LC_ALL=C sort -u > "$work_dir/slices.sorted"
    echo "Distinct triples, whole mapping: $(wc -l < "$work_dir/ref.sorted" | tr -d ' ')"
    echo "Distinct triples, union of the slices: $(wc -l < "$work_dir/slices.sorted" | tr -d ' ')"
    if diff "$work_dir/ref.sorted" "$work_dir/slices.sorted" > "$work_dir/ref_vs_slices.diff"; then
        echo "TEST OK: the slices give the same triples as the whole mapping"
        exit 0
    fi
    echo "TEST FAILED: the slices differ from the whole mapping, see $work_dir/ref_vs_slices.diff"
    head -20 "$work_dir/ref_vs_slices.diff"
    exit 1
fi

# The files of all the parts get the name of their part as prefix, as they all start at _0
for name in "${names[@]}"; do
    for ttl_file in "$work_dir/$name/"*.ttl; do
        [ -e "$ttl_file" ] || continue
        mv "$ttl_file" "$output_dir_path/${name}_${ttl_file##*/}"
    done
    rmdir "$work_dir/$name" 2>/dev/null
done

echo "The RDF data were created and saved as TURTLE files..."
