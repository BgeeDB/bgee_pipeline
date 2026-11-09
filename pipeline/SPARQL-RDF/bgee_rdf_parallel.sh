#!/bin/bash
usage()
{
    echo "usage: easyBgee to RDF triples with several Ontop processes working in parallel.
                The mapping is split by slice of globalExpressionId (see slice_mapping.pl)
                [ {-x | --ontop} <ontop directory> ]
                [ {-m | --mapping-file} <mapping file> ]
                [ {-o | --output-dir-path} <Turtle file output> ]
                [ {-p | --ontop-property-file} <property file>  ]
                [ {-t | --ontology-file} <ontology file> ]
                [ {-i | --id-range} <\"lowest highest\" globalExpressionId of easyBgee> ]
                [ {-n | --slices} <number of slices, i.e. Ontop processes reading globalExpression> (optional, default 12) ]
                [ {-d | --delay} <seconds between the start of two Ontop processes> (optional, default 30) ]
                [ {-j | --java-ontop-args } <Java Ontop memory arguments > (optional) ]
                [ {-l | --test-id-limit} <globalExpressionId> (optional, test mode) ]
                [ {-h | --help} <usage help> ]

   The range of identifiers is retrieved with:
    SELECT MIN(globalExpressionId), MAX(globalExpressionId) FROM globalExpression;

   Example:
    ./bgee_rdf_parallel.sh -m ./conf/genex_adapt.obda -o ./ttl -p ./conf/genex_adapt.properties -x ./ontop-cli-3 -t ./conf/genex_adapt.owl -i '1 1385000000' -n 12 -j '-Xmx16G'

   Test mode: the identifiers are limited to those lower than the value given, e.g. the lowest identifier + 200000.
   The mapping is materialized in one go and in slices, and the triples obtained are compared.
   Nothing is written to the output directory but the working files.
    ./bgee_rdf_parallel.sh -m ... -o ./ttl_test -i '1 1385000000' -l 200001
    "
}
#Variables to be assigned, currently the default is no value.
ontop_dir_path=
obda_file_path=
output_dir_path=
ontop_property_file=
ontology_file=
id_range=
slice_count=12
delay=30
java_ontop_args=
test_id_limit=

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
        -i | --id-range )       shift
                                id_range=$1
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
        -l | --test-id-limit )  shift
                                test_id_limit=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

for required in ontop_dir_path obda_file_path output_dir_path ontop_property_file ontology_file id_range; do
    if [ -z "${!required}" ]; then
        echo "Missing argument: $required"
        usage
        exit 1
    fi
done
read -r min_id max_id <<< "$id_range"
if [ -z "$min_id" ] || [ -z "$max_id" ]; then
    echo "Invalid id range, expected \"lowest highest\": $id_range"
    exit 1
fi

here=$(cd "$(dirname "$0")" && pwd)
work_dir="$output_dir_path/parts"

if [ -z "$test_id_limit" ] && compgen -G "$output_dir_path/*.ttl" > /dev/null; then
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
slicer_args=(-mapping "$obda_file_path" -prefix "$work_dir/mapping" -slices "$slice_count" -min_id "$min_id" -max_id "$max_id")
if [ -n "$test_id_limit" ]; then
    slicer_args+=(-id_limit "$test_id_limit")
fi
mapping_files=$(perl "$here/slice_mapping.pl" "${slicer_args[@]}") || { echo "The mapping file could not be split"; exit 1; }

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
    "$ontop_dir_path/ontop" materialize --separate-files -m "$mapping_file" -f turtle -o "$work_dir/$name" -p "$ontop_property_file" -t "$ontology_file" > "$work_dir/$name.log" 2>&1 &
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

if [ -n "$test_id_limit" ]; then
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
