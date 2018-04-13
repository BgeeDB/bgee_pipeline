#!/bin/bash

DATA=/var/tmp

# Get absolute path for $0, even if symlink, to reach Makefile.Config
DIR=$(dirname $([ -L $0 ] && readlink -f $0 || echo $0))
PASSWD=`grep '^DBPASS ' $DIR/../Makefile.Config | awk '{print $3}'`
LOGIN=`grep '^DBUSER '  $DIR/../Makefile.Config | awk '{print $3}'`
HOST='localhost'
PORT=3306


# Get & initialize parameters/arguments
OLD_CLEANING=0
ALL_DB=0
while getopts "d:ras:" Option
do
    case $Option in
        d     ) DB=$OPTARG;;
        r     ) OLD_CLEANING=1;;
        a     ) ALL_DB=1;;
        s     ) SPECIES=$OPTARG;;
        *     ) echo "\n\tUnimplemented option chosen\n"; exit 1;;   # Default.
    esac
done


REL=${DB#ensembl_compara_}
TAX=${REL%%_*}
TAX=${TAX:-nothing}
if [ $TAX != "metazoa" ] && [ $TAX != "bacteria" ] && [ $TAX != "fungi" ] && [ $TAX != "plants" ] && [ $TAX != "protists" ]; then
    DB=""
fi

FTP="ftp://ftp.ensemblgenomes.org/pub/$TAX/current/mysql/"
DIRS=4



##### Check db name passed as argument
if [ ${#DB} -lt 1 ]; then
    echo -e "\n\tMissing db name as argument\n\te.g. $0 -d ensembl_compara_metazoa_30_83 [-r -a]\n\tOptions:\n\t-r    Remove the previous database release\n\t-a    Install all ensembl core db/tables, not only gene tree requirements\n\t-s    Species list to only include (species names with '_', separated by ',' e.g. drosophila_ananassae,drosophila_yakuba)\n"
    exit 2
fi

CLEAN=0
if [ $OLD_CLEANING -eq 1 ]; then
    let CLEAN=${DB##*_}-1
fi
# All species if species list not defined
if [ -s $SPECIES ]; then
    SPECIES="*"
fi


##### Go to repository dir
if [ ! -d $DATA ]; then
    echo "\n\tRepository directory for MySQL dumps doesn't seem to exist\n"
    exit 1
else if [ ! -w $DATA ]; then
    echo "\n\tRepository directory for MySQL dumps doesn't seem to be writable\n"
    exit 1
    fi
fi
rm -Rf $DATA/$DB
mkdir $DATA/$DB
cd   $DATA/$DB


##### Remove the previous database, if requested
if [ $CLEAN -gt 0 ]; then
    echo "\t__ Removing the previous MySQL databases __"
    for old_db in `find /var/lib/mysql/ -maxdepth 1 -type d -name \*[0-9]_$CLEAN\* -exec basename {} \;`
    do
        mysql -u $LOGIN -p$PASSWD -h $HOST -P $PORT -e "drop database if exists $old_db"
    done
    mysql -u $LOGIN -p$PASSWD -h $HOST -P $PORT -e "flush tables"
fi

# Everything in lowercase
SPECIES=`echo $SPECIES | tr '[:upper:]' '[:lower:]'`
array=(${SPECIES//,/ })


##### Download MySQL schemata and required dumps for species core database
echo "\t__ Downloading ensembl_core    current release MySQL schemata and required dumps __"
for species in "${array[@]}"
do
    if [ $ALL_DB -eq 1 ]; then
        echo "\t___ WARNING ! You must have enough disk space, more than 40 GB ___"
        wget -q -r --cut-dirs=$DIRS -nH -A "*" -R "CHECKSUMS*" "$FTP/${species}_core_${REL#${TAX}_}_*"
        wget -q -r --cut-dirs=$DIRS -nH -A "*" -R "CHECKSUMS*" "$FTP/${species}_funcgen_${REL#${TAX}_}_*"
    else
        echo "\t___ WARNING ! You must have enough disk space, more than 3 GB ___"
        wget -q -r --cut-dirs=$DIRS -nH -A "*_core_*.sql.gz, assembly.*, coord_system.*, dna.*, exon*, meta*, seq_region.*, transcript.*, transcript_stable_id.*, translation.*, translation_stable_id.*" -R "prediction_transcript.*, CHECKSUMS*" "$FTP/${species}_core_${REL#${TAX}_}_*"
        find . -mindepth 1 -type d ! -name \*_core_\* -exec rm -Rf {} 2>/dev/null \;
    fi
done
rm -Rf *_mart_* *.sql_40.gz


##### Download MySQL schemata and required dumps for ensembl_website_X database
# NOTE not in ensembl genomes compara !
#echo "\t__ Downloading ensembl_website current release MySQL dumps __"
#wget -q -r --cut-dirs=$DIRS -nH -A "*" -R "CHECKSUMS.*" "$FTP/ensembl_website_${DB##*_}"


##### Download MySQL dumps from current ensembl_compara release
echo "\t__ Downloading ensembl_compara current release MySQL dumps __"
if [ $ALL_DB -eq 1 ]; then
    wget -q -r --cut-dirs=$DIRS -nH -A "*" -R "CHECKSUMS*" "$FTP/$DB/*"
else
    wget -q -r --cut-dirs=$DIRS -nH -A "ensembl_compara_*.sql.gz, genome_db.*, member.*, meta.*, ncbi_taxa_name.*, ncbi_taxa_node.*, protein_tree_attr.*, protein_tree_member.*, protein_tree_member_score.*, protein_tree_node.*, protein_tree_stable_id.*, protein_tree_tag.*, sequence.*, subset_member.*" -R "domain_member.*, family_member.*, CHECKSUMS*" "$FTP/$DB/*"
fi
du -h $DATA/$DB/


##### Remove useless data: hmmer, hmmer_matches, ortholog, ...?
# Required: meta.txt.gz ncbi_taxa_name.txt.gz ncbi_taxa_node.txt.gz
#    protein_tree_node.txt.gz protein_tree_stable_id.txt.gz protein_tree_tag.txt.gz protein_tree_member.txt.gz
#    sequence.* ?
#    member.txt.gz
#    genome_db.txt.gz ?

# In species core db need:
# meta.txt.gz (ACHTUNG! dat files for tables have the same name between species !)
# meta_coord.txt.gz
# transcript.txt.gz
#       dna.txt.gz ?
# transcript_stable_id.txt.gz
# translation_stable_id.txt.gz
# translation.txt.gz
# seq_region.txt.gz
# coord_system.txt.gz
# exon.txt.gz
# exon_transcript.txt.gz
# exon_stable_id.txt.gz
#       gene.txt.gz ?
#       gene_stable_id.txt.gz ?



##### Extract sql file(s)
echo "\t__ Extracting SQL files __"
find . -type f -name \*.sql.gz -exec gunzip {} \;


##### Use Database $DB & then create tables from sql file(s)
# => NEED an account with create tables rights in MySQL
# => NEED to create db in advance, to avoid to use MySQL root account here
# => NEED GRANT CREATE, SELECT, INSERT, UPDATE, INDEX, ALTER ON $DB TO selectomeROOT@localhost; before !
echo "\t__ Creating MySQL database schemata __"
if [ $ALL_DB -eq 1 ]; then
    echo "\t___ WARNING ! You must have enough disk space: about 150 GB ___"
else
    echo "\t___ WARNING ! You must have enough disk space: about 10 GB ___"
fi
for sql in `find . -type f -name \*.sql`
do
    db=${sql/#.*\//}
    BD=${db/%.sql/}
    mysql -u $LOGIN -p$PASSWD -h $HOST -P $PORT -e "drop database if exists $BD"
    mysql -u $LOGIN -p$PASSWD -h $HOST -P $PORT -e "create database $BD character set utf8"
    mysql -u $LOGIN -p$PASSWD -h $HOST -P $PORT $BD < $sql
done
find . -type f -name \*.sql -exec rm -f {} \;


##### Fill DB with proper table content
for dump in `find . -type f -name \*.txt\* | sort`
do
    echo "\t___ Extracting MySQL data: $dump ___"
    gunzip $dump
    # \*.txt.gz         Ensembl compara 55+,         TreeFam 4     syntax
    # \*.txt.table.gz   Ensembl compara metazoa2 54, TreeFam 5/6/7 syntax

    if [ ! -s ${dump%%.gz} ]; then
        echo "\t\t    WARNING !!! this ${dump%%.gz} file is empty !!!"
    fi

    db=${dump/#.\//}
    BD=${db/%\/*/}

    basename=${dump##*/}
    echo "\t___ Loading    MySQL data: ${dump%%.gz} ___" #To see which file is running if problems
    mysql -u $LOGIN -p$PASSWD -h $HOST -P $PORT $BD -e "load data infile '$PWD/${dump%%.gz}' into table ${basename%%.*}"

    rm -f ${dump%%.gz}
done
mysql -u $LOGIN -p$PASSWD -h $HOST -P $PORT -e "FLUSH tables"


##### Clean tmp dir
echo "\t__ Cleaning __"
find . -type f -print #To check is smthg is staying there


##### Back to original dir
cd -
rm -Rf $DATA/$DB/
exit 0

