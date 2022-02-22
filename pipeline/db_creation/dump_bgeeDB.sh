#!/bin/sh

#NOTE Run it from a mounted DCSR NAS folder to have enough space, and be already available everywhere (cluster, ...)

while getopts "d:h:" Option
do
    case $Option in
        d     ) DB_NAME=$OPTARG;;
        h     ) DB_HOST=$OPTARG;;
        *     ) echo "\n\tUnimplemented option chosen\n"; exit 1;;   # Default.
    esac
done

##### Check db name passed as argument
if [ ${#DB_NAME} -lt 1 ]; then
    echo -e "\n\tMissing db name as argument\n\te.g. $0 -d bgee_v15_0 -h rbioinfo.unil.ch\n"
    exit 2
fi
##### Check db host passed as argument
if [ ${#DB_HOST} -lt 1 ]; then
    echo -e "\n\tMissing db host as argument\n\te.g. $0 -d bgee_v15_0 -h rbioinfo.unil.ch\n"
    exit 2
fi

DIR=$(dirname $([ -L $0 ] && readlink -f $0 || echo $0))
LOGIN=`grep '^DBUSER '  $DIR/../Makefile.Config | awk '{print $3}'`
PASSWD=`grep '^DBPASS ' $DIR/../Makefile.Config | awk '{print $3}'`



# Dump all Bgee DB tables  BUT  the huge globalExpression table
time mysqldump -u $LOGIN -p$PASSWD -h $DB_HOST $DB_NAME --ignore-table=$DB_NAME.globalExpression --ignore-table=$DB_NAME.gene > $DB_NAME-allBut-globalExpression.sql

# Dump the globalExpression table by chunk
## Human
time mysqldump -u $LOGIN -p$PASSWD -h $DB_HOST $DB_NAME globalExpression gene --where='(bgeeGeneId IN (SELECT bgeeGeneId FROM gene WHERE speciesId=9606))'  --skip-add-drop-table > $DB_NAME-globalExpression-Human.sql
## Mouse
time mysqldump -u $LOGIN -p$PASSWD -h $DB_HOST $DB_NAME globalExpression gene --where='(bgeeGeneId IN (SELECT bgeeGeneId FROM gene WHERE speciesId=10090))' --skip-add-drop-table > $DB_NAME-globalExpression-Mouse.sql
## Other species
time mysqldump -u $LOGIN -p$PASSWD -h $DB_HOST $DB_NAME globalExpression gene --where='((bgeeGeneId NOT IN (SELECT bgeeGeneId FROM gene WHERE speciesId=9606)) AND (bgeeGeneId NOT IN (SELECT bgeeGeneId FROM gene WHERE speciesId=10090)))' --skip-add-drop-table > $DB_NAME-globalExpression-otherSpecies.sql

#NOTE because of --skip-add-drop-table think to DROP globalExpression and gene tables if not for a fresh db!

exit 0

