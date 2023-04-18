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
    echo -e "\n\tMissing db name as argument\n\te.g. $0 -d bgee_v15_1 -h rbioinfo.unil.ch\n"
    exit 2
fi
##### Check db host passed as argument
if [ ${#DB_HOST} -lt 1 ]; then
    echo -e "\n\tMissing db host as argument\n\te.g. $0 -d bgee_v15_1 -h rbioinfo.unil.ch\n"
    exit 2
fi

DIR=$(dirname $([ -L $0 ] && readlink -f $0 || echo $0))
LOGIN=`grep '^DBUSER '  $DIR/../Makefile.Config | awk '{print $3}'`
PASSWD=`grep '^DBPASS ' $DIR/../Makefile.Config | awk '{print $3}'`


MYSQL_OPTIONS='--skip-triggers --no-create-info --no-tablespaces --compact --skip-add-drop-table'

#TODO try mysqlpump instead of mysqldump? The --where option is the key point!
#     https://dev.mysql.com/doc/refman/8.0/en/mysqlpump.html
#     --no-tablespaces        ->  mysqlpump does not dump InnoDB CREATE TABLESPACE statements
#     --ignore-table=         ->  --exclude-tables=
#     --skip-triggers=        ->  --exclude-triggers=
#     --no-create-info=       ->  --no-create-info=
#     --skip-add-drop-table=  ->  
#     --compact=              ->  
#     --where=                ->  
#
#     Maybe now usefull options:
#       --defer-table-indexes
#       --watch-progress


# Dump all Bgee DB tables  BUT  the huge globalExpression table
time mysqldump -u $LOGIN -p$PASSWD -h $DB_HOST $DB_NAME --ignore-table=$DB_NAME.globalExpression $MYSQL_OPTIONS > $DB_NAME-allBut-globalExpression.sql

# Dump the globalExpression table by chunk
## Human
time mysqldump -u $LOGIN -p$PASSWD -h $DB_HOST $DB_NAME globalExpression gene --where='(bgeeGeneId IN (SELECT bgeeGeneId FROM gene WHERE speciesId=9606))'  $MYSQL_OPTIONS > $DB_NAME-globalExpression-Human.sql
## Mouse
time mysqldump -u $LOGIN -p$PASSWD -h $DB_HOST $DB_NAME globalExpression gene --where='(bgeeGeneId IN (SELECT bgeeGeneId FROM gene WHERE speciesId=10090))' $MYSQL_OPTIONS > $DB_NAME-globalExpression-Mouse.sql
## D. melanogaster
time mysqldump -u $LOGIN -p$PASSWD -h $DB_HOST $DB_NAME globalExpression gene --where='(bgeeGeneId IN (SELECT bgeeGeneId FROM gene WHERE speciesId=7227))'  $MYSQL_OPTIONS > $DB_NAME-globalExpression-Dmelanogaster.sql
## Other species
time mysqldump -u $LOGIN -p$PASSWD -h $DB_HOST $DB_NAME globalExpression gene --where='(bgeeGeneId IN (SELECT bgeeGeneId FROM gene WHERE speciesId!=9606 AND speciesId!=10090 AND speciesId!=7227))' $MYSQL_OPTIONS > $DB_NAME-globalExpression-otherSpecies.sql

#NOTE because of --skip-add-drop-table think to DROP tables if not for a fresh db!
#NOTE $DB_NAME-allBut-globalExpression.sql should be loaded last to deal with gene table duplication (and issue with gene table lock with chunks)
# is it still with  --single-transaction  and data only?

exit 0

