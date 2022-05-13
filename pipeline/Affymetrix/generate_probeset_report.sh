#!/bin/sh

# Get absolute path for $0, even if symlink, to reach Makefile.Config
#NOTE On Mac OS, readlink is different and does not work like this!
#     Install a kind of Linux readlink from fink or homebrew (in coreutils package)
#     amd put it in the PATH
DIR=$(dirname $([ -L $0 ] && readlink -f $0 || echo $0))
PASSWD=`grep '^DBPASS ' $DIR/../Makefile.Config | awk '{print $3}'`
LOGIN=`grep '^DBUSER '  $DIR/../Makefile.Config | awk '{print $3}'`
DBNAME=bgee_v15

MYSQL="mysql -u $LOGIN -p$PASSWD $DBNAME"

for chipTypeId in `$MYSQL -e 'SELECT DISTINCT chipTypeId FROM affymetrixChip' | grep -v 'chipTypeId'`; do
    fileName=${chipTypeId}.tsv
    echo "Probeset ID" > $fileName
    $MYSQL -e "SELECT affymetrixProbesetId FROM affymetrixProbeset WHERE bgeeAffymetrixChipId = (SELECT bgeeAffymetrixChipId FROM affymetrixChip WHERE chipTypeId = '${chipTypeId}' LIMIT 1) ORDER BY affymetrixProbesetId" | grep -v "affymetrixProbesetId" >> $fileName

    for bgeeAffymetrixChipId in `$MYSQL -e "SELECT bgeeAffymetrixChipId FROM affymetrixChip WHERE chipTypeId = '$chipTypeId'" | grep -v 'bgeeAffymetrixChipId'`; do
        tempFileName=whatever
        $MYSQL -e "SELECT CONCAT(affymetrixChipId, ' (', bgeeAffymetrixChipId, ')') AS customId FROM affymetrixChip WHERE bgeeAffymetrixChipId = $bgeeAffymetrixChipId" | grep -v "customId" > $tempFileName
        $MYSQL -e "SELECT detectionFlag FROM affymetrixProbeset WHERE bgeeAffymetrixChipId = $bgeeAffymetrixChipId ORDER BY affymetrixProbesetId" | grep -v "detectionFlag" >> $tempFileName
        tempFileName2=whatever2
        paste $fileName $tempFileName > $tempFileName2; mv -f $tempFileName2 $fileName
        rm -f $tempFileName; rm -f $tempFileName2
    done
done

