#!/bin/sh

# Get absolute path for $0, even if symlink, to reach Makefile.Config
#NOTE On Mac OS, readlink is different and does not work like this!
#     Install a kind of Linux readlink from fink or homebrew (in coreutils package)
#     amd put it in the PATH
DIR=$(dirname $([ -L $0 ] && readlink -f $0 || echo $0))
PASSWD=`grep '^DBPASS ' $DIR/../Makefile.Config | awk '{print $3}'`
LOGIN=`grep '^DBUSER '  $DIR/../Makefile.Config | awk '{print $3}'`
DBNAME=bgee_v14

MYSQL="mysql -u $LOGIN -p$PASSWD $DBNAME"

for chipTypeId in `$MYSQL -e 'select distinct chipTypeId from affymetrixChip' | grep -v 'chipTypeId'`; do
    fileName=${chipTypeId}.tsv
    echo "Probeset ID" > $fileName
    $MYSQL -e "select affymetrixProbesetId from affymetrixProbeset where bgeeAffymetrixChipId = (select bgeeAffymetrixChipId from affymetrixChip where chipTypeId = '${chipTypeId}' limit 1) order by affymetrixProbesetId" | grep -v "affymetrixProbesetId" >> $fileName

    for bgeeAffymetrixChipId in `$MYSQL -e "select bgeeAffymetrixChipId from affymetrixChip where chipTypeId = '$chipTypeId'" | grep -v 'bgeeAffymetrixChipId'`; do
        tempFileName=whatever
        $MYSQL -e "select concat(affymetrixChipId, ' (', bgeeAffymetrixChipId, ')') as customId from affymetrixChip where bgeeAffymetrixChipId = $bgeeAffymetrixChipId" | grep -v "customId" > $tempFileName
        $MYSQL -e "select detectionFlag from affymetrixProbeset where bgeeAffymetrixChipId = $bgeeAffymetrixChipId order by affymetrixProbesetId" | grep -v "detectionFlag" >> $tempFileName
        tempFileName2=whatever2
        paste $fileName $tempFileName > $tempFileName2; mv -f $tempFileName2 $fileName
        rm -f $tempFileName; rm -f $tempFileName2
    done
done
