#!/bin/bash

# Get absolute path for $0, even if symlink, to reach Makefile.Config
#NOTE On Mac OS, readlink is different and does not work like this!
#     Install a kind of Linux readlink from fink or homebrew (in coreutils package)
#     amd put it in the PATH
DIR=$(dirname $([ -L $0 ] && readlink -f $0 || echo $0))
PASSWD=`grep '^DBPASS ' $DIR/../Makefile.Config | awk '{print $3}'`
LOGIN=`grep '^DBUSER '  $DIR/../Makefile.Config | awk '{print $3}'`
DBNAME=bgee_v15

set -e

  echo 'creating the database ... '
  mysql -u $LOGIN -p$PASSWD -e "CREATE DATABASE $DBNAME"
  mysql -u $LOGIN -p$PASSWD $DBNAME < bgeeSchema.sql
  echo 'done'

  echo 'filling tables with data ... '
  mysql -u $LOGIN -p$PASSWD $DBNAME < /var/bgee/dump_anatomy.sql
  mysql -u $LOGIN -p$PASSWD $DBNAME < /var/bgee/dump_all.sql
  echo 'done'

  echo 'adding indexes ... '
  mysql -u $LOGIN -p$PASSWD $DBNAME < bgeeIndex.sql
  echo 'done'

  echo 'adding constraints ... '
  mysql -u $LOGIN -p$PASSWD $DBNAME < bgeeConstraint.sql
  echo 'done'

  echo 'adding foreign keys ... '
  mysql -u $LOGIN -p$PASSWD $DBNAME < bgeeForeignKey.sql
  echo 'done'

exit 0
