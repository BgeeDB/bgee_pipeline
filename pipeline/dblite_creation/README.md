Goal: Create the Bgee lite database, and insert data source information.

## Details

This step will create a lite version of bgee following the DB schema defined in the file bgeeLiteSchema.sql
The new database is named `bgeelite_vRELEASE` (e.g., `bgeelite_v14`).

## Data generation

* If it is the first time you execute this step in this pipeline run: `make clean`

* Run Makefile: `make`


