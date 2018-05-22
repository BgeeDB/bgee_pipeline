Goal: Create the Bgee database, and insert data source information.

## Details

This step will create a light version of bgee following the DB schema defined in the file bgeeLightSchema.sql 
The new database is named `bgeelight_vRELEASE` (e.g., `bgee_v14`).

## Data generation

* If it is the first time you execute this step in this pipeline run: `make clean`

* Run Makefile: `make`


