**Goal**: create the Bgee database, and insert data source information.

## Details

This step will call in order the following files, after creating a new database:
[bgeeSchema.sql](bgeeSchema.sql), [insert_data_sources.sql](insert_data_sources.sql), [bgeeConstraint.sql](bgeeConstraint.sql), [bgeeIndex.sql](bgeeIndex.sql), [bgeeForeignKey.sql](bgeeForeignKey.sql). The new database is named `bgee_vRELEASE` (e.g., `bgee_v14`).

## Data generation

* If it is the first time you execute this step in this pipeline run: `make clean`

* Run Makefile: `make`

* At the end of the pipeline, update data source information, with version used and types of data pulled into bgee: `make update_data_sources.sql`

## Data verification

* `generated_files/db_creation/step_verification_RELEASE.txt` should contain: the number of tables
in the database, and the number of data sources that were inserted into the `dataSource` table.

## Error handling


## Other notable Makefile targets

* update the data sources with version information and types of data used in Bgee (to be done at the end of the pipeline):
  `make update_data_sources.sql`

* drop the database: `make dropDatabaseBgeeRELEASE` (e.g., `make dropDatabaseBgee14`)

## Notes

To manually install the Bgee database:

1. create the database, e.g.: `mysql -u root -p -e "create database bgee_vXX"`
2. insert the database schema, e.g.: `mysql -u root -p bgee_vXX < bgeeSchema.sql`
3. insert the data, e.g.: `mysql -u root -p bgee_vXX < dump_bgee_vXX.sql`
4. create the constraints, e.g.: `mysql -u root -p bgee_vXX < bgeeConstraint.sql`
5. create the indexes, e.g.: `mysql -u root -p bgee_vXX < bgeeIndex.sql`
6. create the foreign key constraints (that will create required indexes as well), e.g.:
`mysql -u root -p bgee_vXX < bgeeForeignKey.sql`

**Important remark**: altering a table after data insertion, to add indexes and constraints,
is faster than doing it before data insertion,
but can fail if the table is very large, with the error 1206:
"ERROR 1206 (HY000): The total number of locks exceeds the lock table size".
To solve this problem, you have to increase the buffer pool size (see e.g. http://bugs.mysql.com/bug.php?id=9975), or you have to insert the data AFTER indexes and foreign key constraints generation.
You'd rather increase the buffer pool size.

## non Ensembl modifications

* In `insert_data_sources.sql` file. Add the line: 

```
(33, 'nonEnsembl', '', '', '', '', 'Independent genome assemblies and annotations (by Maker and Blast2GO)', 0, '');
```

* In `update_data_sources.sql` file. Add the line: 

```
-- non Ensembl genomes and annotations by Maker + Blast2GO
update dataSource set releaseDate = '2019-04-22', releaseVersion = '',                             displayOrder = 12  where dataSourceId = 33;
```




