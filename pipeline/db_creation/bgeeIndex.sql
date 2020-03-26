-- this file contains the indexes that do not add any constraints, defined solely
-- for performance issues or for FK definitions from other tables
-- (unique indexes are therefore not present in this file, but in bgeeConstraint.sql)

-- index needed for FK constraints from table geneToOma.
-- TODO: geneToOma needs to be updated and is currently not filled, to reevaluate
ALTER TABLE OMAHierarchicalGroup ADD INDEX (OMAGroupId);

-- index needed to improve performances when inserting ranks in affymetrixProbeset table
ALTER TABLE affymetrixProbeset ADD INDEX (bgeeGeneId, bgeeAffymetrixChipId);
