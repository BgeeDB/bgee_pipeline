-- this file contains the indexes that do not add any constraints, defined solely
-- for performance issues (unique indexes are therefore not present in this file,
-- but in bgeeConstraint.sql)

-- index needed to improve performances when inserting ranks in affymetrixProbeset table
ALTER TABLE affymetrixProbeset ADD INDEX (bgeeGeneId, bgeeAffymetrixChipId);
