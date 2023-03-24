-- this file contains the indexes that do not add any constraints, defined solely
-- for performance issues or for FK definitions from other tables
-- (unique indexes are therefore not present in this file, but in bgeeConstraint.sql)

-- index needed to improve performances when inserting ranks in result tables
ALTER TABLE affymetrixProbeset ADD INDEX (bgeeAffymetrixChipId, expressionId, bgeeGeneId, normalizedSignalIntensity);
ALTER TABLE rnaSeqPopulationCaptureSpeciesMaxRank ADD INDEX (rnaSeqLibraryAnnotatedSampleId, expressionId, bgeeGeneId, abundance);

-- index generated to fasten the retrieval of raw data as proposed in the
-- DBA StackExchange issue : https://dba.stackexchange.com/questions/320207/optimization-with-subquery-not-working-as-expected
-- the improvement provided by these index has not been tested
ALTER TABLE affymetrixProbeset ADD INDEX(bgeeAffymetrixChipId, affymetrixProbesetId, bgeeGeneId);
ALTER TABLE affymetrixChip ADD INDEX(conditionId,  bgeeAffymetrixChipId);
ALTER TABLE cond ADD INDEX(speciesId, conditionId);