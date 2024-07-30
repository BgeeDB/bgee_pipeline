-- We create the tables in such a way to avoid system lock,
-- see https://dba.stackexchange.com/a/15479

-- First, we identify all expressions used in the database.
-- We don't use a TEMP TABLE because they cannot be used twice in a same query
CREATE TABLE tempSafeToDropExpressionUsed (PRIMARY KEY(expressionId))
SELECT t1.expressionId FROM expression AS t1 LIMIT 0;
ALTER TABLE tempSafeToDropExpressionUsed ENGINE=InnoDB;

INSERT INTO tempSafeToDropExpressionUsed
SELECT DISTINCT t1.expressionId FROM expression AS t1
WHERE EXISTS (SELECT 1 FROM expressedSequenceTag WHERE expressedSequenceTag.expressionId = t1.expressionId)
OR EXISTS (SELECT 1 FROM affymetrixProbeset WHERE affymetrixProbeset.expressionId = t1.expressionId)
OR EXISTS (SELECT 1 FROM inSituSpot WHERE inSituSpot.expressionId = t1.expressionId)
OR EXISTS (SELECT 1 FROM rnaSeqLibraryAnnotatedSampleGeneResult WHERE rnaSeqLibraryAnnotatedSampleGeneResult.expressionId = t1.expressionId);

-- Then, we delete the conditions unused.
-- But because of how the query is built, we can't directly delete the conditions
-- (we get an error 'Can't specify target table for update in FROM clause'),
-- so we first store the conditions to delete in another table.
CREATE TEMPORARY TABLE expressionToDelete (PRIMARY KEY(expressionId))
SELECT t1.expressionId FROM expression AS t1 LIMIT 0;
ALTER TABLE expressionToDelete ENGINE=InnoDB;

INSERT INTO expressionToDelete
SELECT DISTINCT t1.expressionId FROM expression AS t1
WHERE t1.expressionId NOT IN (SELECT expressionId from tempSafeToDropExpressionUsed);

-- Delete
DELETE t1 FROM expression AS t1 INNER JOIN expressionToDelete AS t2 ON t1.expressionId = t2.expressionId;

-- We're done, drop the temp tables
DROP TABLE tempSafeToDropExpressionUsed;
DROP TABLE expressionToDelete;
