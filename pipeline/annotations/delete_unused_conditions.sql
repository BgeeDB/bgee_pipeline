-- We create the tables in such a way to avoid system lock,
-- see https://dba.stackexchange.com/a/15479

-- First, we identify all conditions used in the database.
-- We don't use a TEMP TABLE because they cannot be used twice in a same query
CREATE TABLE tempSafeToDropCondUsed (PRIMARY KEY(conditionId))
SELECT t1.conditionId FROM cond AS t1 LIMIT 0;
ALTER TABLE tempSafeToDropCondUsed ENGINE=InnoDB;

INSERT INTO tempSafeToDropCondUsed
SELECT DISTINCT t1.conditionId FROM cond AS t1
WHERE EXISTS (SELECT 1 FROM estLibrary WHERE estLibrary.conditionId = t1.conditionId)
OR EXISTS (SELECT 1 FROM affymetrixChip WHERE affymetrixChip.conditionId = t1.conditionId)
OR EXISTS (SELECT 1 FROM inSituSpot WHERE inSituSpot.conditionId = t1.conditionId)
OR EXISTS (SELECT 1 FROM rnaSeqLibraryAnnotatedSample WHERE rnaSeqLibraryAnnotatedSample.conditionId = t1.conditionId)
OR EXISTS (SELECT 1 FROM expression WHERE expression.conditionId = t1.conditionId)
OR EXISTS (SELECT 1 FROM globalCondToCond WHERE globalCondToCond.conditionId = t1.conditionId);

-- Then, we delete the conditions unused.
-- But because of how the query is built, we can't directly delete the conditions
-- (we get an error 'Can't specify target table for update in FROM clause'),
-- so we first store the conditions to delete in another table.
-- Select a condition for deletion if:
CREATE TEMPORARY TABLE condToDelete (PRIMARY KEY(conditionId))
SELECT t1.conditionId FROM cond AS t1 LIMIT 0;
ALTER TABLE condToDelete ENGINE=InnoDB;

INSERT INTO condToDelete
SELECT DISTINCT t1.conditionId FROM cond AS t1
-- 1) the condition is itself not used, and;
WHERE t1.conditionId NOT IN (SELECT conditionId from tempSafeToDropCondUsed)
-- 2) the condition is not the expression table mapping for a condition that is used
-- (Otherwise we would delete the other used condition, because of an 'ON DELETE CASCADE' clause
-- on the `exprMappedConditionId` field)
AND NOT EXISTS (SELECT 1 FROM cond AS t2 WHERE t2.exprMappedConditionId = t1.conditionId
    AND t2.conditionId IN (SELECT conditionId from tempSafeToDropCondUsed));

-- Delete
DELETE t1 FROM cond AS t1 INNER JOIN condToDelete AS t2 ON t1.conditionId = t2.conditionId;

-- We're done, drop the temp tables
DROP TABLE tempSafeToDropCondUsed;
DROP TABLE condToDelete;