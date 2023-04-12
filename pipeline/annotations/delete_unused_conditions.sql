-- We create the tables in such a way to avoid system lock,
-- see https://dba.stackexchange.com/a/15479

-- First, we identify all conditions used in the database.
-- We don't use a TEMP TABLE because they cannot be used twice in a same query
CREATE TABLE tempSafeToDropCondUsed (PRIMARY KEY(conditionId))
SELECT t1.conditionId FROM cond AS t1 LIMIT 0;
ALTER TABLE tempSafeToDropCondUsed ENGINE=InnoDB;

INSERT INTO tempSafeToDropCondUsed
SELECT DISTINCT t1.conditionId FROM cond AS t1
LEFT OUTER JOIN estLibrary AS t2 ON t1.conditionId = t2.conditionId
LEFT OUTER JOIN affymetrixChip AS t3 ON t1.conditionId = t3.conditionId
LEFT OUTER JOIN inSituSpot AS t4 ON t1.conditionId = t4.conditionId
LEFT OUTER JOIN rnaSeqLibraryAnnotatedSample AS t5 ON t1.conditionId = t5.conditionId
LEFT OUTER JOIN expression AS t6 ON t1.conditionId = t6.conditionId
LEFT OUTER JOIN globalCondToCond AS t7 ON t1.conditionId = t7.conditionId
WHERE t2.conditionId IS NOT NULL OR t3.conditionId IS NOT NULL OR t4.conditionId IS NOT NULL
OR t5.conditionId IS NOT NULL OR t6.conditionId IS NOT NULL OR t7.conditionId IS NOT NULL;

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