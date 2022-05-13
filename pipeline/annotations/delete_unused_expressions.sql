-- First, we identify all expressions used in the database.
-- We don't use a TEMP TABLE because they cannot be used twice in a same query
CREATE TABLE tempSafeToDropExpressionUsed (PRIMARY KEY(expressionId))
SELECT DISTINCT t1.expressionId FROM expression AS t1
LEFT OUTER JOIN expressedSequenceTag AS t2 ON t1.expressionId = t2.expressionId
LEFT OUTER JOIN affymetrixProbeset AS t3 ON t1.expressionId = t3.expressionId
LEFT OUTER JOIN inSituSpot AS t4 ON t1.expressionId = t4.expressionId
LEFT OUTER JOIN rnaSeqResult AS t5 ON t1.expressionId = t5.expressionId
LEFT OUTER JOIN scRnaSeqFullLengthResult AS t6 ON t1.expressionId = t6.expressionId
LEFT OUTER JOIN scRnaSeqTargetBasedResult AS t7 ON t1.expressionId = t7.expressionId
WHERE t2.expressionId IS NOT NULL OR t3.expressionId IS NOT NULL OR t4.expressionId IS NOT NULL
OR t5.expressionId IS NOT NULL OR t6.expressionId IS NOT NULL OR t7.expressionId IS NOT NULL;

-- Then, we delete the conditions unused.
-- But because of how the query is built, we can't directly delete the conditions
-- (we get an error 'Can't specify target table for update in FROM clause'),
-- so we first store the conditions to delete in another table.
CREATE TEMPORARY TABLE expressionToDelete (PRIMARY KEY(expressionId))
SELECT DISTINCT t1.expressionId FROM expression AS t1
WHERE t1.expressionId NOT IN (SELECT expressionId from tempSafeToDropExpressionUsed);

-- Delete
DELETE t1 FROM expression AS t1 INNER JOIN expressionToDelete AS t2 ON t1.expressionId = t2.expressionId;

-- We're done, drop the temp tables
DROP TABLE tempSafeToDropExpressionUsed;
DROP TABLE expressionToDelete;
