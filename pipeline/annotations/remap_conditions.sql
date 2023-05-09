-- SQL script to remap conditions and expression calls
-- after the table remapCond has been filled (see script remap_conditions.pl).
--
-- First, for expression table, there are two cases:
-- * the call bgeeGeneId - conditionId using the remapped condition does not exist;
--   in that case we could simply update the call using the old condtionId with the remapped conditionId,
--   and keep the same expression Id. But it's not that easy, for instance if both cond1 and cond2
--   are remapped to cond3, if we have gene1-cond1, gene1-cond2, doing just a update would lead
--   to have two lines gene1-cond3, and the update would fail because of the uniqueness of (bgeeGeneId, conditionId)
--   in the expression table. For that reason, rather than just updating the calls, we are going to insert
--   new calls with the remapping. And we will end up dealing with the second situation described just below
--   in all cases.
-- * the call bgeeGeneId - conditionId using the remapped condition already exists;
--   in that case we need to have a mapping between the old expression ID using the old condition,
--   and the remapped expressionId using the remapped condition, and update
--   all the necessary tables accordingly.
--
-- First case, we create new expression rows:
INSERT INTO expression (bgeeGeneId, conditionId)
SELECT bgeeGeneId, exprMappedConditionId FROM (
    SELECT DISTINCT t4.bgeeGeneId, t3.exprMappedConditionId
    FROM remapCond AS t1
    INNER JOIN cond AS t2 ON t1.incorrectConditionId = t2.conditionId
    INNER JOIN cond AS t3 ON t1.remappedConditionId = t3.conditionId
        AND t2.exprMappedConditionId != t3.exprMappedConditionId
    INNER JOIN expression AS t4 ON t2.exprMappedConditionId = t4.conditionId
    LEFT OUTER JOIN expression AS t5
        ON t4.bgeeGeneId = t5.bgeeGeneId AND t3.exprMappedConditionId = t5.conditionId
    WHERE t5.expressionId IS NULL
) AS expressionToInsert;

-- Second case, and also dealing with the new lines inserted by the previous query:
DELETE FROM remapExpression;
INSERT INTO remapExpression (incorrectExpressionId, remappedExpressionId)
SELECT DISTINCT t4.expressionId AS 'incorrectExpressionId', t5.expressionId AS 'remappedExpressionId'
FROM remapCond AS t1
INNER JOIN cond AS t2 ON t1.incorrectConditionId = t2.conditionId
INNER JOIN cond AS t3 ON t1.remappedConditionId = t3.conditionId
    AND t2.exprMappedConditionId != t3.exprMappedConditionId
INNER JOIN expression AS t4 ON t2.exprMappedConditionId = t4.conditionId
INNER JOIN expression AS t5
    ON t4.bgeeGeneId = t5.bgeeGeneId AND t3.exprMappedConditionId = t5.conditionId;

-- Now we update all tables using an expressionId using the table remapExpression
UPDATE expressedSequenceTag AS t1
INNER JOIN remapExpression AS t2 ON t1.expressionId = t2.incorrectExpressionId
SET t1.expressionId = t2.remappedExpressionId;

UPDATE affymetrixProbeset AS t1
INNER JOIN remapExpression AS t2 ON t1.expressionId = t2.incorrectExpressionId
SET t1.expressionId = t2.remappedExpressionId;

UPDATE inSituSpot AS t1
INNER JOIN remapExpression AS t2 ON t1.expressionId = t2.incorrectExpressionId
SET t1.expressionId = t2.remappedExpressionId;

UPDATE rnaSeqLibraryAnnotatedSampleGeneResult AS t1
INNER JOIN remapExpression AS t2 ON t1.expressionId = t2.incorrectExpressionId
SET t1.expressionId = t2.remappedExpressionId;

-- And now we update all tables using a conditionId (except the table `expression` already updated)
-- using the table remapCond
UPDATE estLibrary AS t1
INNER JOIN remapCond AS t2 ON t1.conditionId = t2.incorrectConditionId
SET t1.conditionId = t2.remappedConditionId;

UPDATE affymetrixChip AS t1
INNER JOIN remapCond AS t2 ON t1.conditionId = t2.incorrectConditionId
SET t1.conditionId = t2.remappedConditionId;

UPDATE inSituSpot AS t1
INNER JOIN remapCond AS t2 ON t1.conditionId = t2.incorrectConditionId
SET t1.conditionId = t2.remappedConditionId;

UPDATE rnaSeqLibraryAnnotatedSample AS t1
INNER JOIN remapCond AS t2 ON t1.conditionId = t2.incorrectConditionId
SET t1.conditionId = t2.remappedConditionId;

-- Other scripts are used to delete the remaining unused conditions and expressions
