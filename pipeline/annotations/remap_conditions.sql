-- SQL script to remap conditions and expression calls
-- after the table remapCond has been filled (see script remap_conditions.pl).
--
-- First, for expression table, there are two cases:
-- * either the call bgeeGeneId - conditionId using the remapped condition does not exist;
--   in that case we simply have to update the call using the old condtionId with the remapped conditionId,
--   and keep the same expression Id. No update of other tables necessary.
-- * or the call bgeeGeneId - conditionId using the remapped condition already exists;
--   in that case we need to have a mapping between the old expression ID using the old condition,
--   and the remapped expressionId using the remapped condition, and update
--   all the necessary tables accordingly.
--
-- First case:
UPDATE expression
INNER JOIN (
    SELECT DISTINCT t4.expressionId, t3.exprMappedConditionId
    FROM remapCond AS t1
    INNER JOIN cond AS t2 ON t1.incorrectConditionId = t2.conditionId
    INNER JOIN cond AS t3 ON t1.remappedConditionId = t3.conditionId
        AND t2.exprMappedConditionId != t3.exprMappedConditionId
    INNER JOIN expression AS t4 ON t2.exprMappedConditionId = t4.conditionId
    LEFT OUTER JOIN expression AS t5
        ON t4.bgeeGeneId = t5.bgeeGeneId AND t3.exprMappedConditionId = t5.conditionId
    WHERE t5.expressionId IS NULL
) AS expressionToUpdate
ON expression.expressionId = expressionToUpdate.expressionId
SET expression.conditionId = expressionToUpdate.exprMappedConditionId;

-- Second case:
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

UPDATE rnaSeqResult AS t1
INNER JOIN remapExpression AS t2 ON t1.expressionId = t2.incorrectExpressionId
SET t1.expressionId = t2.remappedExpressionId;

UPDATE scRnaSeqFullLengthResult AS t1
INNER JOIN remapExpression AS t2 ON t1.expressionId = t2.incorrectExpressionId
SET t1.expressionId = t2.remappedExpressionId;

UPDATE scRnaSeqTargetBasedResult AS t1
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

UPDATE rnaSeqLibrary AS t1
INNER JOIN remapCond AS t2 ON t1.conditionId = t2.incorrectConditionId
SET t1.conditionId = t2.remappedConditionId;

UPDATE scRnaSeqFullLengthLibrary AS t1
INNER JOIN remapCond AS t2 ON t1.conditionId = t2.incorrectConditionId
SET t1.conditionId = t2.remappedConditionId;

UPDATE scRnaSeqTargetBasedLibraryCellPopulation AS t1
INNER JOIN remapCond AS t2 ON t1.conditionId = t2.incorrectConditionId
SET t1.conditionId = t2.remappedConditionId;

UPDATE deaSampleGroup AS t1
INNER JOIN remapCond AS t2 ON t1.conditionId = t2.incorrectConditionId
SET t1.conditionId = t2.remappedConditionId;

-- Other scripts are used to delete the remaining unused conditions and expressions