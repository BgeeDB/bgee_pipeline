-- in that case we only wanted to remap some developmental stages.
-- create a table to define the stage remapping.
CREATE TEMPORARY TABLE remapStage (
    incorrectStageId varchar(20) not null,
    remappedStageId varchar(20) not null,
    speciesId mediumint unsigned not null,
    PRIMARY KEY(incorrectStageId, speciesId, remappedStageId),
    INDEX(remappedStageId, speciesId)
);
INSERT INTO remapStage (incorrectStageId, remappedStageId, speciesId) VALUES
('UBERON:0000113', 'UBERON:0000066', 7955),
('FBdv:00007083', 'DpseDv:0000007', 7237),
('FBdv:00007083', 'DmojDv:0000007', 7230),
('UBERON:0018241', 'UBERON:0000113', 10090),
('UBERON:0018241', 'UBERON:0000113', 10116);


-- now, identify if there are some existings conditions to remap the incorrect conditions to.
CREATE TEMPORARY TABLE remapCond (PRIMARY KEY(incorrectConditionId, remappedConditionId))
SELECT DISTINCT t1.conditionId AS incorrectConditionId, t3.conditionId AS remappedConditionId
FROM cond AS t1
INNER JOIN remapStage AS t2 ON t1.stageId = t2.incorrectStageId AND t1.speciesId = t2.speciesId
LEFT OUTER JOIN cond AS t3 ON t1.speciesId = t3.speciesId AND t1.anatEntityId = t3.anatEntityId
    AND t2.remappedStageId = t3.stageId AND t1.sex = t3.sex AND t1.sexInferred = t3.sexInferred
    AND t1.strain = t3.strain;

-- For incorrect conditions with no existing conditions to remap to,
-- simply update the incorrect field in the condition table.
UPDATE cond AS t1
INNER JOIN remapCond AS t2 ON t2.incorrectConditionId = t1.conditionId
    AND (t2.remappedConditionId IS NULL OR t2.remappedConditionId = 0)
INNER JOIN remapStage AS t3 ON t3.incorrectStageId = t1.stageId AND t3.speciesId = t1.speciesId
SET t1.stageId = t3.remappedStageId;


-- PER DATA TYPE CHANGES: for incorrect conditions with existing conditions to remap to,
-- we need to update each table holding annotations.

-- Conditions that should really be merged with other ones were only affecting RNA-Seq data
-- in that case.
UPDATE rnaSeqLibrary AS t1
INNER JOIN remapCond AS t2 ON t2.incorrectConditionId = t1.conditionId
SET t1.conditionId = t2.remappedConditionId
WHERE (t2.remappedConditionId IS NOT NULL AND t2.remappedConditionId != 0);


-- ********************************************************
-- TO AVOID: this was needed because expression data were computed before the condition corrections.
-- You should normally check conditions before expression insertion.
-- Update of expressionId in rnaSeqResult shoud be done by
-- pipeline/RNA_Seq/3Insertion/insert_rna_seq_expression.pl with a simulatenous insertion of data
-- into rnaSeqExperimentExpression. The following commands would work only if there is no overlap
-- in a same experiment of incorrect condition/remapped conditions, you must make sure of this.

CREATE TEMPORARY TABLE remapExpression (PRIMARY KEY(incorrectExpressionId, remappedExpressionId))
SELECT STRAIGHT_JOIN t2.expressionId AS incorrectExpressionId, t3.expressionId AS remappedExpressionId
FROM remapCond AS t1
INNER JOIN expression AS t2 ON t2.conditionId = t1.incorrectConditionId
INNER JOIN expression AS t3 ON t3.bgeeGeneId = t2.bgeeGeneId AND t3.conditionId = t1.remappedConditionId
WHERE (t1.remappedConditionId IS NOT NULL AND t1.remappedConditionId != 0);

UPDATE rnaSeqResult AS t1
INNER JOIN remapExpression AS t2 ON t2.incorrectExpressionId = t1.expressionId
SET t1.expressionId = t2.remappedExpressionId;

UPDATE rnaSeqExperimentExpression AS t1
INNER JOIN remapExpression AS t2 ON t2.incorrectExpressionId = t1.expressionId
SET t1.expressionId = t2.remappedExpressionId;


-- To do at the end:
DELETE t1
FROM expression AS t1
LEFT OUTER JOIN expressedSequenceTag AS t2 ON t1.expressionId = t2.expressionId
LEFT OUTER JOIN affymetrixProbeset AS t3 ON t1.expressionId = t3.expressionId
LEFT OUTER JOIN inSituSpot AS t4 ON t1.expressionId = t4.expressionId
LEFT OUTER JOIN rnaSeqResult AS t5 ON t1.expressionId = t5.expressionId
WHERE t2.expressionId IS NULL AND t3.expressionId IS NULL AND t4.expressionId IS NULL
    AND t5.expressionId IS NULL;

UPDATE expression AS t1
INNER JOIN remapCond AS t2 ON t2.incorrectConditionId = t1.conditionId
SET t1.conditionId = t2.remappedConditionId
WHERE (t2.remappedConditionId IS NOT NULL AND t2.remappedConditionId != 0);

-- ********************************************************

-- When all done, delete incorrect conditions that were remapped to existing condition.

