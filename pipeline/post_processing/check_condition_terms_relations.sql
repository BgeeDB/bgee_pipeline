-- Checks the is_a part_of relations of the anatomical entities and cell types used in the conditions
-- of the species @speciesId, against the chains of direct relations valid in that species.
-- The propagation of the conditions, and the API when it retrieves the descendants of a term,
-- read the stored relations (direct and indirect) without walking the graph: a relation missing
-- from the closure, or an indirect relation that no chain supports, propagates incorrectly.
-- Each row is a problem:
--   * 'missing': a chain of direct relations reaches the ancestor, no relation is stored;
--   * 'unsupported': a relation is stored, no chain of direct relations reaches the ancestor;
--   * 'cycle': the term is its own ancestor;
--   * 'not connected to the root': no chain of direct relations reaches BGEE:0000000.
-- A relation, or a term, is valid in the species according to its taxon constraints;
-- a relation also requires its source and its target to be valid in the species.
-- Usage: (echo "SET @speciesId = 9606;"; cat check_condition_terms_relations.sql) | mysql ... bgee_vXX
WITH RECURSIVE
terms AS (
    SELECT anatEntityId AS id FROM cond WHERE speciesId = @speciesId
    UNION
    SELECT cellTypeId FROM cond WHERE speciesId = @speciesId AND cellTypeId IS NOT NULL
),
directRel AS (
    SELECT r.anatEntitySourceId AS sourceId, r.anatEntityTargetId AS targetId
    FROM anatEntityRelation AS r
    WHERE r.relationType = 'is_a part_of' AND r.relationStatus = 'direct'
        AND r.anatEntitySourceId != r.anatEntityTargetId
        AND EXISTS (SELECT 1 FROM anatEntityRelationTaxonConstraint AS rtc
                    WHERE rtc.anatEntityRelationId = r.anatEntityRelationId
                    AND (rtc.speciesId IS NULL OR rtc.speciesId = @speciesId))
        AND EXISTS (SELECT 1 FROM anatEntityTaxonConstraint AS stc
                    WHERE stc.anatEntityId = r.anatEntitySourceId
                    AND (stc.speciesId IS NULL OR stc.speciesId = @speciesId))
        AND EXISTS (SELECT 1 FROM anatEntityTaxonConstraint AS ttc
                    WHERE ttc.anatEntityId = r.anatEntityTargetId
                    AND (ttc.speciesId IS NULL OR ttc.speciesId = @speciesId))
),
storedRel AS (
    SELECT DISTINCT r.anatEntitySourceId AS sourceId, r.anatEntityTargetId AS targetId
    FROM anatEntityRelation AS r
    WHERE r.anatEntitySourceId IN (SELECT id FROM terms)
        AND r.relationType = 'is_a part_of' AND r.relationStatus IN ('direct', 'indirect')
        AND r.anatEntitySourceId != r.anatEntityTargetId
        AND EXISTS (SELECT 1 FROM anatEntityRelationTaxonConstraint AS rtc
                    WHERE rtc.anatEntityRelationId = r.anatEntityRelationId
                    AND (rtc.speciesId IS NULL OR rtc.speciesId = @speciesId))
        AND EXISTS (SELECT 1 FROM anatEntityTaxonConstraint AS stc
                    WHERE stc.anatEntityId = r.anatEntitySourceId
                    AND (stc.speciesId IS NULL OR stc.speciesId = @speciesId))
        AND EXISTS (SELECT 1 FROM anatEntityTaxonConstraint AS ttc
                    WHERE ttc.anatEntityId = r.anatEntityTargetId
                    AND (ttc.speciesId IS NULL OR ttc.speciesId = @speciesId))
),
chain (termId, ancestorId) AS (
    SELECT d.sourceId, d.targetId FROM directRel AS d WHERE d.sourceId IN (SELECT id FROM terms)
    UNION
    SELECT c.termId, d.targetId FROM chain AS c INNER JOIN directRel AS d ON d.sourceId = c.ancestorId
)
SELECT @speciesId AS speciesId, 'missing' AS problem, c.termId, c.ancestorId
FROM chain AS c LEFT OUTER JOIN storedRel AS s ON s.sourceId = c.termId AND s.targetId = c.ancestorId
WHERE s.sourceId IS NULL AND c.termId != c.ancestorId
UNION ALL
SELECT @speciesId, 'unsupported', s.sourceId, s.targetId
FROM storedRel AS s LEFT OUTER JOIN chain AS c ON c.termId = s.sourceId AND c.ancestorId = s.targetId
WHERE c.termId IS NULL
UNION ALL
SELECT @speciesId, 'cycle', c.termId, c.ancestorId
FROM chain AS c WHERE c.termId = c.ancestorId
UNION ALL
SELECT @speciesId, 'not connected to the root', t.id, NULL
FROM terms AS t
WHERE t.id != 'BGEE:0000000'
    AND NOT EXISTS (SELECT 1 FROM chain AS c WHERE c.termId = t.id AND c.ancestorId = 'BGEE:0000000')
ORDER BY problem, termId, ancestorId;
