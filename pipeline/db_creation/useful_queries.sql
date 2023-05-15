-- To retrieve list of anatomical entities
-- with their taxon constraints on the same line,
-- including all species:
SELECT t1.anatEntityId, t1.anatEntityName,
  GROUP_CONCAT(DISTINCT t3.speciesId ORDER BY t3.speciesId),
  GROUP_CONCAT(DISTINCT CONCAT(t3.genus, ' ', t3.species) ORDER BY t3.speciesId)
FROM anatEntity AS t1
INNER JOIN anatEntityTaxonConstraint AS t2
  ON t1.anatEntityId = t2.anatEntityId
INNER JOIN species AS t3
  ON IF(t2.speciesId IS NULL, 1=1, t2.speciesId = t3.speciesId)
GROUP BY t1.anatEntityId
ORDER BY t1.anatEntityId;

-- To retrieve list of anatomical entities that are part of the cell type graph:
SELECT DISTINCT t1.anatEntityId, t1.anatEntityName
FROM anatEntity AS t1
INNER JOIN anatEntityRelation AS t2 ON t1.anatEntityId = t2.anatEntitySourceId
-- 'GO:0005575' is the ID of the root of the cell type graph, as of Bgee 15.0.
-- Change it depending on the version used.
WHERE t2.anatEntityTargetId = 'GO:0005575'
  AND t2.relationType = 'is_a part_of';
