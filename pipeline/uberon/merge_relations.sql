-- the aim of this script is to remap relations going through an old term
-- to make them go through a new term
-- WARNING: this operation is risky, it is advised to produce first a dump
-- of the anatEntityRelation and anatEntityRelationTaxonConstraint tables.
-- Examine precisely the relation outgoing from the oldTargetId before modifications
-- to ensure no errors will be introduced
--
-- Afterwards, you will have to reinfer indirect relations by using the Java API

SET @speciesId    = 7227;
SET @oldSourceId  = 'FBbt:00007002';
SET @newSourceId  = 'CL:0000000';
SET @oldTargetId  = 'FBbt:00007002';
SET @newTargetId = 'CL:0000000';

--
-- Identify relations that would be duplicated
--
CREATE TEMPORARY TABLE duplicatedRels
SELECT anatEntityRelationId AS oldRelId, anatEntityRelationId AS newRelId FROM anatEntityRelation LIMIT 0;
ALTER TABLE duplicatedRels ENGINE=InnoDB;

-- we do it in two queries for performance issues
INSERT INTO duplicatedRels
SELECT t1.anatEntityRelationId, t3.anatEntityRelationId
FROM anatEntityRelation AS t1
INNER JOIN anatEntityRelationTaxonConstraint AS t2
ON t1.anatEntityRelationId = t2.anatEntityRelationId
INNER JOIN anatEntityRelation AS t3
  ON t1.anatEntitySourceId = @oldSourceId AND t3.anatEntitySourceId = @newSourceId
  AND t1.anatEntityTargetId = t3.anatEntityTargetId
  AND t1.relationType = t3.relationType AND t1.relationStatus = t3.relationStatus
WHERE t2.speciesId = @speciesId OR t2.speciesId IS NULL;

INSERT INTO duplicatedRels
SELECT t1.anatEntityRelationId, t3.anatEntityRelationId
FROM anatEntityRelation AS t1
INNER JOIN anatEntityRelationTaxonConstraint AS t2
ON t1.anatEntityRelationId = t2.anatEntityRelationId
INNER JOIN anatEntityRelation AS t3
  ON t1.anatEntityTargetId = @oldTargetId AND t3.anatEntityTargetId = @newTargetId
  AND t1.anatEntitySourceId = t3.anatEntitySourceId
  AND t1.relationType = t3.relationType AND t1.relationStatus = t3.relationStatus
WHERE t2.speciesId = @speciesId OR t2.speciesId IS NULL;

--
-- identify relations for which we need to update taxon constraints
--
CREATE TEMPORARY TABLE duplicatedRelsUpdateTaxonConstraint
SELECT anatEntityRelationId FROM anatEntityRelationTaxonConstraint LIMIT 0;
ALTER TABLE duplicatedRelsUpdateTaxonConstraint ENGINE=InnoDB;

INSERT INTO duplicatedRelsUpdateTaxonConstraint
SELECT DISTINCT anatEntityRelationId FROM anatEntityRelationTaxonConstraint AS t1
INNER JOIN duplicatedRels AS t2 ON t1.anatEntityRelationId = t2.newRelId
WHERE t1.speciesId IS NOT NULL
AND NOT EXISTS (
  SELECT 1 FROM anatEntityRelationTaxonConstraint AS t3
  WHERE t3.anatEntityRelationId = t1.anatEntityRelationId AND t3.speciesId = @speciesId
);

-- Insert the taxon constraint for the new species
INSERT INTO anatEntityRelationTaxonConstraint
SELECT anatEntityRelationId, @speciesId FROM duplicatedRelsUpdateTaxonConstraint;

-- Now we replace with null speciesId if relation exists in all species
CREATE TEMPORARY TABLE relsExistInAllSpecies
SELECT anatEntityRelationId FROM anatEntityRelationTaxonConstraint LIMIT 0;
ALTER TABLE relsExistInAllSpecies ENGINE=InnoDB;

INSERT INTO relsExistInAllSpecies
SELECT anatEntityRelationId FROM anatEntityRelationTaxonConstraint
GROUP BY anatEntityRelationId
HAVING COUNT(DISTINCT speciesId) = (SELECT COUNT(*) FROM species);
-- we delete the taxon constraints and replace them
DELETE t1 FROM anatEntityRelationTaxonConstraint AS t1
INNER JOIN relsExistInAllSpecies AS t2
ON t1.anatEntityRelationId = t2.anatEntityRelationId;
INSERT INTO anatEntityRelationTaxonConstraint (anatEntityRelationId, speciesId)
SELECT anatEntityRelationId, NULL FROM relsExistInAllSpecies;

--
-- Now we remove the relations that would be duplicated once updated
--
DELETE t1 FROM anatEntityRelationTaxonConstraint AS t1
INNER JOIN duplicatedRels AS t2
ON t1.anatEntityRelationId = t2.oldRelId;

DELETE t1 FROM anatEntityRelation AS t1
INNER JOIN duplicatedRels AS t2
ON t1.anatEntityRelationId = t2.oldRelId;

--
-- Now we can proceed with the update of existing relations
--
UPDATE anatEntityRelation AS t1
INNER JOIN anatEntityRelationTaxonConstraint AS t2
ON t1.anatEntityRelationId = t2.anatEntityRelationId
SET t1.anatEntitySourceId = @newSourceId
WHERE t1.anatEntitySourceId = @oldSourceId and t1.relationStatus != 'reflexive'
AND (t2.speciesId = @speciesId OR t2.speciesId IS NULL);

UPDATE anatEntityRelation AS t1
INNER JOIN anatEntityRelationTaxonConstraint AS t2
ON t1.anatEntityRelationId = t2.anatEntityRelationId
SET t1.anatEntityTargetId = @newTargetId
WHERE t1.anatEntityTargetId = @oldTargetId and t1.relationStatus != 'reflexive'
AND (t2.speciesId = @speciesId OR t2.speciesId IS NULL);