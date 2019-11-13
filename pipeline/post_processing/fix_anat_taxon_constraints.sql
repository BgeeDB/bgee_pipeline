-- identify problematic anatEntityTaxonConstraint
CREATE TEMPORARY TABLE anatProblem (PRIMARY KEY(anatEntityId, speciesId))
    select DISTINCT t1.speciesId, t1.anatEntityId
    from cond as t1
    left outer join anatEntityTaxonConstraint as t2 on t1.anatEntityId = t2.anatEntityId
        and (t2.speciesId is null or t1.speciesId = t2.speciesId)
    where t2.anatEntityId is null;

-- insert them into anatEntityTaxonConstraint
INSERT INTO anatEntityTaxonConstraint (anatEntityId, speciesId)
    SELECT anatEntityId, speciesId FROM anatProblem;

-- Some anat. entities might now be defined as existing in all species,
-- in which case we use a NULL speciesId in the table anatEntityTaxonConstraint,
-- rather than listing all species. We make the replacement if needed.
-- First, identify anat. entities existing in all species using a list of species
-- rather than a NULL speciesId.
CREATE TEMPORARY TABLE anatInAllSpecies (PRIMARY KEY(anatEntityId))
    SELECT anatEntityId FROM anatEntityTaxonConstraint
    WHERE speciesId IS NOT NULL
    GROUP BY anatEntityId
    HAVING COUNT(DISTINCT speciesId) = (SELECT COUNT(*) FROM species);
-- Insert the NULL speciesId constraints
INSERT INTO anatEntityTaxonConstraint (anatEntityId, speciesId)
    SELECT anatEntityId, NULL FROM anatInAllSpecies;
-- Delete the constraints incorrectly using a list of species
DELETE t1 FROM anatEntityTaxonConstraint AS t1
INNER JOIN anatInAllSpecies AS t2 ON t1.anatEntityId = t2.anatEntityId
WHERE t1.speciesId IS NOT NULL;


-- Now, identify relations related to the problematic anat. entities
-- This table will allow to only insert anatEntityRelation for which "linked" anatEntity is present in the problematic species
CREATE TABLE anatRelProblem (PRIMARY KEY(anatEntityRelationId, speciesId))
    SELECT DISTINCT t1.anatEntityRelationId, t2.speciesId
    FROM anatEntityRelation AS t1
    INNER JOIN anatProblem AS t2
        ON (t1.anatEntitySourceId = t2.anatEntityId OR t1.anatEntityTargetId = t2.anatEntityId)
-- join for filtering purpose only: to retrieve only relations where the "linked" entity exists in the problematic species
    INNER JOIN anatEntityTaxonConstraint AS aetc
        -- part of the ON clause to retrieve taxon constraints of the "linked" entity
        ON (IF(t1.anatEntityTargetId = t2.anatEntityId, t1.anatEntitySourceId, t1.anatEntityTargetId) = aetc.anatEntityId AND
        -- part of the ON clause to check that the "linked" entity does exist in the problematic species
            (t2.speciesId = aetc.speciesId OR aetc.speciesId IS NULL))
-- TODO: modify for not retrieving reflexive is_a part_of
    WHERE t1.anatEntityTargetId != t1.anatEntitySourceId AND
-- Just a safe guard for not retrieving relations already existing in the problematic species
-- (should never happen, but it doesn't hurt)
    NOT EXISTS (
        SELECT 1 FROM anatEntityRelationTaxonConstraint AS t3
        WHERE t3.anatEntityRelationId = t1.anatEntityRelationId
            AND (t3.speciesId IS NULL OR t3.speciesId = t2.speciesId));

-- **WARNING, IMPORTANT STEP:**
-- some relations might be valid only in a given taxon. Inspect the relation taxon constraints
-- with the following query, and remove from anatRelProblem the relation taxon constraints
-- that you do NOT want to be inserted.
SELECT t1.anatEntityRelationId, GROUP_CONCAT(DISTINCT t1.speciesId ORDER BY t1.speciesId) AS speciesIds, t2.speciesId, t3.*
FROM anatEntityRelationTaxonConstraint AS t1
INNER JOIN anatRelProblem AS t2 ON t1.anatEntityRelationId = t2.anatEntityRelationId
INNER JOIN anatEntityRelation AS t3 ON t1.anatEntityRelationId = t3.anatEntityRelationId
GROUP BY t1.anatEntityRelationId;

-- **REMOVAL STEPS. CHECK THEM BEFORE DELETING**
-- Queries used to delete instances from the anatRelProblem table for bgee14
-- Check the queries to be sure that they are still correct
-- delete from anatRelProblem all anatEntityRelationId for 7955 (zebrafish) when there is 2 or less
-- already existing taxonConstraints with one being 9606 or 10090
DELETE FROM anatRelProblem WHERE anatEntityRelationId IN
    (SELECT DISTINCT t1.anatEntityRelationId FROM anatEntityRelationTaxonConstraint AS t1
    INNER JOIN (select * from anatRelProblem) AS t2 ON t1.anatEntityRelationId = t2.anatEntityRelationId
    INNER JOIN anatEntityRelation AS t3 ON t1.anatEntityRelationId = t3.anatEntityRelationId
    WHERE t2.speciesId = 7955 GROUP BY t1.anatEntityRelationId
    HAVING COUNT(DISTINCT t1.speciesId) <= 2
    	AND ((GROUP_CONCAT(DISTINCT t1.speciesId) LIKE "%9606%") OR (GROUP_CONCAT(DISTINCT t1.speciesId) LIKE "%10090%")));
-- delete from anatRelProblem all anatEntityRelationId for 13616 (Monodelphis domestica) when there is only 1
-- already existing taxonConstraint for 9606
DELETE FROM anatRelProblem WHERE anatEntityRelationId IN
    (SELECT DISTINCT t1.anatEntityRelationId FROM anatEntityRelationTaxonConstraint AS t1
    INNER JOIN (select * from anatRelProblem) AS t2 ON t1.anatEntityRelationId = t2.anatEntityRelationId
    INNER JOIN anatEntityRelation AS t3 ON t1.anatEntityRelationId = t3.anatEntityRelationId
    WHERE t2.speciesId = 13616 GROUP BY t1.anatEntityRelationId
    HAVING COUNT(DISTINCT t1.speciesId) = 1
    	AND (GROUP_CONCAT(DISTINCT t1.speciesId) LIKE "%9606%"));
-- delete from anatRelProblem all anatEntityRelationId for 10090 (Mus musculus) when the only already
-- already existing taxonConstraint is for 9606
DELETE FROM anatRelProblem WHERE anatEntityRelationId IN
    (SELECT DISTINCT t1.anatEntityRelationId FROM anatEntityRelationTaxonConstraint AS t1
    INNER JOIN (select * from anatRelProblem) AS t2 ON t1.anatEntityRelationId = t2.anatEntityRelationId
    INNER JOIN anatEntityRelation AS t3 ON t1.anatEntityRelationId = t3.anatEntityRelationId
    WHERE t2.speciesId = 10090 GROUP BY t1.anatEntityRelationId
    HAVING COUNT(DISTINCT t1.speciesId) = 1
    	AND (GROUP_CONCAT(DISTINCT t1.speciesId) LIKE "%9606%"));
-- delete from anatRelProblem all anatEntityRelationId for 10116 (Rattus) when the only already
-- existing taxonConstraint is for 9606
DELETE FROM anatRelProblem WHERE anatEntityRelationId IN
    (SELECT DISTINCT t1.anatEntityRelationId FROM anatEntityRelationTaxonConstraint AS t1
    INNER JOIN (select * from anatRelProblem) AS t2 ON t1.anatEntityRelationId = t2.anatEntityRelationId
    INNER JOIN anatEntityRelation AS t3 ON t1.anatEntityRelationId = t3.anatEntityRelationId
    WHERE t2.speciesId = 10116 GROUP BY t1.anatEntityRelationId
    HAVING COUNT(DISTINCT t1.speciesId) = 1
    	AND (GROUP_CONCAT(DISTINCT t1.speciesId) LIKE "%9606%"));
-- delete from anatRelProblem all anatEntityRelationId for 9031 (Gallus gallus) when the only already
-- existing taxonConstraint is for 9606
DELETE FROM anatRelProblem WHERE anatEntityRelationId IN
    (SELECT DISTINCT t1.anatEntityRelationId FROM anatEntityRelationTaxonConstraint AS t1
    INNER JOIN (select * from anatRelProblem) AS t2 ON t1.anatEntityRelationId = t2.anatEntityRelationId
    INNER JOIN anatEntityRelation AS t3 ON t1.anatEntityRelationId = t3.anatEntityRelationId
    WHERE t2.speciesId = 9031 GROUP BY t1.anatEntityRelationId
    HAVING COUNT(DISTINCT t1.speciesId) = 1
    	AND (GROUP_CONCAT(DISTINCT t1.speciesId) LIKE "%7955%"));


-- Insert missing relation taxon constraints
INSERT INTO anatEntityRelationTaxonConstraint (anatEntityRelationId, speciesId)
    SELECT anatEntityRelationId, speciesId FROM anatRelProblem;

-- Insert reflexive "is_a part_of" anatEntityRelationTaxonConstraints for problematic anatEntityTaxonConstraint
INSERT INTO anatEntityRelationTaxonConstraint (anatEntityRelationId, speciesId)
    SELECT distinct t2.anatEntityRelationId, t1.speciesId
    FROM anatProblem AS t1
    INNER JOIN anatEntityRelation AS t2
        ON (t1.anatEntityId = t2.anatEntitySourceId AND t1.anatEntityId = t2.anatEntityTargetId)
    LEFT OUTER JOIN anatEntityRelationTaxonConstraint AS t3 ON (t2.anatEntityRelationId = t3.anatEntityRelationId AND (t1.speciesId = t3.speciesId OR t3.speciesId IS NULL))
    WHERE t3.anatEntityRelationId IS NULL;

-- Some relations might now be defined as existing in all species,
-- in which case we use a NULL speciesId in the table anatEntityRelationTaxonConstraint,
-- rather than listing all species. We make the replacement if needed.
-- First, identify relations existing in all species using a list of species
-- rather than a NULL speciesId.
CREATE TEMPORARY TABLE anatRelInAllSpecies (PRIMARY KEY(anatEntityRelationId))
    SELECT anatEntityRelationId FROM anatEntityRelationTaxonConstraint
    WHERE speciesId IS NOT NULL
    GROUP BY anatEntityRelationId
    HAVING COUNT(DISTINCT speciesId) = (SELECT COUNT(*) FROM species);
-- Insert the NULL speciesId constraints
INSERT INTO anatEntityRelationTaxonConstraint (anatEntityRelationId, speciesId)
    SELECT anatEntityRelationId, NULL FROM anatRelInAllSpecies;
-- Delete the constraints incorrectly using a list of species
DELETE t1 FROM anatEntityRelationTaxonConstraint AS t1
INNER JOIN anatRelInAllSpecies AS t2 ON t1.anatEntityRelationId = t2.anatEntityRelationId
WHERE t1.speciesId IS NOT NULL;
