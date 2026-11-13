-- TODO enable multispecies groups when necessary
-- TODO fix insertion of download files

-- query to get the single species data groups for easier generation of this file:
SELECT speciesDisplayOrder,
        IF(speciesCommonName IS NOT NULL AND speciesCommonName != '',
           speciesCommonName, CONCAT(genus, ' ', species)),
        speciesDisplayOrder
FROM species ORDER BY speciesDisplayOrder;

-- query to get member of a given taxon, to help build the multi-species groups
SELECT t3.* from taxon as t1 inner join taxon as t2 
    ON t2.taxonLeftBound >= t1.taxonLeftBound
    AND t2.taxonRightBound <= t1.taxonRightBound
INNER JOIN species AS t3 on t3.taxonId = t2.taxonId
WHERE t1.taxonId = 1437010;

DELETE FROM downloadFile;
DELETE from speciesDataGroup;


INSERT INTO speciesDataGroup (speciesDataGroupId, speciesDataGroupName, speciesDataGroupDescription,
speciesDataGroupOrder) VALUES
-- single-species data groups
(1, 'human', NULL, 1),
(2, 'mouse', NULL, 2),
(3, 'zebrafish', NULL, 3),
(4, 'fruit fly', NULL, 4),
(5, 'nematode', NULL, 5),
(6, 'dog', NULL, 6),
(7, 'cat', NULL, 7),
(8, 'horse', NULL, 8),
(9, 'pig', NULL, 9),
(10, 'cattle', NULL, 10),
(11, 'Goat', NULL, 11),
(12, 'Sheep', NULL, 12),
(13, 'rabbit', NULL, 13),
(14, 'guinea pig', NULL, 14),
(15, 'chicken', NULL, 15),
(16, 'Wild turkey', NULL, 16),
(17, 'platypus', NULL, 17),
(18, 'White-tufted-ear marmoset', NULL, 18),
(19, 'Sooty mangabey', NULL, 19),
(20, 'Crab-eating macaque', NULL, 20),
(21, 'macaque', NULL, 21),
(22, 'Pig-tailed macaque', NULL, 22),
(23, 'Olive baboon', NULL, 23),
(24, 'gorilla', NULL, 24),
(25, 'bonobo', NULL, 25),
(26, 'chimpanzee', NULL, 26),
(27, 'Gray mouse lemur', NULL, 27),
(28, 'Green monkey', NULL, 28),
(29, 'Malayan pangolin', NULL, 29),
(30, 'rat', NULL, 30),
(31, 'Naked mole rat', NULL, 31),
(32, 'opossum', NULL, 32),
(33, 'African clawed frog', NULL, 33),
(34, 'western clawed frog', NULL, 34),
(35, 'green anole', NULL, 35),
(36, 'Spotted gar', NULL, 36),
(37, 'European freshwater eel', NULL, 37),
(38, 'Blind cave fish', NULL, 38),
(39, 'Northern pike', NULL, 39),
(40, 'Atlantic salmon', NULL, 40),
(41, 'Atlantic cod', NULL, 41),
(42, 'Guppy', NULL, 42),
(43, 'Japanese rice fish', NULL, 43),
(44, 'Eastern happy', NULL, 44),
(45, 'Fairy cichlid', NULL, 45),
(46, 'Turbot', NULL, 46),
(47, 'Three-spined stickleback', NULL, 47),
(48, 'Turquoise killifish', NULL, 48),
(49, 'Common lancelet', NULL, 49),
(50, 'Coelacanth', NULL, 50),
(51, 'Drosophila pseudoobscura', NULL, 51),
(52, 'Drosophila simulans', NULL, 52);

INSERT INTO speciesToDataGroup (speciesDataGroupId, speciesId) VALUES 
-- single-species data groups
(1, 9606),
(2, 10090),
(3, 7955),
(4, 7227),
(5, 6239),
(6, 9615),
(7, 9685),
(8, 9796),
(9, 9823),
(10, 9913),
(11, 9925),
(12, 9940),
(13, 9986),
(14, 10141),
(15, 9031),
(16, 9103),
(17, 9258),
(18, 9483),
(19, 9531),
(20, 9541),
(21, 9544),
(22, 9545),
(23, 9555),
(24, 9593),
(25, 9597),
(26, 9598),
(27, 30608),
(28, 60711),
(29, 9974),
(30, 10116),
(31, 10181),
(32, 13616),
(33, 8355),
(34, 8364),
(35, 28377),
(36, 7918),
(37, 7936),
(38, 7994),
(39, 8010),
(40, 8030),
(41, 8049),
(42, 8081),
(43, 8090),
(44, 8154),
(45, 32507),
(46, 52904),
(47, 69293),
(48, 105023),
(49, 7740),
(50, 7897),
(51, 7237),
(52, 7240);

