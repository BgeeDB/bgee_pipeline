Information for work related to GenoFish project with Tina Begum

## Retrieval of homologous anatomical entities

Tina provided a list of tissues with data in fishes, need to identify homologous tissue in Euteleostomi.

### Original request

Tissue data in zebrafish, cavefish, medaka, spotted gar, and northern pike: brain, embryo, gills,
heart, bones, intestine, kidney, liver, muscle, ovary and testis.
Need to find homologous organs in human, mouse, cow, rat, macaque, gorilla and chimpanzee.

### Identifiers retrieved for the provided list of tissues:

* brain: UBERON:0000955
* embryo: UBERON:0000922
* pharyngeal gill: UBERON:0000206
* heart: UBERON:0000948
* bone element: UBERON:0001474
* intestine: UBERON:0000160
* kidney: UBERON:0002113
* liver: UBERON:0002107
* muscle organ: UBERON:0001630
* ovary: UBERON:0000992
* testis: UBERON:0000473

### Output

Two files are generated: one containing the homologous entities existing in at least one of the provided species;
one containing the list of anatomical entities for which no homology relation existed, or that did not exist
in any of the requested species.

See [output directory](../../../generated_files/collaboration/genofish/).

### Pipeline details

Retrieving homologous tissues: using Java code `org.bgee.pipeline.expression.GenoFishProject`
in `bgee-pipeline` module. See `homology` step in Makefile. List of species used:
`7955,9606,10090,9913,10116,9544,9593,9598` (currently we don't have in Bgee the species cavefish,
medaka, spotted gar, and northern pike).