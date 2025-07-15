## AllianceMine
* For now, the Alliance/AllianceMine only has expression data for wild-type fish under standard/generic control conditions.
* The Alliance/AllianceMine doesn't have Images yet (it does include references/publications and links to ZFIN figure pages).
* Instead of preserving start/end stages as curated by ZFIN, the Alliance opted to expand out the stage information in an expression row in ZFIN from something like this: 1-4 somites to 10-13 somites to this:
    - 1-4 somites
    - 5-9 somites
    - 10-13 somites
* AllianceMine data is loaded from the data that individual model organism databases submit to the Alliance. That data can be found here:
    https://www.alliancegenome.org/downloads



## ZebrafishMine DOES NOT EXIST ANYMORE
* Why do we use a Python script for using intermine on ZFIN?
    => We were unable to make it work with the Perl scripts :/
* Queries with crossReference fields do not work
* "None" is the value for an empty field!
* sex information: "there will be very few direct annotations to ''female organism'' or ''male organism'' for mRNA _in situ_ hybridization experiments."  Leyla Ruzicka <leyla@zfin.org>
* "Wild-type (or standard) zebrafish is used by default. We use the ''WT'' designation if authors say they used a wild-type line but if they don't specify which wild-type line they used. As you know, some of the standard ''wild-type'' lines are not actually wild-type (for example TL https://zfin.org/action/genotype/view/ZDB-GENO-990623-2)."  Leyla Ruzicka <leyla@zfin.org>
