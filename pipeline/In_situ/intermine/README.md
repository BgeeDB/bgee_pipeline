Scripts try to be simpler and smarter than templates provided by InterMine:
  * No need of     no warnings ('uninitialized');  # Silence warnings when printing null fields
    because remove warnings for the data AND the code!
    Better to check each data type to replace them by an empty string if uninitialized (// '').
    Should be done automatically in the InterMine API !!!
  * No need of     $, = "\t";  # Set the output field separator as tab
    because it add an extra tab at the end of each output lines.
    Avoid using this global variable!
  * No need to specify a species
    We need taxid for insertion so better to filter by species during insertion
    particularly if several species in queried db.


Could try to use InterMineR package ?!?


**FlyMine**
* Queries with crossReference fields do not work, as well as with mRNAExpressionResults.images.url
  * Do we need mRNAExpressionResults.images.url ???
  * Anyway we don't really need crossReference.source.name CONTAINS 'Ensembl' because
    Ensembl uses FBgn identifiers provided as primaryIdentifier by FlyMine.



**WormMine**
* Queries with crossReference fields do not work
