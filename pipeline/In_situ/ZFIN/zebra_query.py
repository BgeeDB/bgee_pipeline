#!/usr/bin/env python3

# This is an automatically generated script to run your query
# to use it you will require the intermine python client.
# To install the client, run the following command from a terminal:
#
#     sudo easy_install intermine
#
# For further documentation you can visit:
#     http://www.intermine.org/wiki/PythonClient

# The following two lines will be needed in every python script:
from intermine.webservice import Service
service = Service("http://zebrafishmine.org/service")

#NOTE "None" is the value for an empty field!
#NOTE sex information: "there will be very few direct annotations to ''female organism'' or ''male organism'' for mRNA in-situ hybridization experiments."  Leyla Ruzicka <leyla@zfin.org>

# Get a new query on the class (table) you will be querying:
query = service.new_query("Gene")

# The view specifies the output columns
query.add_view(
    "primaryIdentifier", "expressions.expressionFound", "symbol",
    "expressions.anatomy.name", "expressions.anatomy.identifier",
    "expressions.startStage.name", "expressions.startStage.identifier",
    "expressions.endStage.name", "expressions.endStage.identifier",
    "expressions.figures.images.primaryIdentifier",
    "expressions.figures.images.label", "organism.taxonId", "expressions.assay",
    "expressions.probe.ThisseCloneRating",
    "expressions.figures.primaryIdentifier",
    "expressions.publication.primaryIdentifier"
)

# Uncomment and edit the line below (the default) to select a custom sort order:
query.add_sort_order("primaryIdentifier", "ASC")

# You can edit the constraint values below
query.add_constraint("expressions.assay", "=", "mRNA in situ hybridization", code = "A")
query.add_constraint("expressions.environment.StandardEnvironment", "=", "true", code = "B")
query.add_constraint("expressions.fish.wildtype", "=", "true", code = "C")                                      #NOTE wild type lines only
#query.add_constraint("expressions.crossReferences.source.name", "CONTAINS", "Ensembl", code = "D")             #FIXME crossReference queries look broken as in FlyMine
#query.add_constraint("expressions.environment.conditions.conditionName", "ONE OF", ["TODO"], code = "E")       # Stable conditions or not
#query.add_constraint("expressions.environment.conditions.primaryIdentifier", "ONE OF", ["TODO"], code = "F")   # Condition identifiers (e.g. ZDB-EXP-041102-1 == (one kind of) normal conditions)
#query.add_constraint("expressions.fish.genotype.primaryIdentifier", "ONE OF", ["TODO"], code = "G")            # Wild type line identifiers

# Uncomment and edit the code below to specify your own custom logic:
# query.set_logic("A and B and C")

# Outer Joins
# (display properties of these relations if they exist,
# but also show objects without these relationships)
#query.outerjoin("expressions.startStage")
#query.outerjoin("expressions.endStage")
query.outerjoin("expressions.figures.images")
query.outerjoin("expressions.probe")
query.outerjoin("expressions.publication")

print("\t".join(map(str, ["#primaryIdentifier", "expressions.expressionFound", "symbol", "expressions.anatomy.name", "expressions.anatomy.identifier", "expressions.startStage.name", "expressions.startStage.identifier",  "expressions.endStage.name", "expressions.endStage.identifier", "expressions.figures.images.primaryIdentifier", "expressions.figures.images.label", "expressions.figures.primaryIdentifier", "expressions.publication.primaryIdentifier", "organism.taxonId", "expressions.assay", "probe.quality"])))
for row in query.rows():
    # Join list of strings and int, stringified, by tab
    print("\t".join(map(str, [row["primaryIdentifier"], row["expressions.expressionFound"], row["symbol"], \
        row["expressions.anatomy.name"], row["expressions.anatomy.identifier"], \
        row["expressions.startStage.name"], row["expressions.startStage.identifier"], \
        row["expressions.endStage.name"], row["expressions.endStage.identifier"], \
        row["expressions.figures.images.primaryIdentifier"], \
        row["expressions.figures.images.label"], row["expressions.figures.primaryIdentifier"], \
        row["expressions.publication.primaryIdentifier"], \
        row["organism.taxonId"], row["expressions.assay"], \
        row["expressions.probe.ThisseCloneRating"]])))

exit()

