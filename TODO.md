# TODO for Bgee 15

## Bgee lite

* Add expression rank and expression score info
* Add propagation information:
  * See email 19.11.19 02:04, Tom Conlin
  * We could add the values from `PropagationState` (see `org.bgee.model.expressiondata.Call.ExpressionCall#getDataPropagation()`)
  
## RNA-Seq

* Do not produce absent calls for some gene biotypes, depending on the library type
* Same for the ranks: for now, we consider that all genes that have received
  at least one read in any library are all always accessible to rank computation in all libraries.

## scRNA-Seq

Integrate pipeline code from Sara.
  
