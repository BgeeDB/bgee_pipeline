# really small script that update the values of the displayOrder column in the species table.
# it first order lines by species common name. Then it updates the first lines based on a list 
# of hardcoded species ids. Finally, it update the displayOrder value based on the number of the line
# and write the ordered list of species with displayOrder

cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}
## checking if all necessary arguments were passed....
command_arg <- c("bgeeSpeciesFile", "outputFile")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

species_file <- read.table(bgeeSpeciesFile, sep="\t", h=T, comment.char="", quote="")
#order species_file per common name
species_file_ordered <- species_file[order(species_file$speciesCommonName), ]
# now update top species to keep model organisms on top
manually_ordered_species <- c(9606, 10090, 7955, 7227, 6239)
for (line_to_move in length(manually_ordered_species):1) {
    print(manually_ordered_species[line_to_move])
    row_to_move <- which(as.integer(species_file_ordered$speciesId) == as.integer(manually_ordered_species[line_to_move]))
    moved_row <- species_file_ordered[row_to_move, ]
    species_file_ordered <- species_file_ordered[-row_to_move, ]
    species_file_ordered <- rbind(moved_row, species_file_ordered)
}

#finally update the value of the displayOrder column
species_file_ordered$displayOrder <- 1:nrow(species_file_ordered)

# write the bgeeSpecies file with displayOrder updated
write.table(species_file_ordered, file = outputFile, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
