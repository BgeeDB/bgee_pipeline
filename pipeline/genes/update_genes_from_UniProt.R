library(UniProt.ws)
library(RMySQL)
library(foreach)
library(doParallel)

## Session info
print(sessionInfo())

## Reading in arguments provided in command line
cmd_args = commandArgs(TRUE)
print (cmd_args)
if (length(cmd_args)==0){
  stop("no arguments provided\n")
} else {
  for (i in 1:length(cmd_args)){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed in command line
command_arg <- c("species_file", "bgee", "number_cores", "output_file")

for (c_arg in command_arg){
  if (!exists(c_arg)){
    stop(paste(c_arg, "  command line argument not provided\n"))
  } else {
    cat(paste(c_arg,":\t", eval(c_arg),sep=""),"\n")
  }
}

################################# FUNCTIONS ######################################

hack_taxon_id <- function(bgeeSpecies) {
  bgeeSpecies$uniprotTaxId <- bgeeSpecies$speciesId
  for (i in seq_len(nrow(bgeeSpecies))) {
    if(as.numeric(bgeeSpecies$speciesId[i]) == 9593){
      bgeeSpecies$uniprotTaxId <- 9595
    } else if (as.numeric(bgeeSpecies$speciesId[i]) == 7237){
      bgeeSpecies$uniprotTaxId <- 46245
    } 
  }
  return(bgeeSpecies)
}
connect_to_db <- function(bgee) {
  splited_bgee <- as.data.frame(strsplit(x = bgee, "__"))
  parameters <- NULL
  for (i in seq_len(nrow(splited_bgee))) {
    parameters <- rbind(parameters,as.data.frame(t(as.data.frame(strsplit(splited_bgee[i,],"=")))))
  }
  user <- toString(parameters[parameters[,1] == "user",][2])
  pwd <- toString(parameters[parameters[,1] == "pass",][2])
  host <- toString(parameters[parameters[,1] == "host",][2])
  name <- toString(parameters[parameters[,1] == "name",][2])
  return(dbConnect(MySQL(), user=user, password=pwd, dbname=name, host=host))
}

select_bgee_genes <- function(speciesId, mydb) {
  query <- paste0("SELECT DISTINCT t1.geneId FROM gene AS t1 WHERE NOT EXISTS (SELECT 1 FROM 
  geneXRef AS t2 WHERE t1.bgeegeneId = t2.bgeegeneId AND t2.dataSourceId IN (4,5)) 
  AND t1.speciesId = ",speciesId)
  return(dbGetQuery(mydb, query))
}

retrieve_uniprot_mapping <- function(speciesId, dataSourceId, geneIDs) {
  #init uniprot object specific to the species
  uniprot_object <- NULL
  # Sometimes there is an error when trying to create the uniprot object
  # because of a timeout at uniprot. The tryCatch allows to try 10 times before crashing
  already_tried <- 0
  while(is.null(uniprot_object) & already_tried < 10) {
    tryCatch( {
      uniprot_object <- UniProt.ws(taxId = speciesId)
    },
    error=function(cond) {
      Sys.sleep(5)
      already_tried <- already_tried + 1
    })
  }
  
  # init variable of the uniprot query taking into consideration the 
  # dataSource of the species
  query_columns <- c("UNIPROTKB", "UNIPROTKB_ID", "REVIEWED","GENENAME")
  if (dataSourceId == ENSEMBL_DS_ID) {
    query_keytype <- "ENSEMBL"
  } else if (dataSourceId == ENSEMBLMETAZOA_DS_ID) {
    query_keytype <- "ENSEMBL_GENOMES"
  } else if (dataSourceId == NCBI_DS_ID) {
    query_keytype <- "GENEID"
  } else {
    stop("dataSrouce not recognized")
  }
  # hardcoded for now but could be an attribute of the script
  chunk <- 99
  
  current <- 1
  species_mapping <- NULL
  while (current < nrow(geneIDs)) {
    tryCatch(
      {
        current_values <- select(x = uniprot_object, 
          keys = geneIDs$geneId[current:(current+chunk)], 
          columns = query_columns, keytype = query_keytype)
        species_mapping <-rbind(species_mapping, current_values)
      },
      error=function(cond) {
        warning("no information for the chunk")
      },
      finally={
        if(current + chunk > nrow(geneIDs)) {
          current <- nrow(geneIDs)
        }else {
          current <- current + chunk
        }
      }
    )
  }
  if(!is.null(species_mapping)) {
    colnames(species_mapping)[1] <- "GeneID"
    return(na.omit(species_mapping))
  } else {
    species_mapping <- data.frame(matrix(vector(), 0, 5,
      dimnames=list(c(), c("GeneID", query_columns))),
      stringsAsFactors=F)
    return(species_mapping)
  }
}

############################## MAIN CODE ##################################

# create constants that could have to be updated if dataSourceId change 
# in the bgeeSpecies.tsv file. Not exported as attribute of the script
# because they are hardcoded in the file describing bgee species so they 
# should be stable
NCBI_DS_ID = 37
ENSEMBL_DS_ID = 2
ENSEMBLMETAZOA_DS_ID = 24

#read species file
bgee_species <- read.table(file = species_file, sep = "\t", header = TRUE)

# Some NCBI taxonId are not the same in Uniprot than in Bgee. We created a function
# to hack these mismatches and allows retrieval of Uniprot XRefs 
bgee_species <- hack_taxon_id(bgee_species)

start_time <- Sys.time()

#initialize the number of cores used to download data from uniprot
cluster_object <- parallel::makeCluster(number_cores)
# register the cluster for using foreach
doParallel::registerDoParallel(cluster_object)

# retrieve data from UniProt for each species
all_species <- foreach(x = iter(bgee_species, by="row"), .combine = "rbind") %do% {
  library(RMySQL)
  library(UniProt.ws)
  mydb <- connect_to_db(bgee)
  gene_ids <- select_bgee_genes(x["speciesId"], mydb)
  dbDisconnect(mydb)
  species_mapping <- suppressWarnings(retrieve_uniprot_mapping(speciesId = as.numeric(x["uniprotTaxId"]), 
    dataSourceId = as.numeric(x["dataSourceId"]), geneIDs = gene_ids))
  if(nrow(species_mapping) > 0) {
    species_mapping$SpeciesId <- as.numeric(x["speciesId"])
    species_mapping$DataSourceId <- as.numeric(x["dataSourceId"])
  } else {
    species_mapping <- data.frame(matrix(vector(), 0, 7,
        dimnames=list(c(), c(colnames(species_mapping), "SpeciesId", "DataSourceId"))),
        stringsAsFactors=F)
  }
  return(species_mapping)
}
parallel::stopCluster(cluster_object)

end_time <- Sys.time()
cat(end_time - start_time, " seconds")

#write output file in order to insert data in the database using a perl script
write.table(x=all_species, file = output_file, 
  sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)


