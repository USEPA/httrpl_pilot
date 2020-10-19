# Functions specific to working with the httr_probe collection schema

# NOTE: Code for opening DB collections is in mongo.R - that needs to be loaded as well and the open connection passed here
# Make sure to source this script with chdir=TRUE so it can load other modules from this repository
# Shouldn't need to do this if all code is combined into an R package?
if(!exists("openMongo")) {
  if(file.exists("mongo.R")) {
    source("mongo.R", chdir=T, local=T)
  } else {
    stop("Could not find Rlib/db/mongo.R - general DB access functions are required here.\n")
  }
}

# Function to pull out all probe info as a data.frame
# This excludes the transcript field in order to get a flat table, but includes everything else
# The table is also sorted on the index field to ensure the row order of the original manifest
# 
# Parameters:
# DB (mongo object) = If specified, use this open connection to httr_probe collection, ignore db_host, db_name, collection
# db_host (character) = If DB is NULL, this specifies the host to open a connection to
# db_name (character) = If DB is NULL, this specifies the DB name to connect to
# collection (character) = If DB is NULL, this specifies the collection to connect to, default: httr_probe
#
# Return Value:
# (data.frame) = Data frame with complete manifest 
getProbeManifest <- function(DB=NULL, db_host=NULL, db_name=NULL, collection="httr_probe") {
  if(is.null(DB)) {
    stopifnot(!is.null(db_host))
    stopifnot(!is.null(db_name))
    DB = openMongo(host=db_host, db=db_name, collection=collection)
  }
  probeManifest <- DB$find(fields='{"transcripts":0}')
  probeManifest <- probeManifest[order(probeManifest[,"index"], decreasing = F),]
  return(probeManifest)
}
