# Functions specific to working with the httr_well collection schema

# NOTE: Code for opening DB collections is in mongo.R - that needs to be loaded as well and the open connection passed here
# Make sure to source this script with chdir=TRUE so it can load other modules from this repository
if(!exists("openMongo")) {
  if(file.exists("mongo.R")) {
    source("mongo.R", chdir=T, local=T)
  } else {
    stop("Could not find Rlib/db/mongo.R - general DB access functions are required here.\n")
  }
}

# Function to pull out counts and treatment info as data.frames for a set of sample IDs
# 
# Parameters:
# DB (mongo object) = If specified, use this open connection to httr_well collection, ignore db_host, db_name, collection
# sample_ids (character vector) = Sample IDs to extract data for, required, cannot be empty
# db_host (character) = If DB is NULL, this specifies the host to open a connection to
# db_name (character) = If DB is NULL, this specifies the DB name to connect to
# collection (character) = If DB is NULL, this specifies the collection to connect to, default: httr_well
# debug (logical) = Whether to print debug messages, default: FALSE, overridden by options(debug=...)
#
# Return Value:
# (list) = Has two members, both data.frames: $treatments has rows = samples, $counts has rows = probes, columns = samples
getWellCounts <- function(DB=NULL, sample_ids=character(0), 
                          db_host=NULL, db_name=NULL, collection="httr_well", 
                          debug=getOption("debug",default=FALSE)
) {
  # Make sure sample_ids defined
  if(length(sample_ids)==0) {
    stop("sample_ids must be specified as a character vector of length >= 1.")
  }
  # Open DB connection if an open connection object was not provided
  if(is.null(DB)) {
    stopifnot(!is.null(db_host))
    stopifnot(!is.null(db_name))
    DB = openMongo(host=db_host, db=db_name, collection=collection)
  }
  well_query <- mongoQuery(sample_id=sample_ids)
  if(debug) {cat("Setting well_query = '", well_query, "' in getWellCounts().\n", sep="")}
  # Run the query with iterate and then process each returned value
  # There might be a way to handle this more elegantly by creating a handler for the find function
  well_iter <- DB$iterate(well_query)
  treatments <- list()
  counts <- list()
  rm_fields <- character(0)
  while(!is.null(well_data <- well_iter$one())) {
    # Put probe_cnts in one list, all other singular fields in another
    well_sid <- well_data$sample_id
    counts[[well_sid]] <- well_data$probe_cnts
    well_data$probe_cnts <- NULL
    # Remove any other non-singular fields
    well_field_sz <- unlist(lapply(well_data, length))
    well_rm_fields <- names(well_field_sz)[well_field_sz > 1]
    for(field in well_rm_fields) {
      well_data[[field]] <- NULL
      rm_fields <- union(rm_fields, field)
    }
    treatments[[well_sid]] <- well_data
  }
  if(length(rm_fields) > 0) {
    warning("Dropped ", length(rm_fields), " fields with nested data: ", paste(rm_fields, collapse=", "), "\n")
  }
  # Get union of field names and types in all treatments tables - infers type from first observed case of each field
  field_types <- character(0)
  for(i in 1:length(treatments)) {
    for(j in names(treatments[[i]])) {
      if(!(j %in% names(field_types))) {
        field_types[j] <- class(treatments[[i]][[j]])
      }
    }
  }
  # Fill in missing values with NA of correct type and re-order to match across all entries in treatments
  treatments <- lapply(treatments, function(trt){
    missing_fields <- setdiff(names(field_types), names(trt))
    for(j in missing_fields) {
      if(field_types[j] == "character") {
        trt[[j]] <- as.character(NA)
      } else if(field_types[j] == "integer") {
        trt[[j]] <- as.integer(NA)
      } else if(field_types[j] == "numeric") {
        trt[[j]] <- as.numeric(NA)
      } else if(field_types[j] == "logical") {
        trt[[j]] <- as.logical(NA)
      } else {
        warning("Don't know how to fill in missing values of type: ", field_types[j], " - will convert to character\n")
        trt[[j]] <- as.character(NA)
      }
    }
    trt[names(field_types)]
  })
  # Collapse the counts and treatments lists into data frames
  treatments <- do.call(rbind.data.frame, treatments)
  # NOTE: This may not properly handle cases where probes are missing in some sample
  counts <- do.call(rbind.data.frame, counts)
  # Transpose counts to probes x samples
  counts <- as.data.frame(t(counts))
  # Make sure counts and treatments both have same sample_ids
  stopifnot(nrow(treatments)==ncol(counts))
  stopifnot(all(row.names(treatments)==colnames(counts)))
  # warn if not all sample_ids found
  if(!all(sample_ids %in% row.names(treatments))) {
    warning("Database only returned data for ", sum(sample_ids %in% row.names(treatments)), " of ", length(sample_ids), " sample_ids.\n")
    # Subset sample_ids to those in both tables
    sample_ids <- intersect(sample_ids, row.names(treatments))
  }
  # Re-order rows and columns to match sample_ids order
  treatments <- treatments[sample_ids,]
  counts <- counts[,sample_ids]
  # Combine into a list and return
  return(list(treatments=treatments, counts=counts))
}

# Function to pull out treatment info ONLY from httr_well as a data.frame for any relevant query
# 
# Parameters:
# DB (mongo object) = If specified, use this open connection to httr_well collection, ignore db_host, db_name, collection
# db_host (character) = If DB is NULL, this specifies the host to open a connection to
# db_name (character) = If DB is NULL, this specifies the DB name to connect to
# collection (character) = If DB is NULL, this specifies the collection to connect to, default: httr_well
# debug (logical) = Whether to print debug messages, default: FALSE, overridden by options(debug=...)
# ... = Any additional parameters are passed to mongoQuery to constrain the query
#
# Return Value:
# (data.frame) = Table of well treatment info from httr_well, currently excludes _id, probe_cnts, raw_id, counts_id, trt_id, sorts by sample_id
getWellInfo <- function(DB=NULL, db_host=NULL, db_name=NULL, collection="httr_well", 
                        debug=getOption("debug",default=FALSE), ...
) {
  # Open DB connection if an open connection object was not provided
  if(is.null(DB)) {
    stopifnot(!is.null(db_host))
    stopifnot(!is.null(db_name))
    DB = openMongo(host=db_host, db=db_name, collection=collection)
  }
  # Pass ... to query
  well_query <- mongoQuery(...)
  if(debug) {cat("Setting well_query = '", well_query, "' in getWellCounts().\n", sep="")}
  # Setup the field filter: exclude _id, probe_cnts, raw_id, counts_id, trt_id
  well_filter <- mongoQuery('_id'=0, probe_cnts=0, raw_id=0, counts_id=0, trt_id=0)
  # Run the query with iterate so we can handle missing data and any other nested fields
  # There might be a way to handle this more elegantly by creating a handler for the find function
  well_iter <- DB$iterate(query = well_query, fields = well_filter)
  treatments <- list()
  rm_fields <- character(0)
  while(!is.null(well_data <- well_iter$one())) {
    well_sid <- well_data$sample_id
    # No probe_cnts because of the filter above
    # Remove any other non-singular fields
    well_field_sz <- unlist(lapply(well_data, length))
    well_rm_fields <- names(well_field_sz)[well_field_sz > 1]
    for(field in well_rm_fields) {
      well_data[[field]] <- NULL
      rm_fields <- union(rm_fields, field)
    }
    treatments[[well_sid]] <- well_data
  }
  if(length(rm_fields) > 0) {
    warning("Dropped ", length(rm_fields), " fields with nested data: ", paste(rm_fields, collapse=", "), "\n")
  }
  # Get union of field names and types in all treatments tables - infers type from first observed case of each field
  field_types <- character(0)
  for(i in 1:length(treatments)) {
    for(j in names(treatments[[i]])) {
      if(!(j %in% names(field_types))) {
        field_types[j] <- class(treatments[[i]][[j]])
      }
    }
  }
  # Fill in missing values with NA of correct type and re-order to match across all entries in treatments
  treatments <- lapply(treatments, function(trt){
    missing_fields <- setdiff(names(field_types), names(trt))
    for(j in missing_fields) {
      if(field_types[j] == "character") {
        trt[[j]] <- as.character(NA)
      } else if(field_types[j] == "integer") {
        trt[[j]] <- as.integer(NA)
      } else if(field_types[j] == "numeric") {
        trt[[j]] <- as.numeric(NA)
      } else if(field_types[j] == "logical") {
        trt[[j]] <- as.logical(NA)
      } else {
        warning("Don't know how to fill in missing values of type: ", field_types[j], " - will convert to character\n")
        trt[[j]] <- as.character(NA)
      }
    }
    trt[names(field_types)]
  })
  # Collapse the treatments list into data.frame
  treatments <- do.call(rbind.data.frame, treatments)
  # Sort by sample_id, use sample_id as row.names, then return
  treatments <- treatments[order(treatments[,"sample_id"], decreasing = F),]
  row.names(treatments) <- treatments[,"sample_id"]
  return(treatments)
}
