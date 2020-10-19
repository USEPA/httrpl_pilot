# Functions specific to working with the httr_deg collection schema

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

# getAnlName - Convert analysis options to standardized anl_name for httr_deg docs
getAnlName <- function(mean_cnt=5, plate_effect=F, shrinkage="normal", ...) {
  extra_opts <- list(...)
  anl_name <- sprintf("meanncnt0_%i-plateteffect_%i-shrinkage_%s", mean_cnt, plate_effect, shrinkage)
  for(x in names(extra_opts)) {
    anl_name <- paste(anl_name, paste(x, extra_opts, sep="_"), sep="-")
  }
  return(anl_name)
}

# insertOneDEG - Insert a doc into httr_deg collection, matched by combination of trt_grp_id and anl_name
#
# This also automatically fills in the run_time
# 
# Parameters:
# DB (mongo object) = If specified, use this open connection to httr_deg collection, ignore db_host, db_name, collection
# deg_doc (list) = Single document corresponding to httr_deg schema for insert
# db_host (character) = If DB is NULL, this specifies the host to open a connection to
# db_name (character) = If DB is NULL, this specifies the DB name to connect to
# collection (character) = If DB is NULL, this specifies the collection to connect to, default: httr_deg
# init_note (character) = If this is first insert into DB, set update_notes.note to this value
# update_note (character) = If this is subsequent insert into DB, set last update_notes.note to this value (defaults to: "update on MM/DD/YY")
# rerun (logical) = Whether to re-run the insertion of this data and override existing doc with matching trt_grp_id, anl_name
# debug (logical) = Whether to print debug messages, default: FALSE, overridden by options(debug=...)
#
# Return Value:
# (logical) = TRUE if inserted successfully, FALSE otherwise
insertOneDEG <- function(DB=NULL, deg_doc=NULL, 
                         db_host=NULL, db_name=NULL, collection="httr_deg", 
                         init_note="initial run", update_note=NULL,
                         rerun=FALSE, debug=getOption("debug",default=FALSE)
) {
  # Stop if deg_doc is NULL
  if(is.null(deg_doc)) {
    stop("deg_doc cannot be NULL.\n")
  }
  # Make sure deg_doc has trt_grp_id and anl_name fields
  if(!all(c("trt_grp_id","anl_name") %in% names(deg_doc))) {
    stop("Cannot insert into httr_deg without trt_grp_id and anl_name fields.\n")
  }
  # Open DB connection if an open connection object was not provided
  if(is.null(DB)) {
    stopifnot(!is.null(db_host))
    stopifnot(!is.null(db_name))
    DB = openMongo(host=db_host, db=db_name, collection=collection)
  }
  # Mark timestamp of insert for run_time fields
  ins_time <- Sys.time()
  # Check if this doc is already present
  deg_query <- mongoQuery(deg_doc[c("trt_grp_id","anl_name")])
  if(debug) { cat("Checking for existing DEG entry with query:", deg_query, "in insertOneDEG.\n") }
  match_count <- DB$count(deg_query)
  if(match_count > 1) {
    warning("Multiple documents in httr_deg matched query: ", deg_query, "\n")
  }
  if(match_count > 0) {
    if(rerun) {
      # Get the original update_notes and append to the current update_notes
      prev_update_notes <- DB$iterate(deg_query, fields = mongoQuery('_id'=0, update_notes=1))$one()
      if(is.null(prev_update_notes)) {
        warning("Matching document to be udpated is missing update_notes field - will be added.\n")
        prev_update_notes <- list()
      } else {
        prev_update_notes <- prev_update_notes$update_notes
        # If prev_update_notes is a single-level list (one entry), convert to list of dict type
        if(!is.list(prev_update_notes[[1]])) {
          prev_update_notes <- list(prev_update_notes)
        }
        # Preserve correct datetime formating because mongolite handling of datetimes is a NIGHTMARE!!!
        for(i in 1:length(prev_update_notes)) {
          prev_update_notes[[i]]$run_time <- as.mongo.date(as.POSIXct(prev_update_notes[[i]]$run_time))
        }
      }
      # If update_note is NULL, auto_fill
      if(is.null(update_note)) {
        update_note <- paste("update on", format(ins_time, "%D"))
      }
      prev_update_notes[[length(prev_update_notes)+1]] <- list(
        run_time = as.mongo.date(ins_time),
        note = update_note
      )
      deg_doc[["update_notes"]] <- prev_update_notes
      # Doc exists, and rerun=TRUE, so drop it
      # remove apparently does not return useful info the way insert does
      DB$remove(deg_query)
      if(debug) { cat("Removed", match_count, "documents matching", deg_query, "in httr_deg to be overwritten.\n") }
    } else {
      # Return FALSE
      if(debug) { cat("Document with matching trt_grp_id and anl_name already exists in database, no insert.\n") }
      return(FALSE)
    }
  }
  # Add run_time and update_notes fields if they don't exist already
  if(!("run_time" %in% names(deg_doc))) {
    deg_doc[["run_time"]] <- as.mongo.date(ins_time)
  }
  if(!("update_notes" %in% names(deg_doc))) {
    deg_doc[["update_notes"]] <- list(list(
      run_time = as.mongo.date(ins_time),
      note = init_note
    ))
  }
  # Doc does not exist in DB yet or was removed, go ahead and insert
  ins_count <- DB$insert(deg_doc, auto_unbox=T, na="string")$nInserted
  if(ins_count != 1) {
    warning("nInserted = ", ins_count, " when inserting one document with insertOneDEG; expected 1.\n")
  }
  return(ins_count > 0)
}

# insertManyDEG - Insert multiple docs into httr_deg collection, matched by combination of trt_grp_id and anl_name
# 
# Parameters:
# DB (mongo object) = If specified, use this open connection to httr_deg collection, ignore db_host, db_name, collection
# deg_docs (list of lists) = Each top-level entry should be a document (list) corresponding to httr_deg schema
# db_host (character) = If DB is NULL, this specifies the host to open a connection to
# db_name (character) = If DB is NULL, this specifies the DB name to connect to
# collection (character) = If DB is NULL, this specifies the collection to connect to, default: httr_deg
# rerun (logical) = Whether to re-run the insertion of this data and override existing doc with matching trt_grp_id, anl_name
# debug (logical) = Whether to print debug messages, default: FALSE, overridden by options(debug=...)
# ... = Additional params passed to insertOneDEG
#
# Return Value:
# (integer) = Number of documents successfully inserted
insertManyDEG <- function(DB=NULL, deg_docs=NULL, 
                          db_host=NULL, db_name=NULL, collection="httr_deg", 
                          rerun=FALSE, debug=getOption("debug",default=FALSE), ...
) {
  # Stop if deg_docs is NULL
  if(is.null(deg_docs)) {
    stop("deg_docs cannot be NULL.\n")
  }
  # Make sure deg_docs is a list of httr_deg doc structures, not a single top-level doc
  if(any(c("trt_grp_id","anl_name") %in% names(deg_docs))) {
    stop("Trying to insert a single httr_deg doc with insertManyDEG - use insertOneDEG instead, or pass a list of docs to insertManyDEG.\n")
  }
  # Open DB connection if an open connection object was not provided
  if(is.null(DB)) {
    stopifnot(!is.null(db_host))
    stopifnot(!is.null(db_name))
    DB = openMongo(host=db_host, db=db_name, collection=collection)
  }
  if(debug) { cat("Inserting", length(deg_docs), "docs into httr_deg collection.\n") }
  # Use lapply to insert each doc
  doc_success <- lapply(deg_docs, function(x){
    insertOneDEG(DB=DB, deg_doc=x, rerun=rerun, debug=debug, ...)
  })
  insert_count <- sum(unlist(doc_success))
  if(debug) { cat("Inserted", insert_count, "of", length(deg_docs), "httr_deg docs into database.\n") }
  return(insert_count)
}

# degFrame - Take just the degs member of a httr_deg doc and reformat as a data.frame with appropriate classes for each column
# This a low-level function, is automatically applied by getDEGs when pulling from database.
#
# Parameters:
# degs (list) - The raw degs member of a doc in httr_deg collection, as returned directly from mongo
#
# Return Value;
# (data.frame) - Table with 7 columns matching the underlying deg fields
degFrame <- function(degs) {
  # First make sure to fill in missing pvalue and padj values and convert NAs to numeric type
  degs <- lapply(degs, function(x){
    if(is.null(x$pvalue)) {
      x$pvalue <- as.numeric(NA)
    }
    if(!is.na(x$pvalue) && (x$pvalue == "NA")) {
      x$pvalue <- as.numeric(NA)
    }
    if(is.null(x$padj)) {
      x$padj <- as.numeric(NA)
    }
    if(!is.na(x$padj) && (x$padj == "NA")) {
      x$padj <- as.numeric(NA)
    }
    return(x)
  })
  degs <- do.call(rbind.data.frame, degs)
  degs[,"probe_id"] <- as.character(degs[,"probe_id"])
  for(j in c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")) {
    degs[,j] <- as.numeric(degs[,j])
  }
  return(degs)
}

# flattenDEG - Convert a single httr_deg schema doc to a data.frame
# This is a low-level function that is automatically called by getDEGs whenever flatten=T
# This takes all other fields besides degs and makes them columns in the data.frame
# anl_opt and update_notes are flattened using jsonlite::toJSON
#
# Parameters:
# deg_doc (list) = Single mongo doc conforming to httr_deg schema
# 
# Return Value:
# (data.frame) = Flatted version of deg_doc
flattenDEG <- function(deg_doc) {
  require(jsonlite)
  # If deg_doc$degs is still a list, pass through degFrame first
  if(class(deg_doc$degs)=="list") {
    deg_doc$degs <- degFrame(deg_doc$degs)
  }
  # Append standard fields to front columns
  degTbl <- cbind(
    trt_grp_id=deg_doc$trt_grp_id, 
    anl_name=deg_doc$anl_name,
    anl_opt=as.character(toJSON(deg_doc$anl_opt, auto_unbox = T)),
    run_time=deg_doc$run_time,
    update_notes=as.character(toJSON(deg_doc$update_notes, auto_unbox=T)),
    deg_doc$degs
  )
  # Handle any extra columns that are atomic for now
  extra_cols <- setdiff(names(deg_doc), c("degs", colnames(degTbl)))
  for(j in extra_cols) {
    if(length(deg_doc[[j]])==1) {
      degTbl[,j] <- deg_doc[[j]]
    } else {
      warning("Skipping conversion of httr_deg::", j, " non-atomic field in flattenDEG.\n")
    }
  }
  return(degTbl)
}

# getDEGs - Pull out DESeq2 results from httr_deg table
# 
# Parameters:
# DB (mongo object) = If specified, use this open connection to httr_deg collection, ignore db_host, db_name, collection
# db_host (character) = If DB is NULL, this specifies the host to open a connection to
# db_name (character) = If DB is NULL, this specifies the DB name to connect to
# collection (character) = If DB is NULL, this specifies the collection to connect to, default: httr_deg
# single (logical) = Whether to return a single doc - if TRUE and query matches multiple docs, then only the first match is returned
# flatten (logical) = Whether to flatten all results into a single data frame instead of returning in list structure closer to mongo doc structure
# warn_count (integer) = If query is going to return more than X documents, issue a warning (set to Inf to silence), default: 100
# warn_stype (logical) = When TRUE: If query returns documents with a mix of stype values, issue a warning
# warn_anl (logical) = When TRUE: If query returns documents with a mix of anl_name values, issue a warning
# debug (logical) = Whether to print debug messages, default: FALSE, overridden by options(debug=...)
# ... = All other params are passed to mongoQuery to build the filtering query
#
# Return Value:
# (list or data.frame) = When single & !flatten, this is a list matching httr_deg schema; when !single & !flatten, this is a list of lists with multiple httr_deg docs; when flatten, this is a data.frame
getDEGs <- function(DB=NULL, db_host=NULL, db_name=NULL, collection="httr_deg", 
                    single=F, flatten=F, warn_count=100, warn_stype=T, warn_anl=T,
                    debug=getOption("debug",default=FALSE),
                    ...
) {
  # If flatten & !single, data.table package needed for rbindlist function
  if(flatten && !single) {
    require(data.table)
  }
  
  # Open DB connection if an open connection object was not provided
  if(is.null(DB)) {
    stopifnot(!is.null(db_host))
    stopifnot(!is.null(db_name))
    DB = openMongo(host=db_host, db=db_name, collection=collection)
  }
  # Construct query
  deg_query <- mongoQuery(...)
  if(debug) { cat("Searching for httr_deg docs matching query:", deg_query, "\n") }
  # Count matches to query - warn if multiple matches and single=TRUE
  deg_count <- DB$count(deg_query)
  if(single & (deg_count > 1)) {
    warning("Found ", deg_count, " httr_deg documents matching query in getDEGs(), but single=TRUE, will only return first match.\n")
  } else if(deg_count > warn_count) {
    cat("Pulling", deg_count, "httr_deg documents from DB, this may take some time...\n")
  }
  if(single) {
    deg_result <- DB$iterate(deg_query)$one()
    deg_result$degs <- degFrame(deg_result$degs)
    if(flatten) {
      deg_result <- flattenDEG(deg_result)
    }
    return(deg_result)
  } else {
    # Allow multiple matches, always return a list of docs if !flatten
    # Also time the query process here
    start_time <- Sys.time()
    deg_results <- list()
    deg_iter <- DB$iterate(deg_query)
    while(!is.null(deg_doc <- deg_iter$one())) {
      deg_doc$degs <- degFrame(deg_doc$degs)
      deg_results[[length(deg_results)+1]] <- deg_doc
    }
    end_time <- Sys.time()
    if(debug) { cat("Query returned", length(deg_results), "documents in", as.numeric(difftime(end_time, start_time, units = "secs")), "secs.\n") }
    # Check for multiple stype and anl_name values if those warn options are true
    if(warn_stype) {
      all_stype <- unique(unlist(lapply(deg_results, "[[", "stype")))
      if(length(all_stype) > 1) {
        warning("Query returned ", length(deg_results), " documents with a mix of ", length(all_stype), " different sample types.\n",
                " If this was not expected, try refining your query, e.g. by adding the intended stype.\n",
                " If this was the intended result, you can silence this warning by setting warn_stype = FALSE.\n")
      }
    }
    if(warn_anl) {
      all_anl_name <- unique(unlist(lapply(deg_results, "[[", "anl_name")))
      if(length(all_anl_name) > 1) {
        warning("Query returned ", length(deg_results), " documents with a mix of ", length(all_anl_name), " different analysis configurations.\n",
                " If this was not expected, try adding anl_name or anl_opt fields to your query.\n",
                " If this was the intended result, you can silence this warning by setting warn_anl = FALSE.\n")
      }
    }
    # When flatten, flatten each doc and then rbind together
    if(flatten) {
      # Time this process, report when debug = TRUE
      start_time <- Sys.time()
      n_docs <- length(deg_results)
      deg_results <- lapply(deg_results, flattenDEG)
      # rbindlist is faster than do.call(rbind.data.frame, ...) and the options used here make it tolerant to missing fields in some of the list members
      deg_results <- as.data.frame(rbindlist(deg_results, use.names=T, fill=T))
      end_time <- Sys.time()
      if(debug) { cat("Flattened and merged", n_docs, "httr_deg documents in", as.numeric(difftime(end_time, start_time, units = "secs")), "secs.\n") }
    }
    return(deg_results)
  }
}

# getFCmatrix - Pull out DESeq2 results for multiple contrasts from httr_deg collection and reshape into wide matrix
# 
# Parameters:
# db_host (character) = Specifies the host to open a connection to
# db_name (character) = Specifies the DB name to connect to
# collections (named character vector) = Optional vector to re-map collection names used, where name=default collection; value=new collection name
# anl_name (character) = Which analysis configuration to pull results for, default: based on mean_cnt, plate_effect, and shrinkage params, but setting here takes precedence
# mean_cnt (integer) = passed to getAnlName to set anl_name, default: 5
# plate_effect (logical) = passed to getAnlName to set anl_name, default: T
# shrinkage (character) = passed to getAnlName to set anl_name, default: "normal"
# stype (character) = Filter for sample type, set to NULL to remove this filter, default: "test sample"
# debug (logical) = Whether to print debug messages, default: FALSE, overridden by options(debug=...)
# ... = All other params are passed to mongoQuery to provide additional filtering (if chem_id is not specific enough, this will generate an error)
#
# Return Value:
# (list of data.frame) = First member $fc is matrix of contrast (rows) x probe (col), Second member $cmp is matrix of contrast (rows) x meta-data columns, both in same row order. Missing probes will be NA.
getFCmatrix <- function(db_host, db_name, collections=character(0), 
                        anl_name=getAnlName(mean_cnt=mean_cnt, plate_effect=plate_effect, shrinkage=shrinkage), 
                        mean_cnt=5, plate_effect=T, shrinkage="normal",
                        stype="test sample", debug=getOption("debug",default=FALSE),
                        ...
) {
  # In the future could parallelize the foreach constructs in this function, and add parallelization to getDEGs and other functions above
  require(foreach)
  # Check if any collections are re-mapped, fill in defaults for everything else
  # This could be useful to put in a separate function for all functions that need to access multiple collections
  def_collections <- c("httr_deg", "httr_chem")
  def_collections <- setdiff(def_collections, names(collections))
  names(def_collections) <- def_collections
  collections <- c(collections, def_collections)
  
  # Pull all DEG results matching anl_name, stype, and any additional query terms specified
  deg_dump <- getDEGs(db_host=db_host, db_name=db_name, collection=collections["httr_deg"], 
                      single=F, flatten=F, debug=debug,
                      anl_name=anl_name, stype=stype, ...)
  
  # Result is a list reflecting the individual matching documents in httr_deg
  # First construct the $cmp table with meta-data for each document:
  cmp_cols <- c('trt_grp_id', 'anl_name', 'chem_id', 'ctrl', 'media', 'timeh', 'conc', 'conc_unit', 'dose_level', 'stype', 'pg_id', 'block_id', 'rna_src')
  cmp <- foreach(x=deg_dump, .combine = 'rbind') %do% {
    x_cols <- intersect(names(x), cmp_cols)
    x_row <- as.data.frame(x[x_cols])
    missing_cols <- setdiff(cmp_cols, x_cols)
    for(j in missing_cols) {
      x_row[,j] <- NA
    }
    x_row[,cmp_cols]
  }
  
  # Drop cols that are all NA
  drop_cols <- c()
  for(j in colnames(cmp)) {
    if(all(is.na(cmp[,j]))) {
      cmp[,j] <- NULL
      drop_cols <- c(drop_cols, j)
    }
  }
  if(debug && (length(drop_cols) > 0)) {
    cat(length(drop_cols), "meta-data columns returned no data:", paste(drop_cols, collapse=", "), "\n")
  }
  
  # Check that trt_grp_id are all unique and make these the row.names as well
  if(sum(duplicated(cmp[, "trt_grp_id"])) == 0) {
    row.names(cmp) <- cmp[, "trt_grp_id"]
  } else {
    warning("Database returned multiple deg documents with the same trt_grp_id, may need to refine query parameters to prevent this.\n")
  }
  
  # Get the complete set of probe names (union) across all documents
  all_probes <- foreach(x=deg_dump, .combine='c') %do% {
    x$degs$probe_id
  }
  all_probes <- sort(unique(all_probes))
  if(debug) {
    cat("Collating log2fc matrix for", length(all_probes), "across", length(deg_dump), "contrasts.\n")
  }
  
  # Build a table of contrast (trt_grp_id/doc) x probe of log2fc values - fill NA whenever the probe is missing for a contrast
  fc <- foreach(x=deg_dump, .combine='rbind') %do% {
    # TO DO: Make the error messages more informative here:
    stopifnot(!is.null(x$trt_grp_id))
    stopifnot(class(x$trt_grp_id) == "character")
    stopifnot(!is.null(x$degs))
    stopifnot(class(x$degs) == "data.frame")
    stopifnot(!is.null(x$degs$probe_id))
    stopifnot(class(x$degs$probe_id) == "character")
    stopifnot(sum(duplicated(x$degs$probe_id)) == 0)
    stopifnot(!is.null(x$degs$log2FoldChange))
    stopifnot(class(x$degs$log2FoldChange) == "numeric")
    
    # Extract named vector of log2fc values:
    x_lfc <- setNames(x$degs$log2FoldChange, x$degs$probe_id)
    
    # Resize and re-order to have all-probes (missing probes will get NA values)
    x_lfc <- x_lfc[all_probes]
    names(x_lfc) <- all_probes
    
    # Convert to data.frame, put trt_grp_id as row.name
    x_lfc_row <- as.data.frame(as.list(x_lfc))
    row.names(x_lfc_row) <- x$trt_grp_id
    return(x_lfc_row)
  }
  
  # Debug reporting: How many NA values total?
  if(debug) {
    fc_na <- is.na(fc)
    cat(sum(fc_na), "out of", length(fc_na), "total fold changes in matrix are NA.\n")
  }
  
  # If chem_id column is present in cmp, try querying for additional info in httr_chem - if any is found, map onto cmp
  if("chem_id" %in% colnames(cmp)) {
    cmp_chem_id <- setdiff(unique(cmp$chem_id), "NA")
    if(debug) {
      cat("Attempting to map", length(cmp_chem_id), "chem_id values to additional info in", collections["httr_chem"], "\n")
    }
    httr_chem <- openMongo(host=db_host, db=db_name, collection=collections["httr_chem"])
    chem_info <- httr_chem$find(mongoQuery(chem_id = cmp_chem_id))
    if(nrow(chem_info) == 0) {
      if(debug) {
        cat("No additional info was found for any chem_id values - either", collections["httr_chem"], "does not exist or the chemicals queried are not annotated.\n")
      }
    } else if(any(duplicated(chem_info[, "chem_id"]))) {
      warning(collections["httr_chem"], " returned a table with duplicated chem_id values - mapping of additional chemical info failed as a result.\n")
    } else {
      row.names(chem_info) <- chem_info[,"chem_id"]
      chem_info$chem_id <- NULL
      cat("Found", ncol(chem_info), "annotation columns for", nrow(chem_info), "out of", length(cmp_chem_id), "chemicals in", collections["httr_chem"], "\n")
      # Map onto cmp
      cmp <- cbind(cmp, chem_info[cmp[, "chem_id"], ])
    }
  }
  
  # Need better error reporting here
  stopifnot(all(row.names(fc) == row.names(cmp)))
  
  # Return fc and cmp
  return(list(fc=fc, cmp=cmp))
}
