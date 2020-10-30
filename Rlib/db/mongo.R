# R versions of main mongo access functions in lib/db/mongo.py
# Mongo API here is provided through the mongolite R package, does not have the full functionality of pymongo package
# The main difference is that mongolite connections are opened separately for each collection, whereas pymongo connections are made once per DB
# There are also some issues with type differences in R

# This is needed for the URLencode function, which protects passwords in mongoURL
# This package should be part of the core R installation
library(utils)

# Load the passwords.R module for credentials management:
if(!exists("getCredentials")) {
  if(file.exists("passwords.R")) {
    source("passwords.R", chdir=T, local=T)
  } else {
    stop("Could not find Rlib/db/passwords.R - password management functions are required here.\n")
  }
}

# Just format a mongo URL based on mongoServer, username, password, and database
# This currently does no validation to check for problems
# It also does not currently support anonymous access (no username/password)
#        A generally better approach might be to just have a ... catch-all param for all extended variables on the mongoURL, convert that to param=var&... string
#        And as long as that string is not empty, append it with ? separator
mongoURL <- function(host, user, passwd, db, authSource=NULL, authMechanism=NULL) {
  passwd <- URLencode(passwd, reserved = TRUE)  # Proect any special characters in passwd
  if(is.null(authSource) & is.null(authMechanism)){
    paste0("mongodb://", user, ":", passwd, "@", host, "/", db)
  }else {
    paste0("mongodb://", user, ":", passwd, "@", host, "/", db, "?", "authSource=", authSource, "&", "authMechanism=", authMechanism)
  }
}

# openMongo - open a connection to a specific collection in a mongo db
# NOTE: A major difference between pymongo and mongolite interfaces is that mongolite connects to a specific collection, not the entire DB
openMongo <- function(host=getOption("httrDefaultHost"), user=NULL, passwd=NULL,
                      db=NULL, collection=NULL, authSource = NULL, authMechanism = NULL) {
  require(mongolite)
  if(is.null(user) | is.null(passwd)) {
    myCred <- getCredentials(host=host, db=db)
    if(is.null(user)) {
      user <- myCred$user
    }
    if(is.null(passwd)) {
      passwd <- myCred$passwd
    }
    if(is.null(authSource)) {
      authSource <- myCred$authSource
    }
    if(is.null(authMechanism)) {
      authMechanism <- myCred$authMechanism
    }
  }
  # Open the mongo connection
  mongo(collection=collection, url=mongoURL(host=host, user=user, passwd=passwd, db=db, authSource = authSource, authMechanism = authMechanism), verbose=getOption("verbose"))
}

# getDocIDs - extract a certain ID field (e.g. sample_id) out of a list of documents
# This is designed to mimic the function of the same name from lib/db/mongo.py
# Parameters:
#  docs (list) = List of documents dumped out from MongoDB as a list type
#  id_field (character) = Name of field to pull out of every member of docs
# Return Value:
#  (vector) = Vector of ID values from docs, atomic type depends on the type of ID (usually character)
getDocIDs <- function(docs, id_field="sample_id") {
  return(unlist(lapply(docs, "[[", id_field)))
}

# findByID - Find all documents in a collection with optional ID filtering.
# Given a vector of IDs (typically either _id or sample_id), return an iterable object (mongo iterator or list) corresponding to all 
# matching documents in a collection. If ids is an empty list, return all documents in the collection. Also has options to limit the result 
# to just the id_field or any other subset of fields, and to dump out the data as a list (stored in local memory) instead of returning the 
# pymongo Cursor object.
#
# Parameters:
#  DB (mongo object) = Open connection to specific Mongo Database and Collection
#  ids (vector) = Vector of IDs (usually character) to query on, gets all values if empty, default: empty vector.
#  id_field (character) = Which ID field to match against ids, default: sample_id.
#  id_only (logical) = If True, return only the values in the ID field being matched against, otherwise return the whole document
#  fields (Mongo query, as character) = If id_only = False, this will be used to filter the return fields instead.
#  dump (logical) = If True, loop over the initial return object (pymongo.cursor.Cursor) and extract all results to a new list object.
#  debug (logical) = If True, dump out some extra debugging messages, can also set using options(debug=True)
#
# Return Value:
#  (mongo iterator, list, or vector) = Either:
#   an iterator object to return results from Mongo (dump=False), 
#   a list of all results (dump=True, id_only=False)
#   a vector of IDs (dump=True, id_only=True)
#
findByID <- function(DB, ids=c(), id_field="sample_id", id_only=FALSE, fields='{}', dump=FALSE, debug=getOption("debug",default=FALSE)) {
  require(mongolite)
  require(jsonlite)
  # Set query based on ids, id_field
  if (length(ids) > 0) {
    my_query = list(list('$in'=ids))
    names(my_query)[1] <- id_field[1]
    my_query = toJSON(my_query)
  } else {
    my_query = '{}'
  }
  if(debug) {cat("Setting my_query = '", my_query, "'\n", sep="")}
  # Set filter based on id_only or fields
  if(id_only) {
    if(id_field == "_id") {
      my_fields = paste0('{ "', id_field, '": true }')
    } else {
      my_fields = paste0('{ "_id": false, "', id_field, '": true }')
    }
  } else {
    # TO DO: Might need to set the default for this to '{}'?
    my_fields = fields
  }
  if(debug) {cat("Setting my_fields = '", my_fields, "'\n", sep="")}
  # Run the query with appropriate query and filter, return iterator object
  query_res = DB$iterate(query=my_query, fields=my_fields)
  # Dump to list?
  if(dump) {
    query_dump <- list()
    while(!is.null(x <- query_res$one())) {
      if(id_only) {
        query_dump[[length(query_dump)+1]] <- x[[id_field]]
      } else {
        query_dump[[length(query_dump)+1]] <- x
      }
    }
    if(id_only) {
      query_dump <- unlist(query_dump)
    }
    return(query_dump)
  } else {
    return(query_res)
  }
}

# splitIDs - Split a list of IDs based on which are present/absent in a collection.
# Given a list of IDs, search against a particular collection, and return a dict with two lists keyed by "present" and "absent".
#
# Parameters:
#  DB (mongo object) = Open connection to a Mongo DB and Collection
#  ids (vector) = Vector of IDs (typically character)
#  id_field (character) = Name of ID field, default: sample_id
#
# Return Value:
#  (list) with two members:
#   $present = vector of IDs that were found in the collection
#   $absent  = vector of IDs that were NOT found in collection
#
splitIDs <- function(DB, ids, id_field="sample_id") {
  # First get the set of matching IDs that are in the DB
  found_id <- findByID(DB, ids=ids, id_field=id_field, id_only=TRUE, dump=TRUE)
  # Split by which ones overlap, return as two member list
  my_res <- list(
    present = ids[ids %in% found_id],
    absent = ids[!(ids %in% found_id)]
  )
  return(my_res)
}

# deleteByID - Delete documents from a collection.
# Given a list of IDs, delete all matching documents from a collection, or delete ALL documents. The primary purpose is to remove documents 
# from a collection before replacing them when rerun=True in a particular pipeline step.
#
# Parameters:
#  DB (mongo object) = Open connection to database and collection
#  ids (vector) = Vector (usually character) of IDs for documents to delete, default=c(), which defers to delete_all
#  id_field (character) = Alternate ID field to use, default: "sample_id"
#  delete_all (logical) = Setting to True will delete the whole collection, default: False
#  debug (logical) = Whether to post debug messages
# 
# Return Value:
#  (int) = Number of documents that were deleted, will also warn if != length(ids) when delete_all=False
deleteByID <- function(DB, ids=c(), id_field="sample_id", delete_all=FALSE, debug=getOptions("debug",default=FALSE)) {
  require(mongolite)
  require(jsonlite)
  # If delete_all, query is {}, otherwise construct from ids, id_field
  if(delete_all) {
    my_query = '{}'
  } else {
    # Otherwise, build query for user-specified id list
    my_query = list(list('$in'=ids))
    names(my_query)[1] <- id_field[1]
    my_query = toJSON(my_query)
    # Warn if ids=[]
    if(length(ids) == 0) {
      warning("deleteByID called with ids=c() and delete_all=False, no documents will be deleted.\n")
    }
  }
  # Debug message:
  if(debug) {
    cat("Deleting documents with",id_field,"matching",length(ids),"IDs.\n")
  }
  # Call $remove to do the actual deletion - use $count to make sure the right number of documents were removed
  # TO DO: Does $remove return a list of counts similar to what $insert returns? If so, should use that to get del_cnt and make sure other numbers are 0?
  prev_cnt = DB$count()
  delete_res = DB$remove(query=my_query)
  new_cnt = DB$count()
  del_cnt = prev_cnt - new_cnt
  if((del_cnt != length(ids)) && !delete_all) {
    warning("Called deleteByID with ", length(ids), " IDs but ", del_cnt, " documents were deleted (id_field='",id_field,"').\n")
  }
  return(del_cnt)
}

# insertByID - Insert a list of documents into a MongoDB collection without creating redundant entries.
# Takes a list of dicts and first checks whether any documents with matching IDs are already present. If so, these documents are either 
# skipped for insert (rerun=False) or replaced (rerun=True)
#
# Parameters:
#  DB (mongo object) = Open database/collection connection for inserting
#  docs (list) = List of documents to insert into collection
#  id_field (character) = The field to use for matching up documents, default: sample_id
#  rerun (logical) = Whether to overwrite existing data with same sample_id
#  debug (logical) = Whether to report debug messages
# 
# Return Value:
#  (integer) = The number of documents successfully inserted, 
#  Note: Unlike the pymongo API, mongolite does not return the IDs after a successful insert
#
# NOTE: For this function in particular, it makes sense to have a custom wrapper for each collection, e.g. to make sure there are matching IDs in other collections
insertByID <- function(DB, docs, id_field="sample_id", rerun=FALSE, debug=getOption("debug",default=FALSE)) {
  require(mongolite)
  # TO DO: Convert this python code to R
  # Check if any matching sample_id in the collection already
  id_status = splitIDs(DB, ids=getDocIDs(docs, id_field=id_field), id_field=id_field)
  # Check if any of these sample_ids are already in the DB:
  if(length(id_status$present) > 0) {
    if(rerun) {
      # If rerun = True, drop these IDs from the DB before final insert below
      deleteByID(DB, ids=id_status['present'], id_field=id_field, debug=debug)
    } else {
      # rerun = False, so only insert the 'absent' values
      if(length(id_status$absent) == 0) {
        warning("Collection already contains documents with ", id_field, " matching all ", length(docs), " entries in docs, no DB insert to be done (rerun=False).")
        return(0)
      }
      # Otherwise, filter docs to just those with sample_id matching absent
      subset_docs = sapply(docs, function(x){x[[id_field]] %in% id_status$absent})
      docs = docs[subset_docs]
      warning("Collection already contains ", length(id_status$present), " matching ", id_field, ", only ", length(docs), " new documents will be inserted (rerun=False)")
    }
  }
  if(debug) {
    cat("Inserting",length(docs),"new documents into collection, rerun =", rerun, "\n")
  }
  
  # Can use the list object returned by each insert statement to make sure the right number of docs went in
  # prev_count = DB$count()
  
  # Insert into the database - mongo$insert does NOT handle lists of documents so need to use lapply to insert each document
  # Alternatively, might be able to convert each document to a JSON string first and pass a vector of those strings to insert all at once?
  insert_return <- lapply(docs, function(x){DB$insert(x, auto_unbox=T)})
  # Summarize the insert results
  insert_summary <- list(
    nInserted = sum(unlist(lapply(insert_return, function(x){x$nInserted}))),
    nMatched = sum(unlist(lapply(insert_return, function(x){x$nMatched}))),
    nRemoved = sum(unlist(lapply(insert_return, function(x){x$nRemoved}))),
    nUpserted = sum(unlist(lapply(insert_return, function(x){x$nUpserted}))),
    nWriteErrors = sum(unlist(lapply(insert_return, function(x){length(x$writeErrors)})))
  )
  
  ins_cnt = insert_summary$nInserted
  # TO DO: Should also check that other count fields in insert_summary are all 0 and issue warnings if not?
  
  # One final check to make sure the correct number of docs were inserted
  if(ins_cnt != length(docs)) {
    warning("Attempted to insert ", length(docs), " documents, but collection only increased by ", ins_cnt, " documents.\n")
  }
  return(ins_cnt)
}

# mongoQuery - Build a mongo query from an arbitrary set of params or a list
# Takes a list of named values/vectors and converts to a JSON-format mongo query. Each entry with length > 1 is interpreted
# as a set of potential matching values and converted to the $in keyword for mongo
#
# Parameters:
#  ... (any) = Any set of named parameters OR a single named list of all arguments
# 
# Return Value:
#  (json) = Query ready to pass to mongolite functions
mongoQuery <- function(...) {
  require(jsonlite)
  # Get all args as a list:
  args <- list(...)
  # If singular list argument, elevate that up to args
  if((length(args)==1) && (class(args) == "list") && is.null(names(args))) {
    args <- args[[1]]
  }
  # Drop members that are NULL or have length 0
  for(key in names(args)) {
    if(is.null(args[[key]]) || (length(args[[key]])==0)) {
      args[[key]] <- NULL
    }
  }
  # If args is length 0, return '{}'
  if(length(args)==0) {
    null_query <- '{}'
    class(null_query) <- "json"
    return(null_query)
  }
  # Loop through args and convert any vector args to $in construct
  for(key in names(args)) {
    if((length(args[[key]]) > 1) && (class(args[[key]]) != "list")) {
      args[[key]] <- list('$in'=args[[key]])
    }
  }
  # Convert to JSON and return
  return(toJSON(args, auto_unbox = T))
}

# findIDs - Generic function to get IDs matching any query from any collection
# Given an open DB connection or connection parameters, run a query and get all distinct values of id_field
#
# Parameters:
# DB (mongo object) = If specified, use this open connection to a collection, ignore db_host, db_name, collection
# db_host (character) = If DB is NULL, this specifies the host to open a connection to
# db_name (character) = If DB is NULL, this specifies the DB name to connect to
# collection (character) = If DB is NULL, this specifies the collection to connect to
# id_field (character) = What field in the DB to return IDs from, defaults to the generic mongo _id field
# debug (logical) = Whether to print debug messages, default: FALSE, overridden by options(debug=...)
# ... = All additional params passed to mongoQuery to form the query string
#
# Return Value:
# (character vector) = The ID values returned by DB$distinct call
findIDs <- function(DB=NULL, db_host=NULL, db_name=NULL, collection=NULL, id_field="_id", debug=getOption("debug", default=FALSE), ...) {
  # Open DB connection if an open connection object was not provided
  if(is.null(DB)) {
    stopifnot(!is.null(db_host))
    stopifnot(!is.null(db_name))
    if(debug) {cat("Opening new connection to ", db_host, "/", db_name, ".", collection, "\n", sep="")}
    DB <- openMongo(host=db_host, db=db_name, collection=collection)
  }
  # Construct the query
  my_query <- mongoQuery(...)
  # Call DB$distinct and return the results
  my_ids <- DB$distinct(key=id_field, query=my_query)
  if((class(my_ids)=="list") && (length(my_ids) == 0)) {
    my_ids <- character(0)
  }
  return(my_ids)
}

# as.mongo.date - Convert "POSIXct" to a list object with one member: '$date'=(ISODate str)
#
# Given an R datetime object, e.g. as returned by Sys.time(), convert to a structure appropriate for conversion to JSON and mongo insert
#
# Parameters:
# x (POSIXct) = R datetime object, e.g. as returned by Sys.time()
#
# Return Value:
# (list) = Has one member: '$date'=(ISODate str), converting this to JSON with auto_unbox=T will insert this as a proper date
as.mongo.date <- function(x) {
  list('$date'=strftime(x, "%Y-%m-%dT%H:%M:%SZ", 'UTC'))
}

