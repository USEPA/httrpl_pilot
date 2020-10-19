# Main API functions for R portions of httrpl
# Make sure to source this script with chdir=T so it can find all sub-modules
#
# Usage:
#  source("Rlib/httrpl.R", chdir = T)
#

# Load sub-modules
db_modules <- c("mongo.R", "probe.R", "trtgrp.R", "well.R", "fc.R")
gexp_modules <- c("qc.R", "deseq2.R", "tss.R", "bmdx.R")
modules <- c(paste0("db/", db_modules), paste0("gexp/", gexp_modules) )
for(mod in modules) {
  source(mod, chdir=T, local=T)
}

# Function to run QC functions on all documents in httr_counts and store results in httr_counts_qc
# Also pulls probe info out of httr_probe (should delegate to functions in db/probes.R),
# NOTE: There's no real reason here to do sample_id batching, but could still use config files to store the DB info and collection names?
#
# Parameters:
# db_host (character) = MongoDB host URL
# db_name (character) = MongoDB name
# db_probe (character) = Name of probe collection, default: httr_probe
# db_counts (character) = Name of counts collection, default: httr_counts
# db_collection (character) = Name of counts_qc collection, default: httr_counts_qc
# keep_flags (character vector) = Keep all probes with these flags in the downstream QC (everything else is considered BAD), default: "OK"
# calc_flag (logical) = Whether or not to calculate qc_flag, default True but when set to False all flags will be set to "OK"
# rerun (logical) = Whether to overwrite existing data in httr_counts_qc
# debug (logical) = Whether to report debug messages
qcBatch <- function(db_host, db_name,
  db_probe="httr_probe", db_counts="httr_counts", db_collection="httr_counts_qc", keep_flags="OK", 
  calc_flag=TRUE, rerun=FALSE, debug=getOption("debug", default=FALSE)
) {
  # Get probe manifest
  probeManifest <- getProbeManifest(db_host=db_host, db_name=db_name, collection=db_probe)
  if(debug) cat("[", date(), "] Pulled table of ", nrow(probeManifest), " probes from ", db_probe, "\n", sep="")
  
  # Pull out proben and bad_probes
  bad_probes <- probeManifest[!(probeManifest[,"probe_flag"] %in% keep_flags),"probe_name"]
  proben <- nrow(probeManifest)
  if(debug) cat(length(bad_probes), "out of", proben, "probes will be filtered out as bad.\n")
  
  # Open connection to httr_counts, get all sample_id fields
  httr_counts <- openMongo(host = db_host, db = db_name, collection = db_counts)
  count_samples <- sort(httr_counts$distinct(key="sample_id"))
  if(debug) cat("[", date(), "] ", db_counts, " contains ", length(count_samples), " sample IDs.\n", sep="")
  
  # Open connection to httr_counts_qc, split count_samples by what's in qc collection already
  httr_counts_qc <- openMongo(host = db_host, db = db_name, collection = db_collection)
  if(debug) cat("[", date(), "] ", db_collection, " contains ", httr_counts_qc$count(), " samples before running QC pipeline.\n", sep="")
  count_samples <- splitIDs(DB=httr_counts_qc, ids=count_samples)
  
  # Decide what to do with existing samples based on rerun param
  if(length(count_samples$present) > 0) {
    if(rerun) {
      # Delete the QC data for these sample_id
      warning("Dropping QC data for ", length(count_samples$present), " samples; rerun=TRUE\n")
      deleteByID(DB=httr_counts_qc, ids=count_samples$present, delete_all=FALSE, debug=debug)
      count_samples <- unique(unlist(count_samples))
    } else {
      warning("Skipping ", length(count_samples$present), " samples alread in ", db_collection, "; rerun=FALSE\n")
      count_samples <- count_samples$absent
    }
  } else {
    # Doesn't matter, nothing has been processed yet
    count_samples <- count_samples$absent
  }
  
  # Loop over count_samples and process QC for each one
  cat("Running QC pipeline on", length(count_samples), "samples...\n")
  
  # NOTE: Could use foreach here to parallelize this and return the debug strings from each one to output in linear log
  # However, parallelizing DB I/O probably won't speed things up and may lead to errors, so would need to separate out just the countQC call and parallelize that part
  for(sample_id in count_samples) {
    if(debug) cat("[", date(), "] Running QC on ", sample_id, "\n", sep="")
    counts_doc <- httr_counts$iterate(query=paste0('{"sample_id":"',sample_id,'"}'), fields='{}')$one()
    # Convert probe_cnts to vector
    # TO DO: There should be a standardized access function for pulling a document out of httr_counts that also does this conversion
    counts_doc$probe_cnts <- unlist(counts_doc$probe_cnts)
    # Make sure all probe_cnt in probeManifest
    if(!all(names(counts_doc$probe_cnts) %in% probeManifest[,"probe_name"])) {
      unknown_probes <- setdiff(names(counts_doc$probe_cnts), probeManifest[,"probe_name"])
      warning(sample_id, " contains ", length(unknown_probes), " probe names not in ", db_probe, ": ", paste(head(unknown_probes), collapse=", "), "...\n")
    }
    qc_doc <- countQC(counts_doc=counts_doc, bad_probes=bad_probes, proben=proben, calc_flag=calc_flag, debug=debug)
    httr_counts_qc$insert(qc_doc, auto_unbox=T)
  }
 
  # Make sure httr_counts and httr_counts_qc contain the same number of documents, same sample_ids
  counts_sz <- httr_counts$count()
  qc_sz <- httr_counts_qc$count()
  if(debug) cat("[", date(), "] ", db_collection, " contains ", qc_sz, " documents after QC pipeline.\n", sep="")
  if(counts_sz != qc_sz) {
    warning(db_counts, " has ", counts_sz, " documents, but ", db_collection, " only has ", qc_sz, " documents after running full QC pipeline - these two collections should have equal size.\n")
  }
  count_samples <- sort(httr_counts$distinct(key="sample_id"))
  qc_samples <- sort(httr_counts$distinct(key="sample_id"))
  count_only <- setdiff(count_samples, qc_samples)
  if(length(count_only) > 0) {
    warning(db_counts, " contains ", length(count_only), " samples without matching entry in ", db_collection, " after QC pipeline:", paste(head(count_only), collapse=", "), "...\n")
  }
  qc_only <- setdiff(qc_samples, count_samples)
  if(length(qc_only) > 0) {
    warning(db_collection, " contains ", length(qc_only), " samples without source entry in ", db_counts, " after QC pipeline:", paste(head(qc_only), collapse=", "), "...\n")
  }
}
