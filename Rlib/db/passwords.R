# Module for credentials/password management
# This is analogous to lib/db/passwords.py
# This module does not handle any of the mongo connections, it just provides functions to
# find the appropriate user/pass for a DB connection so these don't need to be exposed in the code

# The following global options can be set using the options function before calling 
# any DB connection functions and they will be recognized here:
#  httrDefaultHost = The URL of the mongo host to use
#  httrDefaultUser = The default username to use for Mongo connections
#  httrDefaultPasswd = The default password to use for Mongo connections

# getCredentials()
# This function replaces the Keychain class in lib/db/passwords.py
# Rather than load entire keychains into memory, it just quickly scans ~/.mongopw and ~/.mngdb/passwd for relevant user,passwd combo
# Returns NULL or a list with members user,passwd (and host,db if specified)
getCredentials <- function(host=getOption("httrDefaultHost"), db=NULL) {
  # TO DO: if only one of host,db specified, generate a warning
  if(!is.null(host) & !is.null(db)) {
    # If both host and db were specified, then try searching ~/.mongopw
    # TO DO: Make the name of ~/.mongopw configurable
    keyfile <- "~/.mongopw"
    if(file.exists(keyfile)) {
      # Use jsonlite package to parse the keychain json file
      require(jsonlite)
      # Read in the keyfile - if formatted properly this should become a dataframe with columns host, db, user, passwd
      keydata <- read_json(keyfile, simplifyVector = T)
      # Check for a row matching host and db (if multiple rows, take first one but throw a warning)
      matchRows <- which((keydata[,"host"] == host) & (keydata[,"db"] == db))
      if(length(matchRows) > 1) {
        warning("Found multiple entries in ", keyfile, " with host=", host, " and db=", db, ", using first such entry.")
      }
      if(length(matchRows) > 0) {
        matchRows <- matchRows[1]
        return(as.list(keydata[matchRows,]))
      }
    }
  }
  # At this point, either no host/db specified, keyfile didn't exist, or no matching entries
  # Check for /.mngdb/passwd file which should just have "user:passwd"
  pwfile <- "~/.mngdb/passwd"
  if(file.exists(pwfile)) {
    pwdata <- readLines(pwfile)
    pwdata <- sub("^[ \t]","",pwdata)
    pwdata <- sub("[ \t]$","",pwdata)
    pwdata <- strsplit(pwdata, split=":", fixed=T)[[1]]
    return(list(user=pwdata[1], passwd=pwdata[2]))
  }
  # At this point, defer to global options httrDefaultUser and httrDefaultPasswd if defined
  if(!is.null(getOption("httrDefaultUser")) & !is.null(getOption("httrDefaultPasswd"))) {
    return(list(user=getOption("httrDefaultUser"), passwd=getOption("httrDefaultPasswd")))
  }
  # At this point, no user,passwd was found - throw warning and return NULL
  warning("Could not find any credential files.")
  return(NULL)
}
