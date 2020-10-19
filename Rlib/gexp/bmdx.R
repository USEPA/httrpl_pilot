# Functions to help run BMD Express and process outputs

# --- Functions to create BMD Express config file for command-line --- #

# See documentation of the config parameters and format here: https://github.com/auerbachs/BMDExpress-2/wiki/Command-Line

# expressionDataConfigs - Creates a list object with the expressionDataConfigs section of a BMD Express config file
# 
# Parameters:
# inputFileName (character) = The file or directory with files to import, required
# outputName (character) = If not specified, the imported filename will be used, defaults to NULL
# platform (character) = The name of the platform, defaults to "EPA_BSP_WholeTranscriptome_170424"
# logTransformation (character) = Whether the data is log transformed, must be one of NONE|BASE2|BASE10|NATURAL, defaults to BASE2
#
# Return Value:
# (list) = Contains each of the named parameters above when specified
expressionDataConfigs <- function(
  inputFileName, 
  outputName=NULL, 
  platform="EPA_BSP_WholeTranscriptome_170424", 
  logTransformation="BASE2"
) {
  # TO DO: Make sure all parameters are single length?
  # Warn if logTransformation set to something weird
  # TO DO: Could also check platform
  if(!(logTransformation %in% c("NONE","BASE2","BASE10","NATURAL"))) {
    warning("logTransformation : ", logTransformation, " is not a known valid value for this parameter.\n")
  }
  config_list <- list(inputFileName=inputFileName)
  if(!is.null(outputName)) {
    config_list$outputName = outputName
  }
  config_list$platform = platform
  config_list$logTransformation = logTransformation
  return(config_list)
}

# preFilterConfigs - Creates a list object with the preFilterConfgis section of a BMD Express config file
# Non-required fields can be excluded
#
# Parameters:
# type (character) = The type of pre-filtering module to use, specified in config list as "@type", must be one of anova|williams|oriogen, defaults to williams
# inputName (character) = Which expression data set to pre-filter, optional
# outputName (character) = Name for the output data set, optional
# pValueCutoff (numeric) = P-value cut-off for whichever test is used, defaults to 1 (no filtering by test)
# useMultipleTestingCorrection (logical) = Whether to apply multiple testing correction to the p-values, defaults to FALSE
# filterOutControlGenes (logical) = Whether to filter out control genes, defaults to FALSE
# useFoldChange (logical) = Whether to also filter by fold-change, defaults to TRUE
# foldChange (numeric) = Minimum fold-change to use if useFoldChange=TRUE, defaults to 2, note this is NOT on log scale even when input data is
# numberOfPermutations (integer) = Number of permutations to use for williams test, defaults to 1000
# initialBootstraps (integer) = Number of initial bootstraps to use for oriogen test, required if type=oriogen
# maxBootstraps (integer) = Number of max bootstraps to use for oriogen test, required if type=oriogen
# s0Adjustment (numeric) = s0 adjustment to use for oriogen test, required if type=oriogen
#
# Return Value:
# (list) = Contains each of the named parameters above when specified
preFilterConfigs <- function(
  type = "williams",
  inputName = NULL,
  outputName = NULL,
  pValueCutoff = 1,
  useMultipleTestingCorrection = FALSE,
  filterOutControlGenes = FALSE,
  useFoldChange = TRUE,
  foldChange = 2,
  numberOfPermutations = 1000,
  initialBootstraps = NULL,
  maxBootstraps = NULL,
  s0Adjustment = NULL
) {
  # TO DO: Make sure all parameters are single length?
  # Make sure type is valid
  if(!(type %in% c("anova","williams","oriogen"))) {
    warning("type : ", type, " is not a known valid value for this parameter.\n")
  }
  config_list <- list('@type' = type)
  if(!is.null(inputName)) {
    config_list$inputName <- inputName
  }
  if(!is.null(outputName)) {
    config_list$outputName <- outputName
  }
  config_list$pValueCutoff <- pValueCutoff
  config_list$useMultipleTestingCorrection <- useMultipleTestingCorrection
  config_list$filterOutControlGenes <- filterOutControlGenes
  config_list$useFoldChange <- useFoldChange
  config_list$foldChange <- foldChange
  if(type == "williams") {
    config_list$numberOfPermutations <- numberOfPermutations
  }
  if(type == "oriogen") {
    config_list$initialBootstraps <- initialBootstraps
    config_list$maxBootstraps <- maxBootstraps
    config_list$s0Adjustment <- s0Adjustment
  }
  return(config_list)
}

# modelConfigs - Creates a list object with the bmdsConfigs.modelConfgis section of a BMD Express config file
# Non-required fields can be excluded
# Can short-hand the exp and poly models, e.g. type="poly2" will automatically set @type="poly", degree=2
#
# Parameters:
# type (character) = The type of bmds model to use, specified in config list as "@type", must be one of poly|power|exp|hill, required
# expModel (integer) = The type of exponential model to use, must be one of 2|3|4|5, required when type="exp", or shorthand with type="exp#"
# degree (integer) = The degree of polynomial model to use, must be one of 1|2|3|4, required when type="poly", or shorthand with type="poly#"
#
# Return Value:
# (list) = Contains each of the named parameters above when specified
modelConfigs <- function(
  type,
  expModel = NULL,
  degree = NULL
) {
  # TO DO: Make sure all parameters are single length?
  # Check for short-hand notations: exp# or poly#
  if(grepl("^exp[0-9]$", type)) {
    expModel <- as.integer(sub("^exp", "", type))
    type <- "exp"
  } else if(grepl("^poly[0-9]$", type)) {
    degree <- as.integer(sub("^poly", "", type))
    type <- "poly"
  }
  # Make sure type is valid
  if(!(type %in% c("poly","power","exp","hill"))) {
    warning("type : ", type, " is not a known valid value for this parameter.\n")
  }
  config_list <- list('@type' = type)
  # Set expModel or degree based on type
  if(type == "poly") {
    if(is.null(degree)) {
      warning("No degree specified for poly model.\n")
    } else if(!(degree %in% 1:4)) {
      warning("degree : ", degree, " is not a known valid value for this parameter.\n")
    }
    config_list$degree <- degree
  }
  if(type == "exp") {
    if(is.null(expModel)) {
      warning("No expModel specified for exp model.\n")
    } else if(!(expModel %in% 2:5)) {
      warning("expModel : ", expModel, " is not a known valid value for this parameter.\n")
    }
    config_list$expModel <- expModel
  }
  return(config_list)
}

# bmdsBestModelSelection - Creates a list object with the bmdsConfgis.bmdsBestModelSelection section of a BMD Express config file
#
# Parameters:
# bestPolyTest (integer) = Whether to use Nested Chi Square (1) or Lowest AIC (2) for selecting best poly test, defaults to 2
# pValueCutoff (numeric) = P-value cut-off for goodness of fit test, defaults to 0.05
# flagHillWithKParameter (logical) = Whether to flag hill models with low k values, defaults to TRUE
# kParameterValue (integer) = How to determine cut-off for flagging low k values in hill models, must be 1|2|3, defaults to 2
# bestModelSelectionWithFlaggedHill (integer) = What to do when best model is hill, but flagged for low k, must be 1|2|3|4|5, defaults to 1 (include flagged hill)
# modifyFlaggedHillWithFractionMinBMD (numeric) = How to modify flagged models, as a fraction of min BMD, defaults to 0.05
# 
# Return Value:
# (list) = Contains each of the named parameters above when specified
bmdsBestModelSelection <- function(
  bestPolyTest = 2,
  pValueCutoff = 0.05,
  flagHillWithKParameter = TRUE,
  kParameterValue = 2,
  bestModelSelectionWithFlaggedHill = 1,
  modifyFlaggedHillWithFractionMinBMD = 0.05
) {
  config_list <- list()
  # TO DO: Make sure all parameters are single length?
  # TO DO: Perform type checking?
  if(!(bestPolyTest %in% 1:2)) {
    warning("bestPolyTest : ", bestPolyTest, " is not a known valid value for this parameter.\n")
  }
  config_list$bestPolyTest <- bestPolyTest
  config_list$pValueCutoff <- pValueCutoff
  config_list$flagHillWithKParameter <- flagHillWithKParameter
  if(!(kParameterValue %in% 1:3)) {
    warning("kParameterValue : ", kParameterValue, " is not a known valid value for this parameter.\n")
  }
  config_list$kParameterValue <- kParameterValue
  if(!(bestModelSelectionWithFlaggedHill %in% 1:5)) {
    warning("bestModelSelectionWithFlaggedHill : ", bestModelSelectionWithFlaggedHill, " is not a known valid value for this parameter.\n")
  }
  config_list$bestModelSelectionWithFlaggedHill <- bestModelSelectionWithFlaggedHill
  config_list$modifyFlaggedHillWithFractionMinBMD <- modifyFlaggedHillWithFractionMinBMD
  return(config_list)
}

# bmdsInputConfig - Creates a list object with the bmdsConfgis.bmdsInputConfig section of a BMD Express config file
#
# Parameters:
# maxIterations (integer) = Number of iterations for each model, defaults to 250
# confidenceLevel (numeric) = Confidence level for each model, defaults to 0.95
# constantVariance (logical) = Whether to assume constant variance across all doses, defaults to TRUE
# restrictPower (logical) = Whether to restrict power models, defaults to TRUE
# bmrFactor (numeric) = Multiplier to compute BMR from standard deviation of controls, defaults to 1.349
# 
# Return Value:
# (list) = Contains each of the named parameters above when specified
bmdsInputConfig <- function(
  maxIterations = 250,
  confidenceLevel = 0.95,
  constantVariance = TRUE,
  restrictPower = TRUE,
  bmrFactor = 1.349
) {
  config_list <- list()
  # TO DO: Make sure all parameters are single length?
  # TO DO: Perform type and bounds checking?
  config_list$maxIterations <- maxIterations
  config_list$confidenceLevel <- confidenceLevel
  config_list$constantVariance <- constantVariance
  config_list$restrictPower <- restrictPower
  config_list$bmrFactor <- bmrFactor
  return(config_list)
}

# bmdsConfigs - Creates a list object with the bmdsConfgis section of a BMD Express config file
# modelConfigs can be given as a vector of short-hand model notations, each one will be passed to modelConfigs function above
#
# Parameters:
# modelConfigs (character vector or list) = Models to include under bmdsConfigs.modelConfigs
# bmdsBestModelSelection (list) = bmdsConfigs.bmdsBestModelSelection, if not specified creates default object with specific parameters from ...
# bmdsInputConfig (list) = bmdsConfigs.bmdsInputConfig, if not specified creates default object with specific parameters from ...
# inputCategory (character) = Which input category to use, must be one of anova|williams|oriogen|expression, defaults to williams
# inputName (character) = Name of input data set to use for BMDS modeling, optional
# outputName (character) = Name for BMDS model outputs, optional
# killTime (integer) = Number of seconds to allow models to run before destorying them, optional (See BMD Express for internal default)
# tmpFolder (character) = Override the temp folder used to read and write bmds model files, useful for parallelization of multiple BMD Express runs, optional
# numberOfThreads (integer) = Number of threads to run at once, defaults to 1
# ... = Any other parameters will be passed to sub-sections, should be named as (subsection).(parameter), e.g. bmdsBestModelSelection.bestPolyTest
# 
# Return Value:
# (list) = Contains each of the named parameters above when specified
bmdsConfigs <- function(
  modelConfigs = c("hill", "exp2", "exp3", "exp4", "exp5", "poly1", "poly2", "power"),
  bmdsBestModelSelection = NULL,
  bmdsInputConfig = NULL,
  inputCategory = "williams",
  inputName = NULL,
  outputName = NULL,
  killTime = NULL,
  tmpFolder = NULL,
  numberOfThreads = 1,
  ...
) {
  # Get ... in a mutable list
  ext_args <- list(...)
  # Start the return object
  config_list <- list()
  # TO DO: Make sure all parameters are single length?
  # TO DO: Perform type and bounds checking?
  # If modelConfigs is a character vector, pass each one to modelConfigs function:
  if(!is.list(modelConfigs)) {
    config_list$modelConfigs <- lapply(modelConfigs, modelConfigs)
  } else {
    # Otherwise just pass directly
    config_list$modelConfigs <- modelConfigs
  }
  # If bmdsBestModelSelection is NULL, get any relevant ext_args and pass to bmdsBestModelSelection function
  if(is.null(bmdsBestModelSelection)) {
    bms_args <- ext_args[grepl("^bmdsBestModelSelection[.]", names(ext_args))]
    # Remove the relevant args from overall list - want to warn if anything left in ext_args at the end
    ext_args <- ext_args[which(!(names(ext_args) %in% names(bms_args)))]
    # Drop the prefix before passing bms_args to bmdsBestModelSelection function
    names(bms_args) <- sub("^bmdsBestModelSelection[.]", "", names(bms_args))
    config_list$bmdsBestModelSelection <- do.call("bmdsBestModelSelection", bms_args)
  } else {
    # If specified, just pass through
    config_list$bmdsBestModelSelection <- bmdsBestModelSelection
  }
  # If bmdsInputConfig is NULL, do the same thing as for bmdsBestModelSelection above
  if(is.null(bmdsInputConfig)) {
    bic_args <- ext_args[grepl("^bmdsInputConfig[.]", names(ext_args))]
    # Remove the relevant args from overall list - want to warn if anything left in ext_args at the end
    ext_args <- ext_args[which(!(names(ext_args) %in% names(bic_args)))]
    # Drop the prefix before passing bic_args to bmdsInputConfig function
    names(bic_args) <- sub("^bmdsInputConfig[.]", "", names(bic_args))
    config_list$bmdsInputConfig <- do.call("bmdsInputConfig", bic_args)
  } else {
    # If specified, just pass through
    config_list$bmdsInputConfig <- bmdsInputConfig
  }
  # All remaining parameters are top-level, single value:
  config_list$inputCategory <- inputCategory
  if(!is.null(inputName)) {
    config_list$inputName <- inputName
  }
  if(!is.null(outputName)) {
    config_list$outputName <- outputName
  }
  if(!is.null(killTime)) {
    config_list$killTime <- killTime
  }
  if(!is.null(tmpFolder)) {
    config_list$tmpFolder <- tmpFolder
  }
  config_list$numberOfThreads <- numberOfThreads
  # Check if anything left in ext_args
  if(length(ext_args) > 0) {
    warning(length(ext_args), " parameters were ignored: ", paste(names(ext_args), collapse=", "), "\n")
  }
  return(config_list)
}

# categoryAnalysisConfigs - Creates a list object with the categoryAnalysisConfigs section of a BMD Express config file
#
# Parameters:
# type (character) = Must be one of defined|go|pathway, required
# inputName (character) = Name of BMDS output data to analyze, optional, if not specified then all bmd analysis is run
# outputName (character) = Name of category analysis output data, optional, if not specified then default name is assigned
# removePromiscuousProbes (logical) = Exclude probes mapping to multiple genes? Defaults to TRUE
# removeBMDGreaterHighDose (logical) = Exclude BMD values above the highest dose? Defaults to TRUE
# identifyConflictingProbeSets (logical) = Identify conflicting probesets? Defaults to FALSE (TO DO: Check NTP Guidelines and GUI doc to figure out what this is?)
# correlationCutoffForConflictingProbeSets (numeric) = Correlation cutoff to determine conflicting probe sets, defaults to 1
# bmdPValueCutoff (numeric) = P-value cutoff for BMDs, optional
# bmdBMDLRatioMin (numeric) = Maximum ratio of BMD to BMDL, optional
# bmduBMDRatioMin (numeric) = Maximum ratio of BMD to BMDU, optional
# bmduBMDLRatioMin (numeric) = Maximum ratio of BMDL to BMDU, optional, note this one is not documented on the Wiki but appears in example config files
# nFoldBelowLowestDose (numeric) = optional
# maxFoldChange (numeric) = optional
# prefilterPValueMin (numeric) = optional
# probeFilePath (character) = File path to probe to gene mapping, required when type="defined"
# categoryFilePath (character) = File path to the gene to category mapping, required when type="defined"
# deduplicateGeneSets (logical) = optional
# goCategory (character) = Main GO category to use, must be one of universal|biological_process|molecular_function|cellular components, required when type="go", defaults to "universal"
# signalingPathway (character) Main pathway category to use, must be one of REACTOME, required when type="pathway", defaults to "REACTOME"
# 
# Return Value:
# (list) = Contains each of the named parameters above when specified
categoryAnalysisConfigs <- function(
  type,
  inputName = NULL,
  outputName = NULL,
  removePromiscuousProbes = TRUE,
  removeBMDGreaterHighDose = TRUE,
  identifyConflictingProbeSets = FALSE,
  correlationCutoffForConflictingProbeSets = 1,
  bmdPValueCutoff = NULL,
  bmdBMDLRatioMin = NULL,
  bmduBMDRatioMin = NULL,
  bmduBMDLRatioMin = NULL,
  nFoldBelowLowestDose = NULL,
  maxFoldChange = NULL,
  prefilterPValueMin = NULL,
  probeFilePath = NULL,
  categoryFilePath = NULL,
  deduplicateGeneSets = NULL,
  goCategory = "universal",
  signalingPathway = "REACTOME"
) {
  config_list <- list()
  # Make sure type is valid
  if(!(type %in% c("defined","go","pathway"))) {
    warning("type : ", type, " is not a known valid value for this parameter.\n")
  }
  config_list[["@type"]] <- type
  if(!is.null(inputName)) {
    config_list$inputName <- inputName
  }
  if(!is.null(outputName)) {
    config_list$outputName <- outputName
  }
  config_list$removePromiscuousProbes <- removePromiscuousProbes
  config_list$removeBMDGreaterHighDose <- removeBMDGreaterHighDose
  config_list$identifyConflictingProbeSets <- identifyConflictingProbeSets
  config_list$correlationCutoffForConflictingProbeSets <- correlationCutoffForConflictingProbeSets
  if(!is.null(bmdPValueCutoff)) {
    config_list$bmdPValueCutoff <- bmdPValueCutoff
  }
  if(!is.null(bmdBMDLRatioMin)) {
    config_list$bmdBMDLRatioMin <- bmdBMDLRatioMin
  }
  if(!is.null(bmduBMDRatioMin)) {
    config_list$bmduBMDRatioMin <- bmduBMDRatioMin
  }
  if(!is.null(bmduBMDLRatioMin)) {
    config_list$bmduBMDLRatioMin <- bmduBMDLRatioMin
  }
  if(!is.null(nFoldBelowLowestDose)) {
    config_list$nFoldBelowLowestDose <- nFoldBelowLowestDose
  }
  if(!is.null(maxFoldChange)) {
    config_list$maxFoldChange <- maxFoldChange
  }
  if(!is.null(prefilterPValueMin)) {
    config_list$prefilterPValueMin <- prefilterPValueMin
  }
  if(type == "defined") {
    if(is.null(probeFilePath)) {
      warning("probeFilePath parameter is required when categoryAnalysisConfigs.type = defined.\n")
    }
    config_list$probeFilePath <- probeFilePath
    if(is.null(categoryFilePath)) {
      warning("categoryFilePath parameter is required when categoryAnalysisConfigs.type = defined.\n")
    }
    config_list$categoryFilePath <- categoryFilePath
  }
  if(!is.null(deduplicateGeneSets)) {
    config_list$deduplicateGeneSets <- deduplicateGeneSets
  }
  if(type == "go") {
    if(is.null(goCategory)) {
      warning("goCategory parameter is required when categoryAnalysisConfigs.type = go.\n")
    }
    config_list$goCategory <- goCategory
  }
  if(type == "pathway") {
    if(is.null(signalingPathway)) {
      warning("signalingPathway parameter is required when categoryAnalysisConfigs.")
    }
    config_list$signalingPathway <- signalingPathway
  }
  return(config_list)
}

# bmdExpressConfig - Creates a list object with all the data necessary to write a BMD Express config file
# Each subsection of the config file can be created separately using the subsection functions, 
# or a default object will be created with any specific parameters specified by (subsection).(parameter)
#
# Parameters:
# bm2FileName (character) = The file name of the bm2 file that will be outputted, optional
# jsonExportFileName (character) = The file name of the json export file, optional
# overwrite (logical) = Whether to overwrite existing output, optional
# basePath (character) = override the default base path for storing meta data and annotation files for bmdexpress2. The default location is in the user home directory under "bmdexpress2." This option allows localization of file input/output for running many instances of the application in massive parallel on many compute nodes. Optional.
# expressionDataConfigs (list) = Unnamed list of expressionDataConfigs sections, as generated by expressionDataConfigs function, optional
# preFilterConfigs (list) = Unnamed list of preFilterConfigs sections, as generated by preFilterConfigs function, optional
# bmdsConfigs (list) = Unnamed list of bmdsConfigs sections, as generated by bmdsConfigs function, optional
# categoryAnalysisConfigs (list) = Unnamed list of categoryAnalysisConfigs sections, as generated by categoryAnalysisConfigs function, optional
# global_path (character) = Top-level path that all input/output files should be under, if specified will be appended to all relevant path parameters
# input_dir (character) = If global_path specified, use this subdir for input files, default = "data"
# output_dir (character) = If global_path specified, use this subdir for bm2/json output files, default = "bmde_bm2"
# base_dir (character) = If global_path specified, use this subdir for basePath setting, default = "base"
# sample_name (character) = If specified, will be used to infer the input file, auto-generate the input/output names of each module, and create a sample-specific subdir under basePath, optional
# ... = Any other parameters will be passed to sub-sections, should be named as (subsection).(parameter), e.g. bmdsBestModelSelection.bestPolyTest
# 
# Return Value:
# (list) = Contains each of the named parameters above when specified
bmdExpressConfig <- function(
  bm2FileName = NULL,
  jsonExportFileName = NULL,
  overwrite = NULL,
  basePath  = NULL,
  expressionDataConfigs = NULL,
  preFilterConfigs = NULL,
  bmdsConfigs = NULL,
  categoryAnalysisConfigs = NULL,
  global_path = NULL,
  input_dir = "data",
  output_dir = "bmde_bm2",
  base_dir = "base",
  sample_name = NULL,
  ...
) {
  # Get ... in a mutable list
  ext_args <- list(...)
  # Start the return object
  config_list <- list()
  # TO DO: Make sure all top-level parameters are single length?
  # TO DO: Perform type and bounds checking?
  if(!is.null(bm2FileName)) {
    config_list$bm2FileName <- bm2FileName
  } else {
    # Attempt to auto-fill based on other parameters
    if(!is.null(sample_name)) {
      config_list$bm2FileName <- paste0(sample_name, ".bm2")
    }
  }
  # If bm2FileName and global_path are both defined, append global_path to front with optional output_dir
  if((!is.null(config_list$bm2FileName)) && (!is.null(global_path))) {
    if(!is.null(output_dir)) {
      config_list$bm2FileName <- file.path(output_dir, config_list$bm2FileName)
    }
    config_list$bm2FileName <- file.path(global_path, config_list$bm2FileName)
  }
  # Set jsonExportFileName
  if(!is.null(jsonExportFileName)) {
    config_list$jsonExportFileName <- jsonExportFileName
  } else {
    # Attempt to auto-fill based on other parameters
    if(!is.null(sample_name)) {
      config_list$jsonExportFileName <- paste0(sample_name, ".json")
    }
  }
  # If jsonExportFileName and global_path are both defined, append global_path to front with optional output_dir
  if((!is.null(config_list$jsonExportFileName)) && (!is.null(global_path))) {
    if(!is.null(output_dir)) {
      config_list$jsonExportFileName <- file.path(output_dir, config_list$jsonExportFileName)
    }
    config_list$jsonExportFileName <- file.path(global_path, config_list$jsonExportFileName)
  }
  # Set overwrite parameter
  if(!is.null(overwrite)) {
    config_list$overwrite <- overwrite
  }
  # Set basePath parameter
  if(!is.null(basePath)) {
    config_list$basePath <- basePath
  } else {
    # Attempt to infer if not specified - if global_path, base_dir, and sample_name are defined, set as global_path/base_dir/sample_name
    if((!is.null(global_path)) && (!is.null(base_dir)) && (!is.null(sample_name))) {
      config_list$basePath <- file.path(global_path, base_dir, sample_name)
    }
  }
  # If expressionDataConfigs is NULL, get any relevant ext_args and pass to expressionDataConfigs function
  if(is.null(expressionDataConfigs)) {
    edc_args <- ext_args[grepl("^expressionDataConfigs[.]", names(ext_args))]
    # Remove the relevant args from overall list - want to warn if anything left in ext_args at the end
    ext_args <- ext_args[which(!(names(ext_args) %in% names(edc_args)))]
    # Drop the prefix before passing edc_args to expressionDataConfigs function
    names(edc_args) <- sub("^expressionDataConfigs[.]", "", names(edc_args))
    # If edc_args is missing inputFileName, attempt to infer from global_path, input_dir, and sample_name
    if(!("inputFileName" %in% names(edc_args))) {
      if(!is.null(sample_name)) {
        edc_args$inputFileName <- paste0(sample_name, ".txt")
        if(!is.null(input_dir)) {
          edc_args$inputFileName <- file.path(input_dir, edc_args$inputFileName)
        }
        if(!is.null(global_path)) {
          edc_args$inputFileName <- file.path(global_path, edc_args$inputFileName)
        }
      }
    }
    # If edc_args is missing outputName, attempt to infer from sample_name
    if(!("outputName" %in% names(edc_args))) {
      if(!is.null(sample_name)) {
        edc_args$outputName <- paste0(sample_name, "-expression1")
      }
    }
    config_list$expressionDataConfigs <- list(do.call("expressionDataConfigs", edc_args))
  } else {
    # If specified, just pass through
    config_list$expressionDataConfigs <- expressionDataConfigs
  }
  # If preFilterConfigs is NULL, do the same thing as for expressionDataConfigs above
  if(is.null(preFilterConfigs)) {
    pfc_args <- ext_args[grepl("^preFilterConfigs[.]", names(ext_args))]
    # Remove the relevant args from overall list - want to warn if anything left in ext_args at the end
    ext_args <- ext_args[which(!(names(ext_args) %in% names(pfc_args)))]
    # Drop the prefix before passing pfc_args to preFilterConfigs function
    names(pfc_args) <- sub("^preFilterConfigs[.]", "", names(pfc_args))
    config_list$preFilterConfigs <- list(do.call("preFilterConfigs", pfc_args))
    # If inputName/outputName are missing, attempt to infer
    if(is.null(config_list$preFilterConfigs[[1]]$inputName)) {
      if((length(config_list$expressionDataConfigs) == 1) && !is.null(config_list$expressionDataConfigs[[1]]$outputName)) {
        config_list$preFilterConfigs[[1]]$inputName <- config_list$expressionDataConfigs[[1]]$outputName
      }
    }
    if(is.null(config_list$preFilterConfigs[[1]]$outputName) && !is.null(config_list$preFilterConfigs[[1]]$inputName)) {
      config_list$preFilterConfigs[[1]]$outputName <- paste(config_list$preFilterConfigs[[1]]$inputName, config_list$preFilterConfigs[[1]][["@type"]], sep="-")
    }
  } else {
    # If specified, just pass through
    config_list$preFilterConfigs <- preFilterConfigs
  }
  # If bmdsConfigs is NULL, do the same thing as for expressionDataConfigs above
  if(is.null(bmdsConfigs)) {
    bc_args <- ext_args[grepl("^bmdsConfigs[.]", names(ext_args))]
    # Remove the relevant args from overall list - want to warn if anything left in ext_args at the end
    ext_args <- ext_args[which(!(names(ext_args) %in% names(bc_args)))]
    # Drop the prefix before passing bc_args to bmdsConfigs function
    names(bc_args) <- sub("^bmdsConfigs[.]", "", names(bc_args))
    config_list$bmdsConfigs <- list(do.call("bmdsConfigs", bc_args))
    # If inputName/outputName were not specified here, attempt to infer
    if(is.null(config_list$bmdsConfigs[[1]]$inputName)) {
      if((length(config_list$preFilterConfigs)==1) && !is.null(config_list$preFilterConfigs[[1]]$outputName)) {
        config_list$bmdsConfigs[[1]]$inputName <- config_list$preFilterConfigs[[1]]$outputName
      }
    }
    if(is.null(config_list$bmdsConfigs[[1]]$outputName) && !is.null(config_list$bmdsConfigs[[1]]$inputName)) {
      config_list$bmdsConfigs[[1]]$outputName <- paste0(config_list$bmdsConfigs[[1]]$inputName, "-bmds")
    }
  } else {
    # If specified, just pass through
    config_list$bmdsConfigs <- bmdsConfigs
  }
  # If categoryAnalysisConfigs is NULL, do the same thing as for expressionDataConfigs above, 
  # but ONLY if at least one relevant parameter specified - otherwise leave it blank
  if(is.null(categoryAnalysisConfigs)) {
    bc_args <- ext_args[grepl("^categoryAnalysisConfigs[.]", names(ext_args))]
    if(length(bc_args) > 0) {
      # Remove the relevant args from overall list - want to warn if anything left in ext_args at the end
      ext_args <- ext_args[which(!(names(ext_args) %in% names(bc_args)))]
      # Drop the prefix before passing bc_args to categoryAnalysisConfigs function
      names(bc_args) <- sub("^categoryAnalysisConfigs[.]", "", names(bc_args))
      config_list$categoryAnalysisConfigs <- list(do.call("categoryAnalysisConfigs", bc_args))
      # If inputName/outputName were not specified here, attempt to infer
      if(is.null(config_list$categoryAnalysisConfigs[[1]]$inputName)) {
        if((length(config_list$bmdsConfigs)==1) && !is.null(config_list$bmdsConfigs[[1]]$outputName)) {
          config_list$categoryAnalysisConfigs[[1]]$inputName <- config_list$bmdsConfigs[[1]]$outputName
        }
      }
      if(is.null(config_list$categoryAnalysisConfigs[[1]]$outputName) && !is.null(config_list$categoryAnalysisConfigs[[1]]$inputName)) {
        config_list$categoryAnalysisConfigs[[1]]$outputName <- paste(config_list$categoryAnalysisConfigs[[1]]$inputName, config_list$categoryAnalysisConfigs[[1]][["@type"]], sep="-")
      }
    }
  } else {
    # If specified, just pass through
    config_list$categoryAnalysisConfigs <- categoryAnalysisConfigs
  }
  # Check if anything left in ext_args
  if(length(ext_args) > 0) {
    warning(length(ext_args), " parameters were ignored: ", paste(names(ext_args), collapse=", "), "\n")
  }
  return(config_list)
}

# bmdExpressInput - Creates a single table in the appropriate format to input into BMDExpress
# This function starts from expression and conc data, and generates all sample-specific input files for BMD Express, including the config JSON file
# Note that expression data should already be normalized and filtered, the processing done in this function is minimal.
#
# Parameters:
# EXPR (data.frame or matrix) = All expression data to input to BMD Express, with probes as rows and samples as columns, should already be normalized and filtered as desired
# concs (numeric vector) = Exposure concentration for each sample, must be in same order as columns in EXPR
# sort_conc (logical) = Should the columns in EXPR be re-ordered to have increasing conc? defaults to TRUE
bmdExpressInput <- function(EXPR, concs, sort_conc = TRUE) {
  # Make sure concs and EXPR have matching dimensions
  stopifnot(ncol(EXPR) == length(concs))
  # If concs have names, make sure they match colnames of EXPR
  if(!is.null(names(concs))) {
    stopifnot(all(names(concs) == colnames(EXPR)))
  }
  # Convert EXPR to data.frame and append concs as first row
  bmdxData <- as.data.frame(EXPR)
  bmdxData <- rbind('Probe.ID'=concs, bmdxData)
  # Optionally re-sort the columns by conc
  if(sort_conc) {
    bmdxData <- bmdxData[,order(as.vector(bmdxData[1,]), decreasing = F)]
  }
  bmdxData <- cbind('Probe.ID'=row.names(bmdxData), bmdxData)
  # Return this table
  return(bmdxData)
}

# bmdExpressSampleSetup - Creates all files for running BMD Express on a single test sample
# This function starts from expression and conc data, and generates all sample-specific input files for BMD Express, including the config JSON file
# Note that expression data should already be normalized and filtered, the processing done in this function is minimal.
#
# Parameters:
# EXPR (data.frame or matrix) = All expression data to input to BMD Express, with probes as rows and samples as columns, should already be normalized and filtered as desired
# concs (numeric vector) = Exposure concentration for each sample, must be in same order as columns in EXPR
# sample_name (character) = Will be used to infer the input file, auto-generate the input/output names of each module, and create a sample-specific subdir under basePath, required
# sort_conc (logical) = Should the columns in EXPR be re-ordered to have increasing conc? defaults to TRUE
# json_dir (character) = Subdirectory (under global_path, if specified) to store config file in JSON format, defaults to json_files
# log_dir (character) = Subdirectory (under global_path, if specified) to store log files from BMDExpress runs, defaults to logs
# rerun (logical) = If output file exists already, will give an error instead of overwriting, defaults to FALSE
# debug (logical) = Whether to report debug messages
# ... = All remaining parameters passed to bmdExpressConfig
#
# Return Value:
# (character) = Path to config file that was written
bmdExpressSampleSetup <- function(EXPR, concs, sample_name, 
                                  sort_conc=TRUE, json_dir="json_files", log_dir="logs", 
                                  rerun=FALSE, debug=getOption("debug", default = FALSE), 
                                  ...
) {
  require(jsonlite)
  
  # First pass sample_name and any other params to generate the config structure:
  my_config <- bmdExpressConfig(sample_name=sample_name, ...)
  
  # Determine the location to save the config file
  json_file <- paste0(sample_name, ".json")
  if(!is.null(json_dir)) {
    json_file <- file.path(json_dir, json_file)
  }
  if("global_path" %in% names(list(...))) {
    global_path <- list(...)$global_path
    json_file <- file.path(global_path, json_file)
    # Also append global_path to log_dir if both defined
    if(!is.null(log_dir)) {
      log_dir <- file.path(global_path, log_dir)
    }
  }
  
  # Location for the input file is in the config object already
  expr_file <- my_config$expressionDataConfigs[[1]]$inputFileName
  
  # Check if either output file exists, and if rerun=FALSE stop there
  if(!rerun) {
    stopifnot(!file.exists(json_file))
    stopifnot(!file.exists(expr_file))
  }
  
  # Create all necessary paths
  bmdx_paths <- c(
    dirname(json_file),
    dirname(expr_file),
    dirname(my_config$bm2FileName),
    dirname(my_config$jsonExportFileName),
    my_config$basePath,
    log_dir
  )
  for(bp in bmdx_paths) {
    if(!dir.exists(bp)) {
      if(debug) {cat("Creating directory:", bp, "\n")}
      dir.create(bp, recursive = T)
    }
  }
  
  # Next combine EXPR and concs into a single table for input into BMD Express
  bmdxData <- bmdExpressInput(EXPR=EXPR, concs=concs, sort_conc=sort_conc)
  # Check if input file exists already
  if(debug) {
    if(file.exists(expr_file)) {
      cat("Overwriting existing file:", expr_file, "\n")
    } else {
      cat("Writing new file:", expr_file, "\n")
    }
  }
  write.table(bmdxData, expr_file, sep="\t", quote=F, row.names=F)
  
  # Write the config data in JSON format
  if(debug) {
    if(file.exists(json_file)) {
      cat("Overwriting existing file:", json_file, "\n")
    } else {
      cat("Writing new file:", json_file, "\n")
    }
  }
  write_json(my_config, json_file, auto_unbox = TRUE, pretty = 1)
  
  # Return the path to the config file
  return(invisible(json_file))
}

# readBMDS - Function to read in the tab-delimited text file of bmds results generated by the BMD-Express Export command
# NOTE: There are cases where this file gets written out with extra white space lines at the top and at the end of each data row, which causes problems with read.table
#       This command should be robust to such problems but has not been tested against all versions of BMDExpress
# NOTE: It might be possible (and better) to read this directly from the JSON file and skip the export step entirely
#       Alternatively, could wrap the export step in an R function to automate that part
#
# Parameters:
# bmds_file (character) = Path to tab-delimited file generated by BMD Express export command
# robust (logical) = Whether to use robust but slower loading procedure that can handle some misformatting, defaults to TRUE
# debug (logical) = Whether to report debug messages
#
# Return Value:
# (data.frame) = Complete table of bmds stage outputs
readBMDS <- function(bmds_file, robust = TRUE, debug = getOption("debug", default = FALSE)) {
  if(robust) {
    require(foreach)
    # Use robust mode to read the table
    # First determine the appropriate number of columns
    bmds_nf <- count.fields(bmds_file, sep="\t", quote="", comment.char="", blank.lines.skip = T)
    # Header line should be the first one with >1 fields (everything before that is a possible header block - BMD Express *should* mark these with a comment character but it doesn't)
    header_line <- min(which(bmds_nf > 1))
    if((header_line > 1) && debug) {
      cat("Detected meta-data block on lines 1 -", header_line-1, "with data table header on line", header_line, "\n")
    }
    if(any(bmds_nf > bmds_nf[header_line]) && debug) {
      cat("Header line has", bmds_nf[header_line], "columns, but data lines contain up to", max(bmds_nf), "columns, extra columns will be excluded.\n")
    }
    bmds_nf <- bmds_nf[header_line]
    bmds_scan <- scan(bmds_file, what = as.list(rep("", times=bmds_nf)), sep = "\t", quote = "", skip = header_line-1, blank.lines.skip = T, comment.char = "", quiet = !debug, flush = T, fill = T, strip.white = T)
    # Use foreach to convert and merge each column, excluding headers
    bmds <- foreach(j=bmds_scan, .combine='cbind') %do% {
      j_data <- type.convert(j[-1], as.is=T)
      j_data <- data.frame(X=j_data, stringsAsFactors = F)
      colnames(j_data) <- j[1]
      return(j_data)
    }
    if(debug) {cat("Final data.frame has", nrow(bmds), "obs. of", ncol(bmds), "variables.\n")}
  } else {
    # When robust = FALSE, use read.table which should be faster
    bmds <- read.table(bmds_file, header=T, sep="\t", quote="", comment.char = "", check.names = F)
  }
  # For the best-fit columns that are supposed to be numeric, convert "none" to NA and then convert to numeric
  fix_cols <- c("Best BMD", "Best BMDL", "Best BMDU", "Best fitPValue", "Best fitLogLikelihood", "Best AIC", "Best adverseDirection", "Best BMD/BMDL", "Best BMDU/BMDL", "Best BMDU/BMD")
  fix_cols <- intersect(colnames(bmds), fix_cols)
  fixed_cols <- c()
  for(j in fix_cols) {
    if((class(bmds[, j]) == "character") && any(bmds[, j] == "none")) {
      none_rows <- which(bmds[, j] == "none")
      bmds[none_rows, j] <- NA
      bmds[, j] <- as.numeric(bmds[, j])
      fixed_cols <- c(fixed_cols, j)
    }
  }
  if("Best adverseDirection" %in% fix_cols) {
    bmds[, "Best adverseDirection"] <- as.integer(bmds[, "Best adverseDirection"])
    fixed_cols <- union(fixed_cols, "Best adverseDirection")
  }
  if(debug) {cat("Converted", length(fixed_cols), "columns to appropriate types after changing 'none' to NA:", paste(fixed_cols, collapse = ", "), "\n")}
  return(bmds)
}

# bmdsAnnotProbes - Add probe-annotations (gene symbols and/or Entrez IDs) to a bmds result table
# 
# Parameters:
# bmds (data.frame) = Table of bmds results exported from BMD Express
# probe (data.frame) = Table of probe annotations extracted from httr database
# bmds_name (character) = Name of column in bmds table to use as probe name
# bmds_entrez (character) = Name of column in bmds table to map the Entrez ID into
# bmds_symbol (character) = Name of column in bmds table to map the genes symbols into
# probe_name (character) = Name of column in probe to match up to bmds_name
# probe_entrez (character) = Name of column in probe to copy entrez IDs out of
# probe_symbol (character) = Name of column in probe to copy gene symbols out of
#
# Return Value:
# (data.frame) bmds with the Entrez and/or Symbol columns filled in from probe
bmdsAnnotProbes <- function(bmds, probe,
                            bmds_name = "Probe ID", 
                            bmds_entrez = "Entrez Gene IDs", 
                            bmds_symbol = "Genes Symbols",
                            probe_name = "probe_name", 
                            probe_entrez = "entrez_id", 
                            probe_symbol = "gene_symbol"
) {
  # Index probe by probe_name
  if(!is.null(probe_name)) {
    stopifnot(sum(duplicated(probe[,probe_name])) == 0)
    stopifnot(sum(is.na(probe[,probe_name])) == 0)
    row.names(probe) <- probe[,probe_name]
  }
  # Map entrez ID
  if(!is.null(bmds_entrez) && !is.null(probe_entrez)) {
    bmds[, bmds_entrez] <- probe[bmds[, bmds_name], probe_entrez]
  }
  # Map gene symbols
  if(!is.null(bmds_symbol) && !is.null(probe_symbol)) {
    bmds[, bmds_symbol] <- probe[bmds[, bmds_name], probe_symbol]
  }
  # Return bmds with annotations filled in
  return(bmds)
}

# bmdsFilterProbes - Function to filter down to only those probes appropriate for pathway analysis according to NTP guidelines
#
# Parameters:
# bmds (data.frame) = Table of bmds results exported from BMD Express and with probe annotations added with bmdsAnnotProbes
# top_dose (numeric) = Highest tested dose, all probes with Best BMD above this dose will be excluded
# annot_cols (character) = One or more column names with probe annotations to be used for mapping to pathways. Only probes with one valid ID in each of these columns will be kept.
# annot_sep (character) = Regular expression for separate character/string used to pack multiple IDs into annotation columns.
# pvalue (numeric) = Best fit P-value cutoff, all probes with p-value below this cutoff will be excluded, defaults to 0.1
# bmdul_ratio (numeric) = Ratio of BMDU/BMDL, all probes with ratios above this value will be excluded, defaults to 40
#
# Return Value:
# (data.frame) = Subset of rows in bmds that pass all filtering criteria and are therefore appropriate for pathway analysis
bmdsFilterProbes <- function(bmds, top_dose,
                             annot_cols = "Genes Symbols", annot_sep = "[/\\+,;]",
                             pvalue = 0.1, bmdul_ratio = 40,
                             debug = getOption("debug", FALSE)
) {
  if(debug) {cat("Starting with", nrow(bmds), "probe BMDS results.\n")}
  
  # Filter out probes missing annotations
  rows_noannot <- rep(FALSE, times=nrow(bmds))
  for(j in annot_cols) {
    rows_noannot <- rows_noannot | is.na(bmds[, j])
  }
  if(debug && any(rows_noannot)) {
    cat("Removing", sum(rows_noannot), "probes without gene annotations.\n")
  }
  rows_annot <- !rows_noannot
  bmds_filtered <- bmds[rows_annot, ]
  
  # Filter out probes with multiple annotations
  rows_multi <- rep(FALSE, times=nrow(bmds_filtered))
  for(j in annot_cols) {
    rows_multi <- rows_multi | grepl(annot_sep, bmds_filtered[, j])
  }
  if(debug && any(rows_multi)) {
    cat("Removing", sum(rows_multi), "probes annotated for multiple genes.\n")
  }
  rows_single <- !rows_multi
  bmds_filtered <- bmds[rows_single, ]
  
  # Filter out probes with no converged p-value
  rows_converged <- which((!is.na(bmds_filtered$`Best BMD`)) & (!is.na(bmds_filtered$`Best BMDL`)) & (!is.na(bmds_filtered$`Best BMDU`)))
  rows_removed <- nrow(bmds_filtered) - length(rows_converged)
  if(debug && (rows_removed > 0)) {
    cat("Removing", rows_removed, "probes without fully converged BMD/U/L in best model.\n")
  }
  bmds_filtered <- bmds_filtered[rows_converged, ]
  
  # NOTE: Top dose needs to be provided per chemical...
  rows_below_top <- which(bmds_filtered$`Best BMD` < top_dose)
  rows_above_top <- nrow(bmds_filtered) - length(rows_below_top)
  if(debug && (rows_above_top > 0)) {
    cat("Removing", rows_above_top, "probes with BMD >= highest dose of", top_dose, "\n")
  }
  bmds_filtered <- bmds_filtered[rows_below_top, ]
  
  # Filter by p-value
  rows_above_pval <- which(bmds_filtered$`Best fitPValue` > pvalue)
  rows_below_pval <- nrow(bmds_filtered) - length(rows_above_pval)
  if(debug && (rows_below_pval > 0)) {
    cat("Removing", rows_below_pval, "probes with p-value <=", pvalue, "\n")
  }
  bmds_filtered <- bmds_filtered[rows_above_pval, ]
  
  # Filter by BMDU/L ratio
  rows_within_bmdul <- which(bmds_filtered$`Best BMDU/BMDL` < bmdul_ratio)
  rows_above_bmdul <- nrow(bmds_filtered) - length(rows_within_bmdul)
  if(debug) {
    cat("Removing", rows_above_bmdul, "probes with BMDU/BMDL >=", bmdul_ratio, "\n")
  }
  bmds_filtered <- bmds_filtered[rows_within_bmdul, ]
  
  if(debug) {cat(nrow(bmds_filtered), "probes remain for pathway analysis.\n")}
  
  return(bmds_filtered)
}
  

# bmdsPathwaySingle - Function to perform category analysis on a single pathway using pre-filtered BMDS table
# 
# Parameters:
# bmds (data.frame) = Table of all bmds results that are valid for pathway analysis (pre-filtered by bmdsFilterProbes)
# genes (character vector) = Vector of gene IDs in the pathway of interest
# gene_col (character) = Name of column in bmds to match gene IDs to
# min_genes (integer) = Minimum number of genes to match in bmds for this pathway analysis to be valid, defaults to 3
# min_covg (numeric) = Minimum fraction of genes in pathway that must be found in bmds in order to be valid, defaults to 0.05
# gene_sep (character) = Separator to use when reporting the found genes, defaults to "|"
# avg_dups (logical) = Should BMDs for probes targeting same gene be averaged together before computing overall pathway median? Defaults to TRUE
#
# Return Value:
# (list) = pathway-level BMD/L/U (all NA if not valid), the number of valid probes/genes, pathway coverage, and found genes
bmdsPathwaySingle <- function(bmds, genes,
                              gene_col = "Genes Symbols", 
                              min_genes = 3, min_covg = 0.05, gene_sep = "|", avg_dups = TRUE
) {
  require(foreach)
  # genes should not have any duplicates
  stopifnot(sum(duplicated(genes)) == 0)
  # First determine how many of the genes are found:
  found_genes <- intersect(genes, bmds[, gene_col])
  found_covg <- length(found_genes) / length(genes)
  if((length(found_genes) < min_genes) || (found_covg < min_covg)) {
    # Return results for invalid pathway
    return(list(
      BMD = as.numeric(NA),
      BMDL = as.numeric(NA),
      BMDU = as.numeric(NA),
      N = length(found_genes),
      Covg = found_covg,
      Genes = paste(sort(found_genes), collapse = gene_sep)
    ))
  }
  # Otherwise, continue by computing median BMD/L/U of overlapping probes
  pathway_bmds <- bmds[bmds[, gene_col] %in% genes, ]
  # Optionally average the BMD/L/U values for any genes captured by multiple probes
  if(avg_dups) {
    # This is not discussed in the NTP guideline report but earlier microarray-based studies did this
    # NOTE: The following code could be put in a separate function and might be faster/easier with tidyverse?
    
    # This version should be faster and give same results:
    multi_genes <- unique(pathway_bmds[duplicated(pathway_bmds[, gene_col]), gene_col])
    single_genes <- setdiff(unique(pathway_bmds[, gene_col]), multi_genes)
    if(length(multi_genes) > 0) {
      single_bmds <- pathway_bmds[pathway_bmds[, gene_col] %in% single_genes, c(gene_col,"Best BMD","Best BMDL","Best BMDU")]
      multi_bmds <- foreach(g = multi_genes, .combine='rbind') %do% {
        x <- apply(pathway_bmds[pathway_bmds[, gene_col] == g, c("Best BMD","Best BMDL","Best BMDU")], 2, mean)
        x <- cbind(gene=g, as.data.frame(t(x)))
        colnames(x)[1] <- gene_col
        return(x)
      }
      pathway_bmds <- rbind(single_bmds, multi_bmds)
    }
    
    # Sanity checks that collapsing worked:
    stopifnot(sum(duplicated(pathway_bmds[, gene_col])) == 0)
    stopifnot(nrow(pathway_bmds) == length(found_genes))
  }
  
  # Compute median of BMD, BMDL, BMDU/BMDL, and return with gene/coverage stats
  pathway_results <- list(
    BMD = median(pathway_bmds[, "Best BMD"]),
    BMDL = median(pathway_bmds[, "Best BMDL"]),
    BMDU = median(pathway_bmds[, "Best BMDU"]),
    N = length(found_genes),
    Covg = found_covg,
    Genes = paste(sort(unique(pathway_bmds[, gene_col])), collapse = gene_sep)
  )
  return(pathway_results)
}

# bmdsPathwayMulti - Function to loop over list of pathways and run function above
# Also compiles the outputs of each analysis and combine results into a single table,
# and optionally filter to only those with minimum number of probes/genes and/or coverage
#
# Parameters:
# bmds (data.frame) = Table of all bmds results that are valid for pathway analysis (pre-filtered by bmdsFilterProbes)
# pathways (list of character vectors) = All pathway data, where names = pathway ID, and each list member is a character vector of gene IDs belonging to that pathway
# filter_invalid (logical) = Whether to remove invalid pathways from the result table, defaults to TRUE
# debug (logical) = Whether to report debug info
# ... = All remaining parameters passed bmdsPathwaySingle
bmdsPathwayMulti <- function(bmds, pathways,
                             filter_invalid = TRUE, debug = getOption("debug", default = FALSE),
                             ...
) {
  if(debug) {
    cat("Performing analysis on", length(pathways), "using valid bmds results for", nrow(bmds), "features.\n")
  }
  pathway_results <- lapply(pathways, function(x){
    bmdsPathwaySingle(bmds, x, ...)
  })
  
  # Convert pathway_results to a data.frame
  pathway_table <- do.call(rbind.data.frame, pathway_results)
  pathway_table <- cbind(Pathway=names(pathway_results), pathway_table)
  if(debug) {
    invalid_pathways <- sum(is.na(pathway_table[, "BMD"]))
    valid_pathways <- sum(!is.na(pathway_table[, "BMD"]))
    cat(valid_pathways, "had sufficient features/covg,", invalid_pathways, "were excluded from analysis.\n")
  }
  
  # Optionally filter null pathways and return
  if(filter_invalid) {
    pathway_table <- pathway_table[!is.na(pathway_table[,"BMD"]), ]
  }
  return(pathway_table)
}
