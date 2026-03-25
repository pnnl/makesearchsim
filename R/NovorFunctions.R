
#'Convert mzML to MGF using MSconvert and a config file
#'
#'\code{mzMLtoMGF} is a wrapper for a command line call to MSConvert which
#'converts the input mzML file to an MGF file with scan titles described in the
#'configFile. \code{msconvert.exe} must be in you PATH.
#'
#'
#'@param mzMLname path to the mzML file (with mzML extension)
#'@param outputDir path to the desired output directory, defaults to the current
#'  working directory.
#'@param configFile MSConvert config file containing MGF scan title format. A
#'  default file is included with the package which can be accessed using
#'  \code{system.file("extdata","convertMGF_config.txt",
#'  package="MakeSearchSim", mustWork = TRUE)}.
#'
#'@return Returns the path to the mgf file which will be writtin in the
#'  \code{outputDir} directory.
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#' @seealso \code{\link{convert_raw_file}}
#'
#'@importFrom assertthat assert_that
#'@export
mzMLtoMGF <- function(mzMLname, outputDir = ".", configFile)
{
  #checking inputs
  assertthat::assert_that(file.exists(mzMLname),
                          msg = paste(mzMLname, "not found."))
  assertthat::assert_that(file.exists(configFile),
                          msg = paste(configFile, "not found."))

  assertthat::assert_that(dir.exists(outputDir),
                          msg = paste0(outputDir,
                                       " for function mzMLtoMGF not found."))

  mgfname <- gsub("mzML", "mgf", mzMLname, fixed = TRUE)

  cmdstring <-
    paste0(
      "msconvert ",
      "\"",
      mzMLname,
      "\"",
      " -o ",
      "\"",
      outputDir,
      "\"",
      " --mgf --outfile ",
      mgfname,
      " -c ",
      "\"",
      configFile,
      "\""
    )

  system(cmdstring)

  outpath <-
    paste(normalizePath(outputDir), basename(mgfname), sep = "\\")

  assertthat::assert_that(file.exists(outpath),
                          msg = paste(outpath, "not found."))

  return(mgfname)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#'Run Novor given an MGF and parameter file
#'
#'\code{runNovor} is a wrapper function to run Novor on the command line. It
#'takes an MGF and Novor parameter file in the current working directory, runs
#'Novor, and returns the name of the Novor result file. The Novor result file
#'is the name of the input MGF file (minus extension) with "_Novor.csv" added as
#'a suffix.
#'
#'Two Novor parameter files are included in the package and can be accessed
#'using \code{system.file("extdata",
#'"NovorParams_fTol002_pTol15_MetOx_noAlk_Novor_params.txt",
#'package="MakeSearchSim", mustWork = TRUE)}.
#'
#'\code{NovorParams_fTol002_pTol15_MetOx_Alk_Novor_params.txt} has a static
#'carbamidomethylation modification and a variable methionine oxidation
#'modification. \code{NovorParams_fTol002_pTol15_MetOx_noAlk_Novor_params.txt}
#'only contains a variable methionine oxidation modification.
#'
#'@param MGF_File name of an MGF file in the working directory
#'@param NovorParams path to a Novor parameter file
#'@param Novor_Path Windows style path to the Novor .bat file within the Novor
#'  installaion folder. Backslashes must be escaped by a second backslash. Ex:
#'  `C:\\\\Novor\\\\win\\\\novor.bat`
#'
#'@return Returns the name of the Novor result file
#'
#'@family Novor functions
#'
#'@author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#'@examples
#'\dontrun{
#'runNovor(MGF_File = "test.mgf",
#'         NovorParams = system.file("extdata",
#'              "NovorParams_fTol002_pTol15_MetOx_noAlk_Novor_params.txt",
#'               package="MakeSearchSim", mustWork = TRUE),
#'         Novor_Path = "C:\\Novor\\win\\novor.bat")
#'}
#'
#'@importFrom assertthat assert_that has_extension
#'@export
runNovor <- function(MGF_File, NovorParams, Novor_Path)
{
  #making sure input files exist
  assertthat::assert_that(file.exists(MGF_File),
                          msg = paste(MGF_File, "not found."))

  assertthat::has_extension(MGF_File, "mgf")

  assertthat::assert_that(file.exists(NovorParams),
                          msg = paste(NovorParams, "not found."))

  assertthat::has_extension(NovorParams, "txt")

  assertthat::assert_that(file.exists(Novor_Path),
                          msg = paste(Novor_Path, "not found."))

  assertthat::has_extension(Novor_Path, "bat")


  #creating output file name and path
  outputname <- gsub(".mgf", "_Novor.csv", MGF_File, fixed = TRUE)

  #C:\\Novor\\win\\novor.bat

  #making command line string and quoting all file paths
  cmdstring <-
    paste0(
      "\"",
      Novor_Path,
      "\"",
      " -p ",
      "\"",
      NovorParams,
      "\"",
      " -o ",
      "\"",
      outputname,
      "\"",
      " -f ",
      "\"",
      MGF_File,
      "\""
    )

  system(cmdstring)

  assertthat::assert_that(file.exists(outputname),
                          msg = paste("Novor Output:",
                                      outputname,
                                      "not found."))

  assertthat::assert_that(file.size(outputname) > 0,
                          msg = paste("Novor Output:",
                                      outputname,
                                      "has size of 0."))

  return(outputname)
}



#' Convert mzML and run Novor using a dataTracker row
#'
#' \code{runNovorSingle} calls \code{\link{mzMLtoMGF}} then
#' \code{\link{runNovor}} using the mzML file in the "SimFileName" column of
#' trackingDF_row and the msConvert and Novor parameter files.
#'
#' @param trackingDF_row dataframe with a single row. Must contain columns
#'   "SimFileName" and
#' @param msConvert_config name of msConvert config file in the working
#'   directory
#' @param novor_config name of the Novor parameter file in the working directory
#' @param mgf_output_dir path to the directory where MGF files should be written
#'   by \code{\link{mzMLtoMGF}} . Defaults to the current working directory.
#'
#' @inheritParams runNovor
#'
#' @return Returns a dataframe with a single row. It is the input trackingDF_row
#'   with extra columns for the MGF filename, the Novor output filename, and the
#'   Novor parameter filename
#'
#' @family Novor Functions
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#' @seealso \code{\link{runNovor}} ; \code{\link{mzMLtoMGF}}
#'
#' @importFrom assertthat assert_that
#' @export
runNovorSingle <- function(trackingDF_row,
                           msConvert_config,
                           novor_config,
                           Novor_Path,
                           mgf_output_dir = ".")
{
  #checking mzml file exists
  assertthat::assert_that(
    file.exists(trackingDF_row$SimFileName),
    msg = paste0("Simulated file: ", trackingDF_row$SimFileName,
                 " not found.")
  )

  #checking the mgf_output_dir exists
  assertthat::assert_that(dir.exists(mgf_output_dir),
                          msg = paste0("mgf_output_dir: ",
                                       mgf_output_dir, " not found."))

  #checking that the config files exist
  assertthat::assert_that(
    file.exists(msConvert_config),
    msg = paste0("MSConvert config file: ",
                 msConvert_config, " not found.")
  )

  assertthat::assert_that(file.exists(novor_config),
                          msg = paste0("Novor config file: ",
                                       novor_config, " not found."))

  #running msConvert
  trackingDF_row$MGF_name <-
    mzMLtoMGF(
      mzMLname = trackingDF_row$SimFileName,
      configFile = msConvert_config,
      outputDir = mgf_output_dir
    )

  #running Novor
  trackingDF_row$Raw_Novor_Filename <-
    runNovor(
      MGF_File = trackingDF_row$MGF_name,
      NovorParams = novor_config,
      Novor_Path = Novor_Path
    )

  #adding info about Novor params
  trackingDF_row$Novor_ParamFile <- novor_config

  #returning trackingDF_row
  return(trackingDF_row)
}





#'Filter Novor results using ppm error, Novor score, and peptide length
#'
#'This is a helper function which filters raw Novor results by user specified
#'cutoffs of ppm error, Novor Score, and peptide length. Peptide length is
#'calculated within the function based on the "cleanseq" column and added as
#'"length" column.
#'
#'@param rawdata dataframe of novor results of the type returned by
#'  \code{importSingleNovor}, at minimum must have columns "cleanseq",
#'  "NovorScore", and "ppmError"
#'@param ppmErrorMin minimum allowed precursor error, default -15
#'@param ppmErrorMax maximum allowed precursor error, default +15
#'@param minNovorScore minimum allowed Novor score, default 70
#'@param minLength minimum allowed peptide length, default 6
#'
#'@return Returns a dataframe containing rows that pass the filter parameters
#'
#'@family Novor Functions
#'
#'@author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#'@export
#'@importFrom assertthat assert_that
filterNovor <- function(rawdata,
                        ppmErrorMin = -15,
                        ppmErrorMax = 15 ,
                        minNovorScore = 70,
                        minLength = 6)
  {

    assertthat::assert_that("cleanseq" %in% colnames(rawdata),
                            msg = "Column 'cleanseq' is missing from the rawdata data frame")
    assertthat::assert_that("NovorScore" %in% colnames(rawdata),
                            msg = "Column 'NovorScore' is missing from the rawdata data frame")
    assertthat::assert_that("ppmError" %in% colnames(rawdata),
                            msg = "Column 'ppmError' is missing from the rawdata data frame")

    #adding length
    rawdata$length <- nchar(rawdata$cleanseq)

    #filtering by user parameters
    filt <- rawdata[rawdata$NovorScore >= minNovorScore &
                      rawdata$length >= minLength &
                      rawdata$ppmError <= ppmErrorMax &
                      rawdata$ppmError >= ppmErrorMin, ]

    return(filt)
  }



#' Extract position, scan number, and original filename from an MGF file
#'
#' \code{SingleMGFIDtoScan} reads an MGF file and returns a dataframe which
#' correlates the position of the scan in the MGF file to the scan number and
#' original file name contained in each title line. It expects the scan titles to
#' be in the form made by MSConvert using this filter: filter="titleMaker
#' <RunId>.<ScanNumber>.<ScanNumber>.<ChargeState> File:"<SourcePath>",
#' NativeID:"<Id>""
#'
#' This function is optimized for speed at the cost of memory by reading the mgf
#' file into memory all at once. Reading line by line ends up being either the
#' same speed or slower, I think it's because read/write speed is limited by the
#' hard drive.
#'
#' @param filepath path to an MGF file
#' @param chunkSize number of lines to read from MGF at a time
#'
#' @return Returns a dataframe with columns representing the scan position
#'  (MGForder), scan number (ScanNum), Original dataset (OriginDataset), and
#'  input filename w/o extension (MGFDataset)
#'
#' @family Novor Functions
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}; Anthony Valente
#'
#' @export
#' @importFrom stringi stri_match
#' @importFrom assertthat assert_that
SingleMGFIDtoScan <- function(filepath, chunkSize = 1e6)
{
  assertthat::assert_that(is.numeric(chunkSize),
                          msg = "'chunkSize' must be numeric")
  # read MGF file with at most chunkSize lines at a time
  con = file(filepath, "r")
  rawTitleLines <- list()
  while (TRUE) {
    chunk <- readLines(con, n = chunkSize)
    if (length(chunk) == 0)
      break

    rawTitleLines <- list(rawTitleLines,
                          chunk[grepl("TITLE", chunk)])
  }
  rawTitleLines <- unlist(rawTitleLines, use.names = FALSE)
  close(con)

  #extracting scan numbers
  scans <- stri_match(rawTitleLines, regex = "scan=(\\d+)")

  #extracting name of original raw file
  OriginDatasets <-
    stri_match(rawTitleLines, regex = "File:\"(.*)\\..*\",")

  #creating dataframe with int representing the i-th scan in the file
  #and the scan number
  index <- cbind.data.frame(
    ScanNum = scans[, 2],
    OriginDataset = OriginDatasets[, 2],
    stringsAsFactors = FALSE
  )

  index$ScanNum <- as.numeric(index$ScanNum)
  index$MGForder <- as.numeric(rownames(index))

  #adding dataset name by removing ".mgf"
  dataset_name <-
    gsub(".mgf", "", basename(filepath), fixed = TRUE)

  index$MGFDataset <- dataset_name

  return(index)
}



#' Import Novor results and add scan numbers
#'
#' \code{importSingleNovor} imports Novor results, cleans up the column names,
#' adds a modification free peptide sequence (cleanseq). The function can
#' accommodate Novor files with or without the standard metadata header, however
#' the column names are expected to be the same as in the raw Novor output. If
#' the Novor output does not have scan numbers in the scanNum column (as
#' happens with MGF files made with MSConvert that have no `SCANS=` line`) then
#' \code{\link{SingleMGFIDtoScan}} is called to add correct scan numbers as well
#' as each scan's origin filename. It assumes modification info is always of the
#' pattern "(ModSymbol)".
#'
#' @param NovorOutputFile path to a raw Novor output file. Can be missing the
#'   metadata header but the column names are expected to be the same as in the
#'   raw Novor output.
#'
#' @param MGF_file path to the MGF file corresponding to the `NovorOutputFile`
#'   (necessary if scan numbers not present in Novor output file). Default
#'   `NULL`
#'
#' @return Returns a dataframe of unfiltered Novor results with scan numbers,
#'   clean peptide sequence, and names of the MGF and raw Novor output file
#'
#' @family Novor Functions
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#' @importFrom assertthat assert_that
#' @importFrom utils read.csv
#' @export
importSingleNovor <- function(NovorOutputFile, MGF_file = NULL)
{
  assertthat::assert_that(file.exists(NovorOutputFile),
                          msg = paste(NovorOutputFile, "not found."))

  first_line <- readLines(NovorOutputFile, n = 1)

  #testing if the Novor header is present in case the output files have been
  #previously processed to remove the header
  if (first_line == "#=================== N-O-V-O-R ===================")
  {
    #reads an individual file and skips the metadata the top
    #adds a column for file name
    novordata <- read.csv(
      NovorOutputFile,
      stringsAsFactors = FALSE,
      header = TRUE,
      comment.char = "",
      skip = 19,
      skipNul = TRUE,
      strip.white = TRUE
    )
  } else {
    #reads an file assuming no metadata header is present
    #adds a column for file name
    novordata <- read.csv(
      NovorOutputFile,
      stringsAsFactors = FALSE,
      header = TRUE,
      comment.char = "",
      skipNul = TRUE,
      strip.white = TRUE
    )
  }



  #removing ".csv" to make dataset name
  dataset_name <- gsub(".csv", "", basename(NovorOutputFile))

  #adding filename as column
  novordata$NovorDataset <- dataset_name

  #removing empty column
  novordata$X <- NULL

  #fixing column names
  names <- colnames(novordata)
  names[names == "ppm.1e6.err..mz.z.."] <- "ppmError"
  names[names == "X..id"] <- "NovorID"
  names[names == "mz.data."] <- "mzPrecursor"
  names[names == "z"] <- "ChargePrecursor"
  names[names == "pepMass.denovo."] <- "pepMass_denovo"
  names[names == "score"] <- "NovorScore"
  colnames(novordata) <- names


  #Only calling SingleMGFIDtoScan if no scan numbers present
  if (length(unique(novordata$scanNum)) == 1) {

    #verifying MGF_file exists
    assertthat::assert_that(!is.null(MGF_file),
      msg = "MGF_file required as no scan numbers are present in Novor output file.")

    assertthat::assert_that(file.exists(MGF_file),
                            msg = paste0(MGF_file,
                                         " not found."))

    #calling SingleMGFIDtoScan to add scan
    index <- SingleMGFIDtoScan(MGF_file)

    #merging index with novordata
    export <- merge(
      novordata,
      index,
      by.x = c("NovorID"),
      by.y = c("MGForder"),
      all.x = TRUE,
      sort = FALSE
    )

    #removing empty scan number column
    export$scanNum <- NULL

  } else {
    #adding MGFDataset and OriginDataset column created by SingleMGFIDtoScan
    export <- novordata
    export$OriginDataset <- dataset_name
    export$MGFDataset <- dataset_name


    #renaming scan number column
    colnames(export)[colnames(export) == "scanNum"] <- "ScanNum"

  }

  #adding clean sequence without modification info
  export$cleanseq <- gsub("\\([^\\)]+\\)", "", export$peptide)

  return(export)
}


#' Import/Filter Novor results using dataframe style input
#'
#' \code{\link{importFilter_Novor}} calls \code{\link{importSingleNovor}} and
#' \code{\link{filterNovor}} using a single-row dataframe of the style created
#' by \code{\link{MakeSimFile}} and \code{\link{runNovorSingle}}. As such it is
#' meant to be used with dlply or \code{\link[plyr]{ddply}} on a dataframe
#' containing input parameters. It imports the Novor results, adds metadata from
#' the MGF file, then filters the results using user specifed parameters.
#'
#' @param trackingDF_row At minimum must have columns "Raw_Novor_Filename" and
#'   "MGF_name"
#' @inheritParams filterNovor
#'
#' @return Returns a dataframe of filtered Novor results with added metadata
#'
#' @family Novor Functions
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#' @importFrom assertthat assert_that
#' @export
importFilter_Novor <- function(trackingDF_row,
                               ppmErrorMin = -15,
                               ppmErrorMax = 15 ,
                               minNovorScore = 70,
                               minLength = 6)
  {
    #making sure we have columns with correct names
    assertthat::assert_that("Raw_Novor_Filename" %in% colnames(trackingDF_row),
                     msg = "Column 'Raw_Novor_Filename' missing from trackingDF_row")

    assertthat::assert_that("MGF_name" %in% colnames(trackingDF_row),
                            msg = "Column 'MGF_Name' missing from trackingDF_row")

    #making sure input files exists
    assertthat::assert_that(
      file.exists(trackingDF_row$Raw_Novor_Filename),
      msg = paste0(
        "Novor output file: ",
        trackingDF_row$Raw_Novor_Filename,
        " not found."
      )
    )

    assertthat::assert_that(
      file.exists(trackingDF_row$MGF_name),
      msg = paste0("Simulated MGF file: ",
                   trackingDF_row$MGF_name,
                   " not found.")
    )

    assertthat::assert_that(is.numeric(ppmErrorMin),
                            msg = "'ppmErrorMin' must be numeric")
    assertthat::assert_that(length(ppmErrorMin) == 1,
                            msg = "'ppmErrorMin' must be of length 1")

    assertthat::assert_that(is.numeric(ppmErrorMax),
                            msg = "'ppmErrorMax' must be numeric")
    assertthat::assert_that(length(ppmErrorMax) == 1,
                            msg = "'ppmErrorMax' must be of length 1")

    assertthat::assert_that(is.numeric(minNovorScore),
                            msg = "'minNovorScore' must be numeric")
    assertthat::assert_that(length(minNovorScore) == 1,
                            msg = "'minNovorScore' must be of length 1")

    assertthat::assert_that(is.numeric(minLength),
                            msg = "'minLength' must be numeric")
    assertthat::assert_that(length(minLength) == 1,
                            msg = "'minLength' must be of length 1")

    #import Novor and adding scan info
    rawNovor <- importSingleNovor(trackingDF_row$Raw_Novor_Filename,
                                  trackingDF_row$MGF_name)

    #filtering results
    filtNovor <- filterNovor(
      rawdata = rawNovor,
      ppmErrorMin = ppmErrorMin,
      ppmErrorMax = ppmErrorMax,
      minNovorScore = minNovorScore,
      minLength = minLength
    )

    return(filtNovor)
  }


#' Precursor error by Novor score hexbin plot
#'
#' \code{plotNovorHexBin} makes a hexbin plot of precursor error by Novor score
#' after filtering using the input parameters. Title is made using first value
#' in MGFDataset column.
#'
#' @param NovorDF a dataframe of Novor PSMs. Must have columns: "ppmError" and
#'   "NovorScore", if not specified dataset name assumed to be in "MGFDataset"
#'   column.
#' @param dataset_column name of the column containing the dataset name, first
#'   entry will be used in the graph title.
#' @param ppmMin Minimum precursor error to plot
#' @param ppmMax Maximum precursor error to plot
#' @param scoreMin Minimum Novor score to plot
#' @return Returns a ggplot2 plot object
#'
#' @family Novor Functions
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#' @export
#' @import ggplot2
#' @import hexbin
#' @importFrom assertthat assert_that
plotNovorHexBin <-
  function(NovorDF,
           dataset_column = "MGFDataset",
           ppmMin = -15,
           ppmMax = 15,
           scoreMin = 20
  )
{

  #checking inputs
  assertthat::assert_that(
    dataset_column %in% colnames(NovorDF),
    msg = paste0("dataset_column: ",
                 dataset_column,
                 " not found in column names.")
  )

  assertthat::assert_that("ppmError" %in% colnames(NovorDF),
                          msg = "ppmError not found in column names.")

  assertthat::assert_that("NovorScore" %in% colnames(NovorDF),
                          msg = "NovorScore not found in column names.")

  assertthat::assert_that(ppmMin < ppmMax,
                          msg = "ppmMin must be less than ppmMax")

  assertthat::assert_that(scoreMin >= 0 &
                            scoreMin <= 100,
                          msg = "scoreMin must be between 0 and 100")

  assertthat::assert_that(is.numeric(ppmMin),
                          msg = "'ppmMin' must be numeric")
  assertthat::assert_that(length(ppmMin) == 1,
                          msg = "'ppmMin' must be of length 1")

  assertthat::assert_that(is.numeric(ppmMax),
                          msg = "'ppmMax' must be numeric")
  assertthat::assert_that(length(ppmMax) == 1,
                          msg = "'ppmMax' must be of length 1")

  assertthat::assert_that(is.numeric(scoreMin),
                          msg = "'scoreMin' must be numeric")
  assertthat::assert_that(length(scoreMin) == 1,
                          msg = "'scoreMin' must be of length 1")


  #filtering dataframe
  filtDF <- NovorDF[NovorDF$ppmError >= ppmMin &
                    NovorDF$ppmError <= ppmMax &
                    NovorDF$NovorScore >= scoreMin,]

  #making plot title
  title <- paste0("Novor Precursor Error vs Score:\n",
                  NovorDF[1, dataset_column])

  #making plot
  graph <- ggplot(data = filtDF, aes(x = filtDF$ppmError, y = filtDF$NovorScore))
  graph <- graph + geom_hex()
  graph <- graph + scale_fill_gradient(low = "gray", high = "red")
  graph <- graph + ggtitle(title) + xlab("Precursor Error (ppm)") + ylab("Novor Score")
  return(graph)
}

#' Save list of ggplot2 objects to single pdf
#'
#' Got this from StackOverflow \url{https://stackoverflow.com/questions/12234248/printing-multiple-ggplots-into-a-single-pdf-multiple-plots-per-page}
#'
#' @param list (list) List of ggplot2 objects.
#' @param filename (chr) What to call the pdf.
#'
#' @return Invisible NULL.
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @export
GG_save_pdf = function(list, filename) {
  #start pdf
  pdf(filename)

  #loop
  for (p in list) {
    print(p)
  }

  #end pdf
  dev.off()

  invisible(NULL)
}

