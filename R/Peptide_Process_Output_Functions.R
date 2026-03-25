#this script contains functions related to processing peptides and writing peptide output files


#' Digest a set of peptides and discard short ones
#'
#'\code{digestPeptides} returns set of peptides with no missed trypsin cut sites
#'created by digesting the input peptides (splitting at R or K) and discarding
#'those with a length less than 'discardLessThan'.
#'
#' @param peptides a character vector containing peptides
#' @param prolineBlock should P after R or K prevent splitting at that site,
#'   defaults to TRUE
#' @param discardLessThan minimum length of output peptides, defaults to 6
#'
#' @seealso \code{\link{makePepFile}}
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#' @return a character vector of peptides
#' @export
digestPeptides <- function(peptides, prolineBlock = TRUE, discardLessThan = 6)
{
  assertthat::assert_that(is.logical(prolineBlock),
                          msg = "'prolineBlock' must be TRUE or FALSE")

  assertthat::assert_that(is.numeric(discardLessThan),
                          msg = "'discardLessThan' must be numeric")

  assertthat::assert_that(length(discardLessThan) == 1,
                          msg = "'discardLessThan' must be of length 1")

  #returns results for a single peptide
  digestSingle <- function(peptide,prolineBlock)
  {
    #strsplit can accept a position between characters as a match
    #anything "included in the result" of the match is discarded
    #so you have to use look ahead and behind
    if (prolineBlock)
    {
      #before must match [RK] and after cannot match [P]
      digested <- strsplit(peptide,'(?<=[RK])(?![P])', perl = TRUE)
    } else
    {
      #before must match [RK]
      digested <- strsplit(peptide,'(?<=[RK])', perl = TRUE)
    }
    return(digested[[1]])
  }

  #calls digestSingle on all input peptides
  digested <- as.vector(unlist(sapply(peptides,digestSingle, prolineBlock)))

  #discards short peptides according to discardLessThan
  digested <- digested[ifelse(nchar(digested) >= discardLessThan, TRUE, FALSE)]

  digested <- unique(digested)
  return(digested)
}


#' Export digested peptides from a PSM dataframe
#'
#' \code{makePepFile} takes a dataframe containing at least a column of peptide
#' sequences and a column containing the dataset name and writes a .csv file
#' containing the unique list of peptides resulting from digesting the peptides
#' with \code{\link{digestPeptides}}. The csv file contains a single column of
#' peptides with "x" as the first entry. The filename is the first entry in the
#' 'datasetCol' column with "_peptides.csv" added. The function returns a
#' dataframe row containing the name of the output file along with the number of
#' peptides contained in the file. You can also select an output directory by
#' supplying a relative file path from the working directory.
#'
#' @param inputdf dataframe containing at least a column of peptide sequences
#'   and a column containing the dataset name
#' @param datasetCol name of the column containing the dataset name, the first
#'   value will be used as the export filename with "_peptides.csv" added
#' @param sequenceCol name of the column containing the clean peptide sequences
#'   (no modification info)
#' @param outputFileDir relative file path to the desired output folder, default
#'   is the current working directory
#'
#' @inheritParams digestPeptides
#'
#' @return a dataframe row containing the name of the output file along with the
#'   number of peptides contained in the file
#'
#' @seealso \code{\link{digestPeptides}}
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#' @export
#' @importFrom assertthat assert_that
#' @importFrom utils write.csv
makePepFile <- function(inputdf,
                        datasetCol,
                        outputFileDir = ".",
                        sequenceCol = "cleanseq",
                        prolineBlock = TRUE,
                        discardLessThan = 6)
{
  assertthat::assert_that(datasetCol %in% colnames(inputdf),
                          msg = "datasetCol not found in inputdf column names")

  assertthat::assert_that(sequenceCol %in% colnames(inputdf),
                          msg = "sequenceCol not found in inputdf column names")

  assertthat::assert_that(dir.exists(outputFileDir),
                          msg = paste("outputFileDir",
                                      outputFileDir,
                                      "not found."))

  assertthat::assert_that(is.logical(prolineBlock),
                          msg = "'prolineBlock' must be TRUE or FALSE")

  assertthat::assert_that(is.numeric(discardLessThan),
                          msg = "'discardLessThan' must be numeric")

  assertthat::assert_that(length(discardLessThan) == 1,
                          msg = "'discardLessThan' must be of length 1")

  torun<-digestPeptides(inputdf[,sequenceCol],
                        prolineBlock=prolineBlock,
                        discardLessThan=discardLessThan)

  assertthat::assert_that(length(torun) > 0,
                      msg = paste0("No digested peptides returned for dataset: ",
                                   inputdf[1, datasetCol]))

  export_filename <- paste0(inputdf[1, datasetCol], "_peptides.csv")

  outputFilePath <-
    paste(normalizePath(outputFileDir), export_filename, sep = "\\")

  write.csv(torun, file = outputFilePath, row.names = FALSE)

  assertthat::assert_that(
    file.exists(outputFilePath),
    msg = paste0("Peptide list file: ",  outputFilePath, " not found.")
  )

  return(data.frame(
    PeptideListFile = export_filename,
    Num_Digested_Peptides = length(torun)
  ))
}


