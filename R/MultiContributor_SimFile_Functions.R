#This contains functions to make simulations from multiple contributors using a long format SimRecipes table


#' Select scans from a long format SimRecipes row
#'
#' \code{SingleDataset_SelectScans} uses parameters stored in the input
#' \code{recipeRow} and uses them to isolate the correct set of MSGF+ results
#' from \code{dbResultList} and call \code{\link{selectScans}}. It then adds
#' columns for what and how many scans, peptides, and proteins were selected
#' from the MSGF+ results.
#'
#'
#' @param recipeRow single row of a long format SimRecipes table (as a
#'   dataframe) must contain columns "SimID", "Type", "Dataset",
#'   "TargetNumContributorProteins", "TargetContributorPep_per_Prot". All other
#'   columns are passed through to output.
#' @param dbResultList list of dataframes containing unrolled MSGF+ results (one
#'   protein per row) one of which must contain the dataset referenced in
#'   recipeRows. The dataset must match to the name of one of the entries in the
#'   list. dataframe must have columns "cleanseq", "Protein", and "ScanNum"
#'
#' @inheritParams selectScans
#'
#' @return Returns the input row with added columns containing selected scans
#'   and what peptides and proteins they corrispond to.
#'
#' @family Multi-contributor Simulations
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#' @importFrom assertthat assert_that
#' @export
SingleDataset_SelectScans <- function(recipeRow,
                                      dbResultList,
                                      rng_version = "R_version_default")
{
  #making sure recipeRow has all the needed columns
  assertthat::assert_that("SimID" %in% colnames(recipeRow),
              msg = "SimID column missing from recipeRow")

  assertthat::assert_that("Dataset" %in% colnames(recipeRow),
              msg = "Dataset column missing from recipeRow")

  assertthat::assert_that("Type" %in% colnames(recipeRow),
              msg = "Type column missing from recipeRow")

  assertthat::assert_that("TargetNumContributorProteins" %in% colnames(recipeRow),
              msg = "TargetNumContributorProteins column missing from recipeRow")

  assertthat::assert_that("TargetContributorPep_per_Prot" %in% colnames(recipeRow),
              msg = "TargetContributorPep_per_Prot column missing from recipeRow")


  #making sure the type column is valid
  assertthat::assert_that(recipeRow$Type %in% c("Primary", "Contributor"),
              msg = paste(recipeRow$Type,
                          " is an invalid Type. Please put Primary or Contributor in the Type column."))

  #making sure SimID is numeric
  assertthat::assert_that(is.numeric(recipeRow$SimID),
              msg = "SimID is not numeirc. Make sure SimID is an int.")


  #making sure SimID is not negative
  assertthat::assert_that(recipeRow$SimID > 0, msg = "SimID must be > 0.")


  #checking type column, if type = primary, skipping scan selection and returning NA in the columns with info about selected scans and stuff
  if (recipeRow$Type == "Contributor")
  {
    #making sure contributor dataset is in dbResultList
    assertthat::assert_that(recipeRow$Dataset %in% names(dbResultList),
                msg = paste0("Dataset: ",recipeRow$Dataset,
                             " was not found in dbResults list names."))

    #finding database results for contributor dataset
    datasetDBresults <-
      dbResultList[[match(recipeRow$Dataset, names(dbResultList))]]

    #stupid user resistance
    assertthat::assert_that(nrow(datasetDBresults) > 0,
                msg = "Number of rows in datasetDBresults must be > 0.")
    assertthat::assert_that(recipeRow$TargetNumContributorProteins > 0,
                msg = "TargetNumContributorProteins must be > 0.")
    assertthat::assert_that(recipeRow$TargetContributorPep_per_Prot > 0,
                msg = "contributorPep_per_Prot must be > 0.")


    #calling selectScans
    ContributorScans <- selectScans(databaseResults = datasetDBresults,
                                  Nprot = recipeRow$TargetNumContributorProteins,
                                  Mpep = recipeRow$TargetContributorPep_per_Prot,
                                  SimID = recipeRow$SimID,
                                  rng_version = rng_version)

    assertthat::assert_that(length(ContributorScans) > 0,
                msg = "selectScans returned no scans")

    #finding what peptides were selected
    ContributorPeptides <- unique(datasetDBresults$cleanseq
                                [datasetDBresults$ScanNum %in%
                                    ContributorScans])

    #finding what proteins were selected
    ContributorProteins <- unique(datasetDBresults$Protein
                                [datasetDBresults$ScanNum %in%
                                    ContributorScans])

    #adding contributor scan information
    recipeRow$ContributorScans <- list(ContributorScans)
    recipeRow$NumContributorScans <- length(ContributorScans)

    #adding contributor peptide information
    recipeRow$ContributorPeptides <- list(ContributorPeptides)
    recipeRow$NumContributorPeptides <- length(ContributorPeptides)

    #adding contributor protein information
    recipeRow$ContributorProteins <- list(ContributorProteins)
    recipeRow$NumContributorProteins <- length(ContributorProteins)

  } else
  {
    #adding contributor scan information
    recipeRow$ContributorScans <- NA
    recipeRow$NumContributorScans <- NA

    #adding contributor peptide information
    recipeRow$ContributorPeptides <- NA
    recipeRow$NumContributorPeptides <- NA

    #adding contributor protein information
    recipeRow$ContributorProteins <- NA
    recipeRow$NumContributorProteins <- NA

  }

  return(recipeRow)
}




#' Use long format dataTracker row to get Novor results
#'
#' \code{subsetSingleNovor} takes single row of a long format dataTracker from
#' the multiple dataset / simulation pipeline and returns the appropriate set of
#' Novor results. It returns all the results for the dataset if Type = Primary
#' and only the subset according to ContributorScans if the Type = Contributor.
#' Adds a column for SimID and for the dataset of origin(OriginDataset) to the
#' exported dataframe.
#'
#' @param dataTracker_row a row of a dataTracker dataframe produced by
#'   \code{\link{SingleDataset_SelectScans}} must contain at least columns
#'   "SimID", "Dataset", "Type", "ContributorScans"
#'
#' @inheritParams MixMultipleNovor
#'
#' @return A dataframe containing the appropriate set of Novor results for the
#'   dataset referenced in dataTracker_row
#'
#' @family Multi-contributor Simulations
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#' @importFrom assertthat assert_that
#' @export
subsetSingleNovor <- function(dataTracker_row, RawNovorList)
{

  #making sure Dataset is the name of one of the entries in RawNovorList
  assertthat::assert_that(dataTracker_row$Dataset %in% names(RawNovorList),
              msg = paste0("Dataset: ",dataTracker_row$Dataset,
                           " was not found as an entry of RawNovorList."))

  #making sure Type is either Contributor or Primary
  assertthat::assert_that(dataTracker_row$Type %in% c("Primary", "Contributor"),
      msg = "Invalid value in Type column. Must be either Primary or Contributor.")

  #if Type=Contributor then returning subset
  #else if Type=Primary then returning all the Novor results for that Dataset

  if (dataTracker_row$Type == "Contributor")
  {
    #making sure scans isn't emtpy
    assertthat::assert_that(length(dataTracker_row$ContributorScans[[1]]) > 0,
                msg = "scans length must be greater than 0")

    #making sure scans is an atomic vector
    assertthat::assert_that(is.atomic(dataTracker_row$ContributorScans[[1]]),
                msg = "scans is not atomic, please make sure it is an atomic vector of integers, not a list")

    #making sure scans is numeric
    assertthat::assert_that(is.numeric(dataTracker_row$ContributorScans[[1]]),
                msg = "scans is not numeric, scans must be an atomic vector of integers")

    #extracting the dataset's Novor results
    rawresults <- RawNovorList[[
                    match(dataTracker_row$Dataset, names(RawNovorList))
                              ]]

    #subsetting the results according to the scan numbers in dataTracker_row
    export <- rawresults[
      rawresults$ScanNum %in% dataTracker_row$ContributorScans[[1]],
      ]


  }else #Type = Primary
  {
    #extracting dataset's Novor results
    export <- RawNovorList[[match(dataTracker_row$Dataset, names(RawNovorList))]]

  }
  #adding column for SimID
  export$SimID <- dataTracker_row$SimID

  #adding column for origin dataset
  export$OriginDataset <- dataTracker_row$Dataset
  return(export)

}


#' Simulated Novor results from multiple datasets
#'
#' \code{MixMultipleNovor} takes the subset of a dataTracker dataframe which
#' corresponds to a single simulation and calls \code{\link{subsetSingleNovor}}
#' using adply on each row and merges the results to create the simulated
#' dataframe of Novor results. It then creates and adds a simulated dataset name
#' as a new column.
#'
#' @param dataTracker_rows subset of rows from a dataTracker dataframe which
#'   correspond to the component datasets of a single simulation
#' @param RawNovorList list of dataframes with each dataframe containing the
#'   Novor results for a particular dataset. The name of each entry in the list
#'   be the dataset name. The dataframe must have a column containing scan
#'   numbers which must be called "ScanNum".
#'
#' @return a dataframe of simulated Novor results created from the datasets and
#'   scans listed in each row in dataTracker_row
#'
#' @family Multi-contributor Simulations
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#' @importFrom assertthat assert_that
#' @importFrom plyr adply
#'
#' @export
MixMultipleNovor <- function(dataTracker_rows, RawNovorList)
{
  assertthat::assert_that(length(unique(dataTracker_rows$SimID)) == 1,
                          msg = "dataTracker_rows must only contain 1 SimID")

  simNovor <- plyr::adply(dataTracker_rows,
                          1,
                          subsetSingleNovor,
                          .expand = FALSE,
                          .id = NULL,
                          RawNovorList = RawNovorList)

  simNovor$SimDataset <- paste0(dataTracker_rows$Project[1],
                                "_NC_",
                                nrow(dataTracker_rows),
                                '_SimID_',
                                dataTracker_rows$SimID[1])

  return(simNovor)
}



