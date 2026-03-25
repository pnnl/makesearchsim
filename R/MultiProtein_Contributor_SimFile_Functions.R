#' MakeMultiSimFile
#' Creates a new mzML file that combines contains scans for
#' specific proteins from a database search, with a primary mzML file. This
#' function can take multiple proteins from multiple database searches, as long
#' as each unique combination of contributor file and protein has a different
#' row.
#'
#' @param multi_DT A dataframe which contains at least the following columns.
#'   Dataset: A unique name for the final mzML file, without the extension;
#'   SimID: Unique identifier; PrimaryDataset: Filename of the primary dataset
#'   without the extension.; ContributorDataset: Filename of the contributor
#'   dataset without the extension; ContributorProteins: Protein identifiers
#'   that exactly match proteins in the MS-GF+ output.  Each row corresponds to
#'   a unique combination of primary file, contributor file, and contributor
#'   protein.
#' @param MSGF_results A dataframe of database search PSMs. The follow columns
#'   are required: Dataset, Scan, and Protein. The Dataset column is expected to
#'   have the exact same dataset name of the one_contrib_DT$ContributorDataset.
#'   Additionally, The Protein column is expected to have the exact string that
#'   is found in the one_contrib_DT$ContributorProteins column.
#' @param working_dir The full path to the mzML files and the location that the
#'   new mzML files will be written to.
#'
#' @return Returns the input dataframe 'multi_DT' with two new columns. One is
#'   "SimFileName", which contains the name of the new mzML file created with
#'   the PrimaryDataset and the scans from the proteins from the
#'   ContributorDataset. The other is "Scans", which is a list of the scans
#'   selected for each protein from the MSGF_results dataframe.
#'
#' @family Multi-contributor Simulations
#' @family Specific Protein Simulations
#'
#' @author Isabelle O'Bryon \email{isabelle.obryon@@pnnl.gov}
#'
#'
#' @export
#' @importFrom plyr ddply
#' @importFrom plyr as.quoted
#' @importFrom plyr rbind.fill
#' @importFrom plyr adply
#' @importFrom assertthat assert_that
#'
#' @examples
#' \dontrun{
#' # In this example, the files "Sample_1.mzML", "Sample_2.mzML", "Sample_3.mzML",
#' # and "Primary_sample.mzML" are in the directory "C:/Users/XXX123/my_mzML_dir".
#' # The output file will be called "SimMulti_1.mzML" and it will also be in
#' # "C:/Users/XXX123/my_mzML_dir". The output file will have the scans that
#' # are in "filtered_MSGF_dataframe" that were identified in "Sample_1" as
#' # being from "FEL1B_FELCA", in "Sample_2" as being from "TRPV6_HUMAN", and
#' # in "Sample_3" as being from "Sample_3" along will all of the scans from
#' # "Primary_sample.mzML".
#' contrib_datasets <- c("Sample_1","Sample_2","Sample_3")
#' contrib_proteins <- c("FEL1B_FELCA", "TRPV6_HUMAN", "CTR2_HUMAN")
#' primary_dataset <- c("Primary_sample", "Primary_sample","Primary_sample")
#' example_multi_taxa_DT <- data.frame("Dataset"= c("SimMulti_1", "SimMulti_1", "SimMulti_1"),
#'                                     "SimID" = c(1, 1, 1),
#'                                     "PrimaryDataset" = primary_dataset,
#'                                     "ContributorDataset" = contrib_datasets,
#'                                     "ContributorProteins" = contrib_proteins,
#'                                     stringsAsFactors = FALSE)
#' mzML_dir <- "C:/Users/XXX123/my_mzML_dir"
#' load("filtered_MSGF_dataframe.RData")
#' example_multi_taxa_DT_results <- MakeMultiSimFile(example_multi_taxa_DT,
#' # filtered_MSGF_dataframe , mzML_dir)
#' }
MakeMultiSimFile <- function(multi_DT, MSGF_results, working_dir){
  # Change the working directory
  starting_wd <- getwd()
  setwd(working_dir)
  # Set the directory back
  on.exit(setwd(starting_wd), add = T)

  searched_files <- unique(MSGF_results$Dataset)

  # Check for an appropriate SimID
  assertthat::assert_that("SimID" %in% colnames(multi_DT),
                          msg = "no 'SimID' column found in multi_DT")
  assertthat::assert_that(is.numeric(multi_DT$SimID),
                          msg = "SimID is not numeirc. Make sure SimID is an int.")
  assertthat::assert_that(length(unique(multi_DT$SimID)) == 1,
                          msg = "More than one unique SimID. Only call MakeMultiSimFile with information for one new simulated mzML file.")
  assertthat::assert_that(unique(multi_DT$SimID) > 0, msg = "SimID must be > 0.")

  # Check for expected MS-GF+ columns
  assertthat::assert_that("Dataset" %in% colnames(MSGF_results),
                          msg = "no 'Dataset' column found in MSGF_results")
  assertthat::assert_that("Protein" %in% colnames(MSGF_results),
                          msg = "no 'Protein' column found in MSGF_results")
  assertthat::assert_that("Scan" %in% colnames(MSGF_results),
                          msg = "no 'Scan' column found in MSGF_results")


  # Call selectProteinScans
  multi_DT_scans <-
    plyr::adply(multi_DT, 1, selectProteinScans, dbResults = MSGF_results)
  # Call mergeMultMzML
  multi_DT_scans <- mergeMultMzML(multi_DT_scans)

  # Set the directory back
  setwd(starting_wd)
  return(multi_DT_scans)
}


#' selectProteinScans
#' This function is called by MakeMultiSimFile to select all
#' of the scans from one dataset for specific proteins.
#'
#' @param one_contrib_DT A dataframe that must have the columns:
#'   ContributorDataset, ContributorProteins. It is expected that there will
#'   only be one unique ContributorDataset. There may be one or more
#'   ContributorProteins. DT stands for data tracker.
#' @param dbResults A dataframe of database search PSMs. The follow columns are
#'   required: Dataset, Scan, and Protein. The Dataset column is expected to
#'   have the exact same dataset name of the one_contrib_DT$ContributorDataset.
#'   Additionally, The Protein column is expected to have the exact string that
#'   is found in the one_contrib_DT$ContributorProteins column.
#'
#' @return A dataframe that is the same as the input dataframe one_contrib_DT
#'   with the added column "Scans". The column "Scans" is a list of the scans
#'   from the database result dataframe that had the given protein identifier in
#'   the Protein column.
#' @export
#' @importFrom plyr ddply
#' @importFrom assertthat assert_that
#'
#' @family Specific Protein Simulations
#' @author Isabelle O'Bryon \email{isabelle.obryon@@pnnl.gov}
#'
#' @examples
#' \dontrun{
#' # This function is not meant to be called by the user, it is a function
#' # needed for the MakeMultiSimFile function.
#' # load("filtered_MSGF_dataframe.RData")
#' merge_sample_1_w_2 <- data.frame("Dataset"=  "SimMulti_1",
#'                                     "SimID" = 1,
#'                                     "PrimaryDataset" = "Sample_1",
#'                                     "ContributorDataset" = "Sample_2",
#'                                     "ContributorProteins" = "FEL1B_FELCA",
#'                                     stringsAsFactors = FALSE)
#' merge_sample_1_w_2_Scans <- selectProteinScans(merge_sample_1_w_2, filtered_MSGF_dataframe)
#' }
selectProteinScans <- function(one_contrib_DT, dbResults){
  one_contrib_DT$ContributorDataset <-
    as.character(one_contrib_DT$ContributorDataset)
  one_contrib_DT$ContributorProteins <-
    as.character(one_contrib_DT$ContributorProteins)
  # Check that the contributor dataset is in dbResults
  assertthat::assert_that(
    unique(one_contrib_DT$ContributorDataset) %in% dbResults$Dataset,
    msg = paste0(
      "contributor Dataset: ",
      unique(one_contrib_DT$ContributorDataset),
      "was not found in dbResults list names."
    )
  )

  # Subset database results for contributor dataset
  contributor_MSGF_results <-
    dbResults[dbResults$Dataset %in% unique(one_contrib_DT$ContributorDataset),]
  one_contrib_DT$Scans <- NA

  # For each of the proteins from this contributor the scans from MS-GF+ are selected if the protein identified
  # is present in the Protein column of the contributor_MSGF_results dataframe

  cur_contributor_MSGF_scans <-
    contributor_MSGF_results$Scan[grepl(one_contrib_DT$ContributorProteins,
                                        contributor_MSGF_results$Protein,
                                        fixed = T) == T]
  assertthat::assert_that(
    min(length(cur_contributor_MSGF_scans)) > 0,
    msg = paste0(
      "contributor Dataset ",
      unique(one_contrib_DT$ContributorDataset),
      " protein ",
      one_contrib_DT$ContributorProteins,
      " was not found in database results."
    )
  )
  one_contrib_DT$Scans <- list(cur_contributor_MSGF_scans)

  return(one_contrib_DT)
}


#' mergeMultMzML
#' This function takes in a dataframe and will merge together the
#' primary mzML with the contributors mzML files.
#'
#' @param multi_contrib_DT A dataframe that contains the contributors, proteins
#'   from the contributors, and the scans.
#'
#' @return Returns the input dataframe 'multi_contrib_DT' with a new column
#'   "SimFileName", which contains the name of the new mzML file created with
#'   the PrimaryDataset and the scans from the proteins from the
#'   ContributorDataset.
#'
#' @export
#' @importFrom plyr adply
#' @importFrom plyr ddply
#' @importFrom plyr dlply
#' @importFrom assertthat assert_that
#'
#' @family Multi-contributor Simulations
#' @family Specific Protein Simulations
#' @author Isabelle O'Bryon \email{isabelle.obryon@@pnnl.gov}
#'
#' @examples
#' \dontrun{
#' # This function is not meant to be called by the user, it is a function needed
#' # for the MakeMultiSimFile function.
#' multi_DT_scans <- ddply(multi_DT, .(ContributorDataset), selectProteinScans,  MSGF_results )
#  mergeMultMzML( multi_DT_scans )
#' }
mergeMultMzML <- function(multi_contrib_DT){
  # making sure input files exist
  primary_mzML <-
    paste0(unique(multi_contrib_DT$PrimaryDataset), ".mzML")
  assertthat::assert_that(file.exists(primary_mzML),
                          msg = paste0(primary_mzML, " was not found."))

  # Creates a new mzML file for each 'ContributorDataset' that only has the scans specified by
  # the 'Scans' column
  subset_mzML_files <-
    plyr::dlply(multi_contrib_DT,
                "ContributorDataset",
                subsetMzMLProteinScans)
  # Creates the file list, needed for the msconvert command
  fileConn <- file("filelist.txt", open = "wt")
  writeLines(primary_mzML, fileConn)
  writeLines(unlist(subset_mzML_files), fileConn)
  close(fileConn)

  # Makes the command string to have msconvert merge together all of the mzML files
  simname <- paste0(unique(multi_contrib_DT$Dataset), ".mzML")
  mergestring <-
    paste0("msconvert -f filelist.txt --outfile ", simname, " --merge")
  system(mergestring)

  multi_contrib_DT$SimFileName <- simname
  return(multi_contrib_DT)
}


#' subsetMzMLProteinScans
#' This function takes in a dataframe that has one unqiue
#' 'ContributorDataset'. The scans listed in the 'Scans' column are selected
#' from the Contributor mzML file and written to a separate mzML file.
#'
#' @param one_contrib_DT A dataframe that contains the contributors, proteins
#'   from the contributors, and the scans
#'
#' @return The name of the new mzML file
#'
#' @family Multi-contributor Simulations
#' @family Specific Protein Simulations
#' @author Isabelle O'Bryon \email{isabelle.obryon@@pnnl.gov}
#'
#' @export
#' @importFrom plyr ddply
#' @importFrom assertthat assert_that
#'
#' @examples
#' \dontrun{
#' # This function is not meant to be called by the user, it is a function needed
#' # for the MakeMultiSimFile function.
#' subset_mzML_files <- ddply(multi_contrib_DT, .(ContributorDataset), subsetMzMLProteinScans )
#' }
subsetMzMLProteinScans <- function(one_contrib_DT){
  contributor_mzML <-
    paste0(unique(one_contrib_DT$ContributorDataset), ".mzML")
  assertthat::assert_that(file.exists(contributor_mzML),
                          msg = paste0(contributor_mzML, " was not found."))
  #making sure there are scans
  assertthat::assert_that(length(unlist(one_contrib_DT$Scans)) > 0,
                          msg = "There are no input contributor scans.")
  scans <- unlist(one_contrib_DT$Scans)
  # making string for filter.txt
  scanstring <- toString(scans)
  filterline <- paste0("filter=scanNumber ", scanstring)
  filter_file_name <-
    paste0(
      "filter_",
      unique(one_contrib_DT$ContributorDataset),
      "_SimID_",
      unique(one_contrib_DT$SimID) ,
      ".txt"
    )

  # writing filter.txt
  fileConn <- file(filter_file_name, open = "wt")
  writeLines(filterline, fileConn)
  close(fileConn)

  # Makes the command string to have msconvert create a new file of the given scans
  subsetfile <-
    paste0(
      unique(one_contrib_DT$ContributorDataset),
      '_SimID',
      unique(one_contrib_DT$SimID),
      "_subset.mzML"
    )

  cmdstring <-
    paste0("msconvert ",
           contributor_mzML,
           " -c ",
           filter_file_name,
           " --outfile ",
           subsetfile)
  system(cmdstring, show.output.on.console = T)
  # Checks that the subset file exists
  assertthat::assert_that(
    file.exists(subsetfile),
    msg = paste0("The contributor subset file, ", subsetfile,
                 ", was not found.")
  )
  # Deletes the tmp filter.txt
  file.remove(filter_file_name)

  return(subsetfile)
}
