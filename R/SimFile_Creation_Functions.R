#this script contains functions needed for Simulated dataset creation

#'Randomly select scans representing n proteins and m peptides
#'
#'\code{selectScans} returns a set of scan numbers associated with a randomly
#'selected number of proteins(Nprot) and randomly selected number of
#'peptides(Mpep) per protein. SimID is used as the seed for the random
#'selection.
#'
#'This funciton first calculates the number of peptides per protein and filters
#'out all proteins with < 2 peptides per protein. Of the remaining proteins,
#'"Nprot" are randomly selected. For each selected protein, a random set of
#'"Mpep" peptides are selected. All scans associated with those peptides are
#'then returned. If a protein has fewer peptides than Mpep, then all peptides
#'from that protein are used. All scans associated with each selected peptide
#'are returned. \strong{Warning:} This function only works with "unrolled" database
#'results, meaning one row for each scan, peptide, and protein combination.
#'
#'@section Random Number Generator (RNG) State:
#'  The random number generation settings are a
#'  \strong{global variable}. Within \code{selectScans} they are changed
#'  according to the \code{rng_verison} parameter then the original state is
#'  restored before the function exits. Version 3.6.0 added a new setting
#'  \code{sample.kind} with values 'Rounding'(< 3.6.0) and 'Rejection'(3.6.0
#'  behavior). See \code{\link[base]{Random}} for more details.
#'
#'@param databaseResults a dataframe containing unrolled MSGF+ results (one
#'  protein per row). Must have columns "cleanseq", "Protein", and "ScanNum"
#'@param Nprot how many proteins to select. Note:  It is possible to end up with
#'  more contributor proteins represented in the output scans if a selected
#'  peptide belongs to more than one protein.
#'@param Mpep how many peptides per protein to select
#'@param SimID value used as the seed of \code{\link[base]{sample}} when
#'  selecting peptides and proteins. Set as the ID of the simulated dataset when
#'  called by \code{\link{MakeSimFile}}
#'
#'@param rng_version R version string (Ex: '3.5.3') or 'R_version_default'
#'  (default). The RNG (random number generator) state within selectScans() is
#'  set to corrispond to the defaults of the R version specified in
#'  \code{rng_version}. If 'R_version_default' then the default settings for the
#'  currently installed R version are used.
#'
#'
#'@return Returns an integer vector containing all scan numbers associated with
#'  the selected peptides and proteins
#'
#'@author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}; Natalie Heller
#'
#'@family Simulation Creation
#'
#'@importFrom assertthat assert_that
#'@importFrom reshape2 dcast
selectScans <- function(databaseResults,
                        Nprot,
                        Mpep,
                        SimID,
                        rng_version = "R_version_default")
{

  #changing RNG state and restoring after function exit
  old_rng_state <- RNGkind()

  on.exit(RNGkind(kind = old_rng_state[1], normal.kind = old_rng_state[2]))

  if (getRversion() >= "3.6.0")
  {
    on.exit(RNGkind(sample.kind = old_rng_state[3]), add = TRUE)
  }

  if (rng_version == "R_version_default")
  {
    suppressWarnings(RNGversion(getRversion()))
  } else
  {
    suppressWarnings(RNGversion(rng_version))
  }


  #checking for correct column names
  assertthat::assert_that("cleanseq" %in% colnames(databaseResults),
              msg = "no 'cleanseq' column found in databaseResults")

  assertthat::assert_that("Protein" %in% colnames(databaseResults),
              msg = "no 'Protein' column found in databaseResults")

  assertthat::assert_that("ScanNum" %in% colnames(databaseResults),
              msg = "no 'ScanNum' column found in databaseResults")


  #creating dataframe holding peptide-protein relationships

  subset <- databaseResults[,c("cleanseq", "Protein")]
  subset <- subset[!duplicated(subset),]

  proteintable <-
    reshape2::dcast(subset,
                    Protein ~ .,
                    value.var = "Protein",
                    fun.aggregate = length)
  colnames(proteintable)[2] <- "NumPeptides"

  #keeping only proteins with >= 2 peptides
  proteintable <- proteintable[proteintable$NumPeptides >= 2,]

  #randomly selecting protein accessions
  set.seed(seed = SimID)
  protindex <-
    suppressWarnings(sample(
      x = nrow(proteintable),
      size = Nprot,
      replace = FALSE
    ))
  selected_prots <- proteintable$Protein[protindex]

  #algorithm to randomly select proteins and peptides
  out_scans <- c()
  out_peptides <- c()
  for (i in selected_prots)
  {
    peps <- subset$cleanseq[subset$Protein == i]
    peps <- unique(peps)

    if (length(peps) <= Mpep)
    {
      #if the number of peptides for the protein is <= Mpep
      #then return all scans for all peps

      #selecting corrisponding scans
      addscans <- unique(databaseResults$ScanNum[databaseResults$cleanseq %in% peps])
      #adding them to output
      out_scans <- c(out_scans,addscans)

    }else
    {
      #taking a random sample of peptides
      set.seed(seed = SimID)
      peps_indx <- suppressWarnings(sample(x = length(peps), size = Mpep, replace = FALSE))
      newpeps <- peps[peps_indx]
      #selecting corrisponding scans
      addscans <- unique(databaseResults$ScanNum[databaseResults$cleanseq
                                                 %in% newpeps])
      #adding them to output
      out_scans <- c(out_scans,addscans)
    }
  }
  out_scans <- unique(out_scans)

  return(out_scans)
}


#' Subset a contributor mzML and merge with a primary mzML
#'
#' \code{SubsetMergeMzmls} first makes a new mzML conatining a specified set of
#' scans from the contributor mzML using the MSConvert filter parameter. It then
#' merges the subset contributor mzML with the primary mzML to make a new
#' simulated mzML with a specified file name. The contributor and primary mzMLs
#' must be in the working directory. MSConvert must be installed and in the
#' computer's PATH.
#'
#' @param scans a integer vector containing scan numbers to select from the
#'   contributor
#' @param contributorName dataset name of the contributor mzML file (no ".mzML"
#'   extension)
#' @param primaryName dataset name of the primary mzML (no ".mzML" extension)
#' @param simname filename of the simulated dataset (\strong{with} the ".mzML"
#'   extension)
#' @param simID ID number of the simulated dataset. Used to name the mzML
#'   containing the subset of scans from the contributor mzML and passed to
#'   \code{\link{selectScans}}
#'
#' @return Returns the name of the simulated dataset that was created
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#' @family Simulation Creation
#'
#' @importFrom assertthat assert_that
#' @export
SubsetMergeMzmls <- function(scans,
                             contributorName,
                             primaryName,
                             simname,
                             simID)
{
  contributor <- paste0(contributorName, ".mzML")
  primary <- paste0(primaryName, ".mzML")

  #making sure input files exist
  assertthat::assert_that(file.exists(contributor),
                          msg = paste0(contributor, " was not found."))

  assertthat::assert_that(file.exists(primary),
                          msg = paste0(primary, " was not found."))

  #making sure there are scans
  assertthat::assert_that(length(scans) > 0,
                          msg = "There are no input contributor scans.")

  #making string for filter.txt
  scanstring <- toString(scans)
  filterline <- paste0("filter=scanNumber ", scanstring)

  #writing filter.txt
  fileConn <- file("filter.txt", open = "wt")
  writeLines(filterline, fileConn)
  close(fileConn)

  #making command line string
  subsetfile <-
    paste0(substr(contributor, 1, nchar(contributor) - 5),
           '_SimID',
           simID ,
           "_subset.mzML")

  cmdstring <- paste0("msconvert ",
                      contributor,
                      " -c filter.txt --outfile ",
                      subsetfile)
  system(cmdstring)


  #merging subset with primary mzml

  #making sure subset file exists
  assertthat::assert_that(
    file.exists(subsetfile),
    msg = paste0("The contributor subset file, ", subsetfile,
                 ", was not found.")
  )

  #making file list
  fileConn <- file("filelist.txt", open = "wt")
  writeLines(primary, fileConn)
  writeLines(subsetfile, fileConn)
  close(fileConn)

  #making command string
  mergestring <- paste0("msconvert -f filelist.txt --outfile ",
                        simname,
                        " --merge")
  system(mergestring)

  #deleting filelist.txt and filter.txt
  file.remove("filelist.txt")
  file.remove("filter.txt")

  #making sure simulated file exists
  assertthat::assert_that(file.exists(simname),
                          msg = paste0(simname, " was not found."))
  assertthat::assert_that(file.size(simname) > 0,
                          msg = paste0(simname, " size is 0."))

  return(simname)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#' Make a simulated dataset from a SimRecipes row and MSGF+ results
#'
#' \code{MakeSimFile} is meant to be used with ddply on a SimRecipes table and a
#' list of contributor MSGF+ results. It uses parameters stored the SimRecipes
#' table and a dataframe of MSGF+ results to call \code{\link{SubsetMergeMzmls}}
#' and \code{\link{selectScans}}. It keeps track of what scans, peptides, and
#' proteins were added from the contributor mzML along with when the simulated
#' dataset was made. This information is added to the input recipeRow dataframe
#' and returned. All input mzMLs must be in the working directory.
#'
#' @param recipeRow dataframe with a single row. At minimum must have columns:
#'   "ContributorDataset", "SimID", "PrimaryDataset", "PrimaryOrganism",
#'   "ContributorOrganism", "TargetNumContributorProteins", and
#'   "TargetContributorPep_per_Prot". Any extra columns will be retained and
#'   returned in output.
#' @param dbResults dataframe containing unrolled MSGF+ results (one protein per
#'   row). Must have columns "cleanseq", "Protein", and "ScanNum"
#'
#' @inheritParams selectScans
#'
#' @return Returns the input recipeRow dataframe with added columns containing
#'   what and how many contributor scans, peptides, and proteins were added.
#'   Added
#'   scans, peptides, and proteins are stored as a list of character/integer
#'   vectors.
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#' @family Simulation Creation
#'
#' @export
#' @importFrom assertthat assert_that
MakeSimFile <- function(recipeRow, dbResults, rng_version = "R_version_default")
{

  #making sure contributor dataset is in dbResults
  assertthat::assert_that(recipeRow$ContributorDataset %in% names(dbResults),
              msg = paste0("contributor Dataset: ",recipeRow$ContributorDataset,
                           " was not found in dbResults list names."))

  #finding database results for contributor dataset
  contributorDBresults <- dbResults[[match(recipeRow$ContributorDataset,
                                         names(dbResults))]]

  #stupid user resistance
  assertthat::assert_that(nrow(contributorDBresults) > 0,
              msg = "Number of rows in contributorDBresults must be > 0.")

  assertthat::assert_that(recipeRow$TargetNumContributorProteins > 0,
              msg = "TargetNumContributorProteins must be > 0.")

  assertthat::assert_that(recipeRow$TargetContributorPep_per_Prot > 0,
              msg = "contributorPep_per_Prot must be > 0.")

  assertthat::assert_that(is.numeric(recipeRow$SimID),
              msg = "SimID is not numeirc. Make sure SimID is an int.")

  assertthat::assert_that(recipeRow$SimID > 0, msg = "SimID must be > 0.")

  #calling selectScans
  ContributorScans <- selectScans(databaseResults = contributorDBresults,
                                Nprot = recipeRow$TargetNumContributorProteins,
                                Mpep = recipeRow$TargetContributorPep_per_Prot,
                                SimID = recipeRow$SimID,
                                rng_version = rng_version)

  assertthat::assert_that(length(ContributorScans) > 0,
                          msg = "selectScans returned no scans")

  #finding what peptides were selected
  ContributorPeptides <- unique(contributorDBresults$cleanseq[
                                  contributorDBresults$ScanNum %in%
                                  ContributorScans])

  #finding what proteins were selected
  ContributorProteins <- unique(contributorDBresults$Protein[
                          contributorDBresults$ScanNum %in% ContributorScans])

  #calling SubsetMergeMzmls
  madeSimFile <- SubsetMergeMzmls(scans = ContributorScans,
                                contributorName = recipeRow$ContributorDataset,
                                primaryName = recipeRow$PrimaryDataset,
                                simID = recipeRow$SimID,
                                simname = recipeRow$SimFileName)

  #creating dataTracker row for output

  #adding contributor scan information
  recipeRow$ContributorScans <- list(ContributorScans)
  recipeRow$NumContributorScans <- length(ContributorScans)

  #adding contributor peptide information
  recipeRow$ContributorPeptides <- list(ContributorPeptides)
  recipeRow$NumContributorPeptides <- length(ContributorPeptides)

  #adding contributor protein information
  recipeRow$ContributorProteins <- list(ContributorProteins)
  recipeRow$NumContributorProteins <- length(ContributorProteins)

  #adding date and time simulated file was made
  recipeRow$Birthday <- Sys.time()

  return(recipeRow)
}

#' Select scans from MSGF+ results with SimRecipe input
#'
#' \code{BatchSelectScans} uses parameters stored recipeRow and a dataframe of
#' MSGF+ results to call \code{\link{selectScans}}. It is meant to be used with
#' \code{\link[plyr]{ddply}} on a SimRecipes table and a list of contributor
#' MSGF+ results. The selected scans, and their associated peptides and
#' proteins, are added to the input recipeRow dataframe and returned.
#'
#' @param recipeRow dataframe with a single row. At minimum must have columns:
#'   "ContributorDataset", "SimID", "TargetNumContributorProteins", and
#'   "TargetContributorPep_per_Prot". Any extra columns will be retained and
#'   returned in output.
#' @param dbResults dataframe containing unrolled MSGF+ results (one protein per
#'   row). Must have columns "cleanseq", "Protein", and "ScanNum"
#' @inheritParams selectScans
#'
#' @return Returns the input recipeRow dataframe with added columns containing
#'   what and how many contributor scans, peptides, and proteins were added.
#'   Added scans, peptides, and proteins are stored as a list of
#'   character/integer vector.
#'
#' @importFrom assertthat assert_that
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#' @family Simulation Creation
#'
#' @export
#' @importFrom assertthat assert_that
BatchSelectScans <- function(recipeRow,
                             dbResults,
                             rng_version = "R_version_default")
{
  #making sure contributor dataset is in dbResults
  assertthat::assert_that(recipeRow$ContributorDataset %in% names(dbResults),
              msg = paste0("contributor Dataset: ",
              recipeRow$ContributorDataset,
              "was not found in dbResults list names."))

  #finding database results for contributor dataset
  contributorDBresults <-
    dbResults[[match(recipeRow$ContributorDataset, names(dbResults))]]

  #stupid user resistance
  assertthat::assert_that(nrow(contributorDBresults) > 0,
              msg = "Number of rows in contributorDBresults must be > 0.")
  assertthat::assert_that(recipeRow$TargetNumContributorProteins > 0,
              msg = "TargetNumContributorProteins must be > 0.")
  assertthat::assert_that(recipeRow$TargetContributorPep_per_Prot > 0,
              msg = "contributorPep_per_Prot must be > 0.")
  assertthat::assert_that(is.numeric(recipeRow$SimID),
              msg = "SimID is not numeirc. Make sure SimID is an int.")
  assertthat::assert_that(recipeRow$SimID > 0, msg = "SimID must be > 0.")

  #calling selectScans
  ContributorScans <- selectScans(databaseResults = contributorDBresults,
                                Nprot = recipeRow$TargetNumContributorProteins,
                                Mpep = recipeRow$TargetContributorPep_per_Prot,
                                SimID = recipeRow$SimID,
                                rng_version = rng_version
                                )

  assertthat::assert_that(length(ContributorScans) > 0,
                          msg = "selectScans returned no scans")

  #finding what peptides were selected
  ContributorPeptides <-
    unique(contributorDBresults$cleanseq[
      contributorDBresults$ScanNum %in% ContributorScans])

  #finding what proteins were selected
  ContributorProteins <-
    unique(contributorDBresults$Protein[
      contributorDBresults$ScanNum %in% ContributorScans])

  #adding contributor scan information
  recipeRow$ContributorScans <- list(ContributorScans)
  recipeRow$NumContributorScans <- length(ContributorScans)

  #adding contributor peptide information
  recipeRow$ContributorPeptides <- list(ContributorPeptides)
  recipeRow$NumContributorPeptides <- length(ContributorPeptides)

  #adding contributor protein information
  recipeRow$ContributorProteins <- list(ContributorProteins)
  recipeRow$NumContributorProteins <- length(ContributorProteins)

  return(recipeRow)
}

#' Merge contributor and primary Novor results
#'
#' \code{MixNovor} extracts the Novor results corrisponding to the scan numbers
#' in ContributorScans and adds them to the Novor results of the primary
#' dataset.
#' \strong{Note:}To avoid duplication of scan numbers, it is reccomended to
#' have a column in the input Novor result dataframes which lists the origin
#' file. This is added if you use \code{\link{importSingleNovor}} to import the
#' Novor results from the csv files.
#'
#' @param PrimaryDataset name of the primary dataset
#' @param ContributorDataset name of the contributor dataset
#' @param NovorResultList list of dataframes containing Novor results. The
#'   column containing scan numbers must be called "ScanNum". PrimaryDataset and
#'   ContributorDataset must match to entry names in NovorResultList.
#' @param ContributorScans vector of scan numbers to take from the contributor
#'   results and add to the primary results
#' @param SimID Simulated dataset ID, will be added as a column to the output
#'
#' @return Returns a dataframe of Novor results consisting of all the primary
#'   results and the contributor results associated with the scan numbers in
#'   ContributorScans. A column containing the input SimID is also added.
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#' @family Simulation Creation
#' @family Novor Simulation
#'
#' @export
#' @importFrom assertthat assert_that
MixNovor <- function(PrimaryDataset,
                     ContributorDataset,
                     NovorResultList,
                     ContributorScans,
                     SimID)
{
  #making sure contributor and primary datasets are in NovorResultList
  assertthat::assert_that(ContributorDataset %in% names(NovorResultList),
              msg = paste0("Contributor Dataset: ",ContributorDataset,
                           " was not found in NovorResultList list names."))

  assertthat::assert_that(PrimaryDataset %in% names(NovorResultList),
              msg = paste0("Primary Dataset: ",PrimaryDataset,
                           " was not found in NovorResultList list names."))

  #making sure ContributorScans isn't emtpy
  assertthat::assert_that(length(ContributorScans) > 0,
              msg = "ContributorScans length must be greater than 0")

  #making sure ContributorScans is a vector of ints
  assertthat::assert_that(is.atomic(ContributorScans),
      msg = "ContributorScans is not atomic, please make sure it is an atomic vector of integers, not a list")

  assertthat::assert_that(is.numeric(ContributorScans),
      msg = "ContributorScans is not numeric, ContributorScans must be an atomic vector of integers")



  #extracting database results for contributor and primary datasets
  contributorNovor <- NovorResultList[[match(ContributorDataset,
                                           names(NovorResultList))]]

  primaryNovor <- NovorResultList[[match(PrimaryDataset,names(NovorResultList))]]

  #subsetting contributor using scan numbers
  contributorSubset <-
    contributorNovor[contributorNovor$ScanNum %in% ContributorScans, ]

  #adding to primary Novor results
  mixedNovor <- rbind(primaryNovor, contributorSubset)

  #adding column for SimID
  mixedNovor$SimID <- SimID

  return(mixedNovor)

}

#' Wrapper of MixNovor() using a dataTracker row as input
#'
#' \code{BatchMixNovor} is a wrapper of \code{\link{MixNovor}} which takes
#' parameters for \code{\link{MixNovor}} from dataTracker_row. This is necessary
#' so you can call \code{\link[plyr]{ddply}} on a dataTracker dataframe to make
#' batches of simulated Novor results. After calling \code{\link{MixNovor}} it
#' adds a column defining the simulated dataset name which will be needed to
#' make peptide list files. It gets this value from the SimDataset column in
#' dataTracker_row. You can use the \code{export_csv} and \code{outputFileDir}
#' parameters to export the simulated raw Novor results to a csv file in
#' addition to returning them as a dataframe. If you do not want the results
#' returned as a dataframe, call \code{BatchMixNovor} using
#' \code{\link[plyr]{d_ply}}.
#'
#' @param dataTracker_row a dataframe containing a single row. At minimum must
#'   contain columns "ContributorDataset", "PrimaryDataset", "ContributorScans",
#'   "SimID", and "SimDataset".
#'
#' @param export_csv Should the simulated raw Novor results be exported to a csv
#'   file? Default FALSE
#' @param outputFileDir path to the desired output folder for raw Novor results.
#'   Defaults to the current working directory. The filename is taken from the
#'   SimDataset column with "_Raw_Novor.csv" appended.
#'
#' @inheritParams MixNovor
#'
#' @return Returns a dataframe of Novor results consisting of all the primary
#'   results and the contributor results associated with the scan numbers in
#'   ContributorScans. A column containing the input SimID and SimDataset is also
#'   added.
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#' @family Simulation Creation
#' @family Novor Simulation
#'
#' @export
#' @importFrom assertthat assert_that
BatchMixNovor <- function(dataTracker_row,
                          NovorResultList,
                          export_csv = FALSE,
                          outputFileDir = ".")
{
  assertthat::assert_that(!is.null(dataTracker_row$SimDataset),
                          msg = "SimDataset column not found in dataTracker_row")

  assertthat::assert_that(!is.null(dataTracker_row$SimID),
                          msg = "SimID column not found in dataTracker_row")

  mixed <- MixNovor(PrimaryDataset = dataTracker_row$PrimaryDataset,
                    ContributorDataset = dataTracker_row$ContributorDataset,
                    NovorResultList = NovorResultList,
                    ContributorScans = dataTracker_row$ContributorScans[[1]],
                    SimID = dataTracker_row$SimID)

  #adding column for Simulated filename
  mixed$SimDataset <- dataTracker_row$SimDataset

  if (export_csv)
  {
    assertthat::assert_that(dir.exists(outputFileDir),
                            msg = paste("outputFileDir",
                                        outputFileDir,
                                        "not found."))

    export_filename <- paste0(dataTracker_row$SimDataset,"_Raw_Novor.csv")

    outputFilePath <- paste(normalizePath(outputFileDir),
                            export_filename, sep = "\\")

    write.csv(mixed, file = outputFilePath, row.names = FALSE)

    assertthat::assert_that(file.exists(outputFilePath),
                            msg = paste0("Raw Novor Result file: ",
                                         outputFilePath, " not found."))
  }


  return(mixed)
}

#'MixNovor, FiltNovor, makePepFile using dataTracker row
#'
#'\code{MixFilterDigest} is meant to be used with large sets of simulated
#'datasets (> 50) to avoid holding more than 50 sets of simulated Novor results
#'in memory. To do this it calls \code{\link{MixNovor}},
#'\code{\link{filterNovor}}, and \code{\link{makePepFile}} givin a single
#'dataTracker row. This bypasses \code{\link{BatchMixNovor}}. The downside of
#'using this function is that you can't interrogate the the raw and filtered
#'mixed novor results.
#'
#'@param dataTracker_row a dataframe containing a single row of the format
#'  returned by \code{\link{BatchSelectScans}}. At minimum must contain columns
#'  "ContributorDataset", "PrimaryDataset", "ContributorScans", "SimID", and
#'  "SimDataset".
#'@param NovorResultList list of dataframes containing Novor results. The column
#'  containing scan numbers must be called "ScanNum". PrimaryDataset and
#'  ContributorDataset must match to entry names in NovorResultList.
#'@inheritParams filterNovor
#'@inheritParams makePepFile
#'@return Returns dataTracker_row with new columns containing the name of the
#'  peptide list file (PeptideListFile) and the number of digested peptides
#'  (Num_Digested_Peptides)
#'
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#' @family Simulation Creation
#' @family Novor Simulation
#'
#'@export
#'@importFrom assertthat assert_that
MixFilterDigest <- function(dataTracker_row,
                            NovorResultList,
                            outputFileDir = ".",
                            ppmErrorMin = -15,
                            ppmErrorMax = 15,
                            minNovorScore = 70,
                            minLength = 6,
                            datasetCol = "SimDataset",
                            sequenceCol = "cleanseq",
                            prolineBlock = TRUE,
                            discardLessThan = 6)
{
  assertthat::assert_that(!is.null(dataTracker_row$SimDataset),
              msg = "SimDataset column not found in dataTracker_row")

  assertthat::assert_that(!is.null(dataTracker_row$SimID),
              msg = "SimID column not found in dataTracker_row")

  #making mixed Novor results
  RawMixed <- MixNovor(PrimaryDataset = dataTracker_row$PrimaryDataset,
                       ContributorDataset = dataTracker_row$ContributorDataset,
                       NovorResultList = NovorResultList,
                       ContributorScans = dataTracker_row$ContributorScans[[1]],
                       SimID = dataTracker_row$SimID)

  #adding column for Simulated filename
  RawMixed$SimDataset <- dataTracker_row$SimDataset

  #filtering Novor results
  filterMixed <- filterNovor(RawMixed, ppmErrorMin = ppmErrorMin,
                             ppmErrorMax = ppmErrorMax,
                             minNovorScore = minNovorScore,
                             minLength = minLength)

  #writing peptide list file
  fileNcount <- makePepFile(filterMixed, outputFileDir = outputFileDir,
                            datasetCol = datasetCol,
                            sequenceCol = sequenceCol,
                            prolineBlock = prolineBlock,
                            discardLessThan = discardLessThan)

  #adding peptide list file name and number of digested peptides
  dataTracker_row$PeptideListFile <- fileNcount$PeptideListFile
  dataTracker_row$Num_Digested_Peptides <- fileNcount$Num_Digested_Peptides

  return(dataTracker_row)
}








