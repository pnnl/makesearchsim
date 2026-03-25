#' Example dataTracker dataframe for simulation via mzML
#'
#' An example dataTracker dataframe which results from the mzML version of the
#' simulated dataset creation pipeline run using SimRecipes_SimID5_8.csv
#'
#' @format A dataframe with 4 rows and 21 variables: \describe{
#'   \item{ContributorDataset}{Filename w/o extension of the contributor mzML, chr}
#'   \item{SimID}{ID number of the simulated dataset, int}
#'   \item{PrimaryDataset}{Filename w/o extension of the primary mzML, chr}
#'   \item{PrimaryOrganism}{Name of the primary organism, chr}
#'   \item{ContributorOrganism}{Name of the contributor organism, chr}
#'   \item{Distance}{taxonomic distance between contributor and primary, num}
#'   \item{TargetNumContributorProteins}{desired number of proteins selected from
#'   contributor, int}
#'   \item{TargetContributorPep_per_Prot}{desired number of
#'   peptides/protein from contributor, int}
#'   \item{SimFileName}{Filename of the simulated dataset w/ extension, chr}
#'   \item{ContributorScans}{contributor scan numbers, vector of ints}
#'   \item{NumContributorScans}{total number of contributor scans, int}
#'   \item{ContributorPeptides}{contributor peptides, vector of chr}
#'   \item{NumContributorPeptides}{total number of contributor peptides, int}
#'   \item{ContributorProteins}{contributor proteins, vector or chr}
#'   \item{NumContributorProteins}{total number of contributor proteins, int}
#'   \item{Birthday}{date and time simulated dataset was created, POSIXct format
#'   yyyy-mm-dd hh:mm:ss}
#'   \item{MGF_name}{Filename of the MGF file used by
#'   Novor, chr}
#'   \item{Raw_Novor_Filename}{Filename of the Novor output file,
#'   chr}
#'   \item{Novor_ParamFile}{Filename of the Novor parameter file, chr}
#'   \item{PeptideListFile}{File containing a list of digested peptides created
#'   by makePepFile from the filtered Novor results, chr}
#'   \item{Num_Digested_Peptides}{number of digested peptides contained in the
#'   PeptideListFile, int} }
#' @source Result from running the mzML simulation workflow using the SimRecipes_SimID5_8 SimRecipes table described in the 'Simulation workflow using simulated mzMLs' vignette.
"dataTracker_SimID5_8_mzMLex"


#' Example dual organism SimRecipes table for SimID 5-8
#'
#' This is an example of a SimRecipes table used as input to the workflows which
#' only make simulations from a single primary and contributor dataset. It is
#' used in the 'Simulation workflow using simulated mzMLs' and 'Make simulated
#' peptide lists by mixing prime and contributor Novor results' vignettes.
#'
#' @format A data frame with 4 rows and 8 variables: \describe{
#'   \item{ContributorDataset}{Filename (w/o extension) of the contributor mzML, chr}
#'   \item{SimID}{ID number of the simulated dataset, int}
#'   \item{PrimaryDataset}{Filename (w/o extension) of the primary mzML, chr}
#'   \item{PrimaryOrganism}{Name of the primary organism, chr}
#'   \item{ContributorOrganism}{Name of the contributor organism, chr}
#'   \item{Distance}{taxonomic distance between contributor and primary, num}
#'   \item{TargetNumContributorProteins}{desired number of proteins selected from
#'   contributor, int}
#'   \item{TargetContributorPep_per_Prot}{desired number of
#'   peptides/protein from contributor, int}
#'   \item{SimFileName}{Name chosen for the simulated dataset described by each row}
#'   }
#' @source Example included in package
"SimRecipes_SimID5_8"


#' Example long format (multi-organism) SimRecipes table for SimID 5-8
#'
#' This is an example of a long format SimRecipes table used as input to the
#' multi-organism simulation pipeline described in the "Make simulations from
#' many organisms using Novor results" vignette. It describes the same
#' simulations as SimRecipes_SimID5_8.
#'
#' @format A data frame with 4 rows and 8 variables: \describe{
#'   \item{SimID}{Numeric identifier unique to each simulation, int}
#'   \item{Dataset}{Name of a dataset, has no extension or "_Novor" suffix , chr}
#'   \item{Type}{Either "Primary" or "Contributor", chr}
#'   \item{Organism}{name of the organism associated with the dataset, chr}
#'   \item{TargetNumContributorProteins}{How many proteins to add from the
#'   contributor dataset  (This should be "NA" for Primary datasets), int}
#'   \item{TargetContributorPep_per_Prot}{How many peptides per protein to add
#'   from each selected contributor protein (This should be "NA" for Primary
#'   datasets), int}
#'   \item{Project}{Name of the project associated with the simulation which is
#'   used to name the output peptide list files, chr}
#'   }
#' @source Example input included in package
"SimRecipes5_8_long"


#' Example dataTracker dataframe from multi-organism workflow
#'
#' This is an example of the dataTracker dataframe which would result from
#' running the multi-organism workflow described in the "Make simulations from
#' many organisms using Novor results" vignette using SimRecipes5_8_long as the
#' input SimRecipes table.
#'
#' @format A data frame with 4 rows and 8 variables: \describe{
#'   \item{SimID}{Numeric identifier unique to each simulation, int}
#'   \item{Dataset}{Name of a dataset, has no extension or "_Novor" suffix , chr}
#'   \item{Type}{Either "Primary" or "Contributor", chr}
#'   \item{Organism}{name of the organism associated with the dataset, chr}
#'   \item{TargetNumContributorProteins}{How many proteins to add from the
#'   contributor dataset  (This should be "NA" for Primary datasets), int}
#'   \item{TargetContributorPep_per_Prot}{How many peptides per protein to add
#'   from each selected contributor protein (This should be "NA" for Primary
#'   datasets), int}
#'   \item{Project}{Name of the project associated with the simulation which is
#'   used to name the output peptide list files, chr}
#'   \item{ContributorScans}{list of scan numbers selected from the dataset (NA if a primary dataset), list of ints}
#'   \item{NumContributorScans}{Number of scans selected from the dataset (NA is a primary dataset), int}
#'   \item{ContributorPeptides}{list of peptides associated with the selected scans (NA if a primary dataset), list of character strings}
#'   \item{NumContributorPeptides}{number of unique peptides assocated with the selected scans (NA if a primary dataset), int}
#'   \item{ContributorProteins}{list of proteins associated with the selected scans (NA if a primary dataset), list of character strings}
#'   \item{NumContributorProteins}{number of unique proteins assocated with the selected scans (NA if a primary dataset), int}
#'   \item{PeptideListFile}{name of the file containing the list of digested peptides for each simulation, chr}
#'   \item{Num_Digested_Peptides}{number of digested peptides for each simulation , duplicated for each row associated with a simulation, int}
#'   }
#' @source Multi-dataset simulation workflow run on SimRecipes5_8_long
"multiDataTracker_SimID5_8"

