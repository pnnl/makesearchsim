context("Novor Functions")

#getting local Novor path
path_to_novor <-
  choose.files(caption = "Select path to novor.bat to use for test_that tests.")

#removing old novor result file if it exists
if (file.exists(test_path("test_Novor.csv"))) {
  file.remove(test_path("test_Novor.csv"))
}

test_that("runNovor input assertions", {
  expect_error(runNovor(MGF_File = "bla.mgf"), "bla.mgf not found.")

  expect_error(
    runNovor(MGF_File = test_path("test.mgf"), NovorParams = "param_not_here.txt"),
    "param_not_here.txt not found."
  )

  expect_error(
    runNovor(
      MGF_File = test_path("test.mgf"),
      NovorParams = test_path("params.txt")
    ),
    "argument \"Novor_Path\" is missing, with no default"
  )
})


test_that("Novor Runs with package params", {
  param_path <-
    system.file(
      "extdata",
      "NovorParams_fTol002_pTol15_MetOx_Alk_Novor_params.txt",
      package = "MakeSearchSim",
      mustWork = TRUE
    )

  novor_result_path <-
    runNovor(
      MGF_File = test_path("test.mgf"),
      NovorParams = param_path,
      Novor_Path = path_to_novor
    )

  expect_true(file.exists(novor_result_path))

})

test_that("importSingleNovor results within tolerance", {

  load(test_path("test_novor_df_truth.RData"))

  test_novor_new <- importSingleNovor(
    NovorOutputFile = test_path("test_Novor.csv"),
    MGF_file = test_path("test.mgf"))


  expect_identical(test_novor_df_truth$ScanNum, test_novor_new$ScanNum)

  expect_identical(test_novor_df_truth$NovorID, test_novor_new$NovorID)

  expect_identical(test_novor_df_truth$mzPrecursor, test_novor_new$mzPrecursor)

  expect_identical(test_novor_df_truth$pepMass_denovo, test_novor_new$pepMass_denovo)

  expect_identical(test_novor_df_truth$RT, test_novor_new$RT)

  # Novor running on different computers causes a small amount of peptide
  # assignments to be different, they have the same overall mass but have amino
  # acid flips and different PTM locations

  #calculating the intersection/union
  peptide_similarity <-
    length(intersect(test_novor_df_truth$peptide, test_novor_new$peptide)) /
    length(union(test_novor_df_truth$peptide, test_novor_new$peptide))

  print(paste0("peptide similarity: ", peptide_similarity))

 #will fail if more than 3% of peptides are different
 expect_equal(peptide_similarity, 1, tolerance = 0.03, scale = 1)

})


#removing Novor output file if it exists
#removing old novor result file if it exists
if (file.exists(test_path("test_Novor.csv"))) {
  file.remove(test_path("test_Novor.csv"))
}

