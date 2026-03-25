

#' Convert raw files using ThermoRawFileParser
#'
#' This function is a wrapper for ThermoRawFileParser.exe which converts raw
#' files to mgf or mzML format. The mgf files it produces contain the `SCANS=`
#' line which allows Novor to include scan numbers in its output files.
#' ThermoRawFileParser.exe must be in the systems PATH such that it can run
#' regardless of its installation location.
#'
#' @param input_file Path to the input RAW file.
#' @param output_dir Path to the desired output directory. If NULL (default) the
#'   output will be in the same directory as input_file.
#' @param output_format Desired output format. Can be `mgf` or `mzML`.
#'
#' @return Returns the path to the output file.
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#' @seealso \code{\link{runNovor}} ; \code{\link{mzMLtoMGF}}
#'
#' @export
#' @importFrom assertthat assert_that
#' @importFrom stringr str_to_lower str_sub
convert_raw_file <- function(input_file,
                             output_dir = NULL,
                             output_format) {
  # checking inputs
  assertthat::assert_that(file.exists(input_file),
                          msg = paste0("input_file: ",
                                       input_file,
                                       " not found."))

  input_extension <- stringr::str_sub(input_file,-3)
  assertthat::assert_that(stringr::str_to_lower(input_extension) == "raw",
                          msg = "input_file must be a .raw file.")

  if (!is.null(output_dir)) {
    assertthat::assert_that(dir.exists(output_dir),
                            msg = paste0("output_dir: ",
                                         output_dir,
                                         " not found."))
  }


  assertthat::assert_that((output_format == "mgf" |
                             output_format == "mzML"),
                          msg = "output_format must be 'mgf' or 'mzML'")

  # assigning output_dir to be same as input_file
  # if output_dir is not specified
  if (is.null(output_dir)) {
    output_dir <- normalizePath(dirname(input_file))
  } else {
    output_dir <- normalizePath(output_dir)
  }




  # getting format value for command line call
  if (output_format == "mgf") {
    format_value <- 0
  } else if (output_format == "mzML") {
    format_value <- 2
  } else {
    stop("output_format must be `mgf`` or `mzML`")
  }

  input_file <- normalizePath(input_file)

  # making command line string
  cmdstring <- paste0(
    "ThermoRawFileParser.exe -i=",
    shQuote(input_file, type = "cmd"),
    " -f=",
    format_value,
    " -o=",
    shQuote(output_dir, type = "cmd")
  )


  # running on command line
  shell(cmd = cmdstring)


  # making output file name by replacing the extension of the input_file
  output_name <-
    gsub(input_extension, output_format, basename(input_file))


  # making output path
  output_path <-
    file.path(normalizePath(output_dir), output_name, fsep = "\\")

  # verifying that output exists
  assertthat::assert_that(file.exists(output_path),
                          msg = paste0("Output file: ",
                          output_path, " not found."))

 return(output_path)
}
