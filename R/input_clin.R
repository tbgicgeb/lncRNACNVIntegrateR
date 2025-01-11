#' Preparing input files for clinical data
#'
#' @param file_path Optional character string. If provided, a custom path to a tab-delimited
#'                  clinical file will be used. If not provided, the default system file
#'                  will be used from the package.
#'
#' @return A data frame of selected clinical data
#'
#' @import utils
#' @import readr
#' @export
#' @examples
#' clin_df <- input_clin()  # Uses the default system file
#' clin_df <- input_clin("/custom/path/to/clinical_data.txt.gz")

input_clin <- function(file_path = NULL) {
  if (is.null(file_path)) {
    # Use the system file if no path is provided
    file_path <- system.file("extdata", "clinical_data_COAD.txt.gz", package = "lncRNACNVIntegrateR")
  }

  if (!file.exists(file_path)) {
    stop("The file does not exist: ", file_path)
  }

  if (grepl("\\.gz$", file_path)) {
    con <- gzfile(file_path, "rt")
    clin_df <- read.table(con, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
    close(con)
  } else {
    clin_df <- read.table(file_path, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
  }

  return(clin_df)  # Return the clinical data frame
}
