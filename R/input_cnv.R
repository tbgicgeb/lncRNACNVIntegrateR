#' Preparing input files for CNV data
#'
#' @param file_path Optional character string. If provided, a custom path to a tab-delimited
#'                  CNV file will be used. If not provided, the default system file
#'                  will be used from the package.
#'
#' @return A data frame of absolute CNV values
#'
#' @import utils
#' @export
#' @examples
#' cnv_df <- input_cnv()  # Uses the default system file
#' cnv_df <- input_cnv("/custom/path/to/CNV_data.txt.gz")

input_cnv <- function(file_path = NULL) {
  if (is.null(file_path)) {
    # Use the system file if no path is provided
    file_path <- system.file("extdata", "CNV_data_COAD.txt.gz", package = "lncRNACNVIntegrateR")
  }

  if (!file.exists(file_path)) {
    stop("The file does not exist: ", file_path)
  }

  if (grepl("\\.gz$", file_path)) {
    con <- gzfile(file_path, "rt")
    cnv_df <- read.table(con, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
    close(con)
  } else {
    cnv_df <- read.table(file_path, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
  }

  return(cnv_df)  # Return the CNV data frame
}
