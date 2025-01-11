#' Preparing input files for expression data
#'
#' @param file_path Optional character string. If provided, a custom path to a tab-delimited
#'                  expression file will be used. If not provided, the default system file
#'                  will be used from the package.
#' @return A data frame of expression values
#'
#' @import utils
#' @export
#' @examples
#' expr_df <- input_expr()  # Uses the default system file
#' expr_df <- input_expr("/custom/path/to/gene_expression.txt.gz")
input_expr <- function(file_path = NULL) {
  if (is.null(file_path)) {
    # Use the system file if no path is provided
    file_path <- system.file("extdata", "gene_expression_COAD.txt.gz", package = "lncRNACNVIntegrateR")
  }

  if (!file.exists(file_path)) {
    stop("The file does not exist: ", file_path)
  }

  if (grepl("\\.gz$", file_path)) {
    con <- gzfile(file_path, "rt")
    expr_df <- read.table(con, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    close(con)
  } else {
    expr_df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }

  return(expr_df)  # Return the expression data frame
}
