#' Reading the lncRNA positions in the human genome
#'
#' @param file_path Optional character string. If provided, a custom path to the
#'                  lncRNA positions file will be used. If not provided, the default
#'                  system file will be loaded.
#' @return A data frame of lncRNA positions
#'
#' @import utils
#' @export
#' @examples
#' lncRNA_positions <- input_lncRNA_positions()
#' lncRNA_positions <- input_lncRNA_positions("/custom/path/to/lncRNA_positions.txt")
input_lncRNA_positions <- function(file_path = NULL) {
  if (is.null(file_path)) {
    # Use the system file if no path is provided
    file_path <- system.file("extdata", "all_lncRNAs_positions.txt", package = "lncRNACNVIntegrateR")
  }

  if (!file.exists(file_path)) {
    stop("The file does not exist: ", file_path)
  }

  lncRNA_positions <- read.table(file_path, header = TRUE, sep = "\t", row.names = 1)
  return(lncRNA_positions)
}
