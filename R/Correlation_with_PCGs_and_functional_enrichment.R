#' Prognostic lncRNA signatures functional enrichment analysis
#'
#' @param preprocessing_and_preparation_result list returned by preprocessing_and_preparation function
#'
#' @param result_riskscore_model resulted dataframe from Risk_score_model_development function
#'
#' @param lncRNA_expression data frame of expression values of identified final lncRNA signtaures
#'
#' @param PCG_expression data frame of lncRNAs corresponding CNV details
#'
#' @return result_enrichment
#'
#' @import ggplot2
#'
#' @import dplyr
#'
#' @import enrichR
#'
#' @importFrom enrichR enrichr listEnrichrSites setEnrichrSite listEnrichrDbs
#'
#' @export
#'
#' @examples
#' # Example 1: using the results from previous functions:
#' # Assuming preprocessing_and_preparation_result and result_riskscore_model are obtained from previous functions
#' result_enrichment <- Correlation_with_PCGs_and_functional_enrichment(preprocessing_and_preparation_result, result_riskscore_model)
#'
#' # Example 2: providing the data frames directly:
#' # Assuming PCG_expression_df and lncRNA_expression_df are available
#' result_enrichment <- Correlation_with_PCGs_and_functional_enrichment(PCG_expression, lncRNA_expression)
Correlation_with_PCGs_and_functional_enrichment <- function(preprocessing_and_preparation_result = NULL, result_riskscore_model = NULL, PCG_expression = NULL, lncRNA_expression = NULL) {
  if (!is.null(preprocessing_and_preparation_result) && !is.null(result_riskscore_model)) {
    PCG_expression <- preprocessing_and_preparation_result$PCG_matrix
    lncRNA_expression <- result_riskscore_model$final_expre_matrix_prognostic_markers
  } else if (!is.null(PCG_expression) && !is.null(lncRNA_expression)) {
  } else {
    stop("Either preprocessing_and_preparation_result & result_riskscore_model or PCG_expression & lncRNA_expression must be provided.")  }
  result_df_Expr_PCGs <- PCG_expression
  prognostic_markers_expression <- lncRNA_expression
  columns_to_exclude <- c("vital_status", "days_to_last_followup")

  # Exclude columns by names
  prognostic_markers_expression <- prognostic_markers_expression[, -which(names(prognostic_markers_expression) %in% columns_to_exclude)]
  prognostic_markers_expression.t <- t(prognostic_markers_expression)

  # Get row names
  row_names <- rownames(prognostic_markers_expression.t)

  # Create a list to store individual data frames
  correlation_dfs <- list()

  # Iterate over each row and calculate correlation
  for (row_name in row_names) {
    row_data <- prognostic_markers_expression.t[row_name, , drop = FALSE]
    row_df <- data.frame(row_data)
    row_df$rowname <- row_name  # Add a column for row name

    # Transpose the PCG matrix and order by colnames
    result_df_Expr_PCGs.t <- t(result_df_Expr_PCGs)
    colnames(result_df_Expr_PCGs.t) <- gsub("-", ".", colnames(result_df_Expr_PCGs.t))
    result_df_Expr_PCGs.t <- result_df_Expr_PCGs.t[, order(colnames(result_df_Expr_PCGs.t))]

    # Calculate correlation
    correlation_values <- cor(t(row_data), t(result_df_Expr_PCGs.t))

    # Create a data frame with correlation values
    correlation_df <- data.frame(lncRNA = row_name, correlation_values)

    # Append to the list
    correlation_dfs[[row_name]] <- correlation_df
  }

  # Remove the first column (lncRNA) from each data frame
  correlation_dfs <- lapply(correlation_dfs, function(df) df[, -1, drop = FALSE])

  for (i in seq_along(correlation_dfs)) {
    correlation_dfs[[i]] <- t(correlation_dfs[[i]])
  }

  correlation_dfs <- lapply(correlation_dfs, function(df) df[order(-df[, 1]), ])

  # Function to get top 50 genes for a given lncRNA
  get_top_50_genes <- function(correlation_values) {
    return(head(correlation_values, 50))
  }

  # Initialize a list to store top 50 genes for each lncRNA
  top_50_genes_list <- list()

  # Loop through each lncRNA and store top 50 genes in the list
  for (lncRNA in names(correlation_dfs)) {
    top_50_genes <- get_top_50_genes(correlation_dfs[[lncRNA]])
    top_50_genes_list[[lncRNA]] <- top_50_genes
  }

  gene_names_list <- list()

  # Get the current working directory
  current_directory <- getwd()

  # Specify the output directory within the current working directory
  output_directory <- file.path(current_directory, "enrichment_result")

  # Create the output directory
  dir.create(output_directory, recursive = TRUE)

  # Loop through each lncRNA and store gene names in the list
  for (lncRNA in names(top_50_genes_list)) {
    gene_names <- names(top_50_genes_list[[lncRNA]])
    gene_names_list[[lncRNA]] <- gene_names
  }

  #### 1. if want a combined Gene ontology enrichment for all lncRNAs ####
  # Display the result
  # print(gene_names_list)
 library(enrichR)
  # Now use enrichr
  listEnrichrSites()
  setEnrichrSite("Enrichr")
  websiteLive <- TRUE
  dbs <- listEnrichrDbs()
  if (is.null(dbs)) websiteLive <- FALSE
  if (websiteLive) head(dbs)

  dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021", "KEGG_2021_Human")

  combined_list <- unlist(gene_names_list)
  combined_list
  enriched <- enrichr(combined_list, dbs)
  # Set dimensions for the plot
  plot_width <- 8  # Replace with the desired width in inches
  plot_height <- 6
  for (i in seq_along(enriched)) {
    enrichR::plotEnrich(enriched[[i]], showTerms = 25, numChar = 80, y = "Count", orderBy = "P.value")
    ggplot2::ggsave(paste0(output_directory, "/", "plot_",  names(enriched)[i], ".png"), device = "png", width = plot_width, height = plot_height)
    write.table(enriched[[i]], file = paste0(output_directory, "/", names(enriched)[i], ".txt"), sep = "\t", quote = FALSE)
  }
  # Return the results
  return(list(
    correlation_dfs = correlation_dfs,
    gene_names_list = gene_names_list
  ))
}
