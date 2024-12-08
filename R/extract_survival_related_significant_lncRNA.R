#' extract_survival_related_significant_lncRNA
#'
#' @param preprocessing_and_preparation_result list returned by preprocessing_and_preparation function
#'
#' @param Correlation_lncRNAs_and_CNVs_results list returned by Correlation_lncRNAs_and_CNVs
#'
#' @param lncRNA_expression data frame of expression values
#'
#' @param clinical_data data frame of absolute CNV values
#'
#' @param New_matrix_with_CNV_counts data frame of lncRNAs corresponding CNV details
#'
#' @return result_surv
#'
#' @import RegParallel
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#' # Example 1: using the results from previous functions:
#' # Assuming preprocessing_and_preparation_result and Correlation_lncRNAs_and_CNVs_results are obtained from previous steps
#' result_surv <- extract_survival_related_significant_lncRNA(preprocessing_and_preparation_result, Correlation_lncRNAs_and_CNVs_results)
#'
#' # Example 2: providing the data frames directly:
#' # Assuming lncRNA_expression_df, clinical_data_df, and New_matrix_with_CNV_counts_df are available
#' result_surv <- extract_survival_related_significant_lncRNA(lncRNA_expression, clinical_data, New_matrix_with_CNV_counts)
extract_survival_related_significant_lncRNA <- function(preprocessing_and_preparation_result = NULL, Correlation_lncRNAs_and_CNVs_results = NULL, lncRNA_expression = NULL, clinical_data = NULL, New_matrix_with_CNV_counts = NULL) {
  if (!is.null(Correlation_lncRNAs_and_CNVs_results)) {

    lncRNA_expression <- Correlation_lncRNAs_and_CNVs_results$final_expre_data
    New_matrix_with_CNV_counts <- Correlation_lncRNAs_and_CNVs_results$final_matrix
    clin_selected <- preprocessing_and_preparation_result$clin_selected
  } else if (!is.null(lncRNA_expressions) && !is.null(New_matrix_with_CNV_counts && !is.null(clin_selected))) {
  } else {
    stop("Either preprocessing_and_preparation_result and Correlation_lncRNAs_and_CNVs_results or lncRNA_expression, New_matrix_with_CNV_counts, and clinical_data must be provided.")  }
  expression_data.t <- t(lncRNA_expression)
  merged_df_fetched_lncRNAs_expression <- merge(expression_data.t, clin_selected, by = "row.names")
  rownames(merged_df_fetched_lncRNAs_expression) <- merged_df_fetched_lncRNAs_expression[, 1]
  merged_df_fetched_lncRNAs_expression <- merged_df_fetched_lncRNAs_expression[, -1]

  # Process vital_status column
  merged_df_fetched_lncRNAs_expression$vital_status <- ifelse(merged_df_fetched_lncRNAs_expression$vital_status == "Dead", 0, 1)
  #merged_df_fetched_lncRNAs_expression$biochemical_recurrence <- ifelse(merged_df_fetched_lncRNAs_expression$biochemical_recurrence == "NO", 0, 1)


  merged_df_fetched_lncRNAs_expression <- as.data.frame(merged_df_fetched_lncRNAs_expression)
  # merged_df_fetched_lncRNAs_expression <- merged_df_fetched_lncRNAs_expression[!(merged_df_fetched_lncRNAs_expression$days_to_last_followup == "[Not Available]"), ]
  merged_df_fetched_lncRNAs_expression$days_to_last_followup <- as.numeric(merged_df_fetched_lncRNAs_expression$days_to_last_followup)
  #merged_df_fetched_lncRNAs_expression <- merged_df_fetched_lncRNAs_expression[, !names(merged_df_fetched_lncRNAs_expression) %in% c("vital_status")]
  merged_df_fetched_lncRNAs_expression$days_to_last_followup <- as.numeric(merged_df_fetched_lncRNAs_expression$days_to_last_followup)
  merged_df_fetched_lncRNAs_expression[is.na(merged_df_fetched_lncRNAs_expression)] <- 0
  variables <- colnames(merged_df_fetched_lncRNAs_expression)[!colnames(merged_df_fetched_lncRNAs_expression) %in% c("days_to_last_followup", "vital_status")]

  res <- RegParallel(
    data = merged_df_fetched_lncRNAs_expression,
    formula = 'Surv(days_to_last_followup, vital_status) ~ [*]',
    FUN = function(formula, data)
      coxph(formula = formula,
            data = data,
            ties = 'breslow',
            singular.ok = TRUE),
    FUNtype = 'coxph',
    variables = colnames(merged_df_fetched_lncRNAs_expression)[1:ncol(merged_df_fetched_lncRNAs_expression)],
    blocksize = 30,
    cores = 1,
    nestedParallel = FALSE,
    conflevel = 95)

  # Extract significant genes related to survival
  res_sig <- res[res$P <= 0.05, ]
  res_sig <- res_sig[-nrow(res_sig), ]
  ### extract the survival related lncRNA data from the final corr data matrix ###
  lncRNA_ids_survival_related_sign <- as.data.frame(res_sig$Variable)
  #View(lncRNA_ids_survival_related_sign)
  colnames(lncRNA_ids_survival_related_sign) <- "Names"
  rownames(lncRNA_ids_survival_related_sign) <- lncRNA_ids_survival_related_sign[, 1]

  ### Now fetch these lncRNAs CNV data from our final generated table with all information ###
  merged_df_neww_with_sign_lncRNA_surv_related <- merge(lncRNA_ids_survival_related_sign, New_matrix_with_CNV_counts, by = "row.names")
  rownames(merged_df_neww_with_sign_lncRNA_surv_related) <- merged_df_neww_with_sign_lncRNA_surv_related[, 1]
  merged_df_neww_with_sign_lncRNA_surv_related <- merged_df_neww_with_sign_lncRNA_surv_related[, -(1:2)]

  ### sort on the basis of correlation and CNV frequency ###
  merged_df_neww_with_sign_lncRNA_surv_related_1 <- merged_df_neww_with_sign_lncRNA_surv_related[with(merged_df_neww_with_sign_lncRNA_surv_related, order(-corr, -CNV_frequency)), ]

  ### also extract the survival related lncRNAs expression ###

  result_df_Expr_lncRNAs_surv_related <- merge(lncRNA_ids_survival_related_sign, lncRNA_expression, by = "row.names")
  rownames(result_df_Expr_lncRNAs_surv_related) <- result_df_Expr_lncRNAs_surv_related[,1]
  result_df_Expr_lncRNAs_surv_related <- result_df_Expr_lncRNAs_surv_related[,-(1:2)]
  # Return both the data frame with survival-related lncRNAs and the list of lncRNA IDs
  return(list(data = merged_df_neww_with_sign_lncRNA_surv_related_1, ids = lncRNA_ids_survival_related_sign, final_expr = result_df_Expr_lncRNAs_surv_related, merged_expr_clinical = merged_df_fetched_lncRNAs_expression))
  #return(merged_df_neww_with_sign_lncRNA_surv_related_1)
  #return(lncRNA_ids_survival_related_sign)
}


