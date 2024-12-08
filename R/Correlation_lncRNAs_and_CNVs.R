#' Finding correlation between lncRNA expression and CNVs
#'
#' @param preprocessing_and_preparation_result list returned by preprocessing_and_preparation function
#'
#' @param lncRNA_expressions data frame of expression values
#'
#' @param lncRNA_CNV_profiles data frame of absolute CNV values
#'
#' @return Correlation_lncRNAs_and_CNVs_results
#'
#' @import devtools
#'
#' @import basicPlotteR
#'
#' @export
#'
#' @examples
#' # Example 1: Using the result from the preprocessing_and_preparation function
#' Correlation_lncRNAs_and_CNVs_results <- Correlation_lncRNAs_and_CNVs(preprocessing_and_preparation_result)
#' # Example 2: Providing lncRNA expression and CNV profiles directly
#' Correlation_lncRNAs_and_CNVs_results <- Correlation_lncRNAs_and_CNVs(lncRNA_expressions, lncRNA_CNV_profiles)
Correlation_lncRNAs_and_CNVs <- function(preprocessing_and_preparation_result = NULL, lncRNA_expressions = NULL, lncRNA_CNV_profiles = NULL) {

  if (!is.null(preprocessing_and_preparation_result)) {
    # Extract the data frames from the result
    lncRNA_expre_df <- preprocessing_and_preparation_result$result_df_Expr_lncRNAs
    lncRNA_CNV_df <- preprocessing_and_preparation_result$result_df_Expr_lncRNAs_with_CNV
  } else if (!is.null(lncRNA_expressions) && !is.null(lncRNA_CNV_profiles)) {
    lncRNA_expre_df <- lncRNA_expressions
    lncRNA_CNV_df <- lncRNA_CNV_profiles
  } else {
    stop("Either preprocessing_and_preparation_result or both lncRNA_expressions and lncRNA_CNV_profiles must be provided.")  }

  df1.t <- as.data.frame(t(lncRNA_expre_df))
  df2.t <- as.data.frame(t(lncRNA_CNV_df))

  # Sort the data frames
  df1_sorted <- df1.t[order(rownames(df1.t)), order(names(df1.t))]
  colnames(df1_sorted) <- gsub("-", ".", colnames(df1_sorted))
  df2_sorted <- df2.t[order(rownames(df2.t)), order(names(df2.t))]

  df2_subset <- as.data.frame(lapply(df2_sorted[1:nrow(df2_sorted), 1:ncol(df2_sorted)], as.numeric))
  rownames(df2_subset) <- rownames(df2_sorted)

  all.pcc=c()
  for(i in 1:nrow(df2_subset)){
    all.pcc=c(all.pcc,cor(as.numeric(df2_subset[i,]),as.numeric(df1_sorted[i,]))[1])
  }
  all.pcc.rand=c()
  rand.cnv.inds=sample(1:nrow(df2_subset),nrow(df2_subset))
  rand.lnc.inds=sample(1:nrow(df1_sorted),nrow(df1_sorted))
  for(i in 1:nrow(df2_subset)){
    all.pcc.rand=c(all.pcc.rand,cor(as.numeric(df2_subset[rand.cnv.inds,][i,]),as.numeric(df1_sorted[rand.lnc.inds,][i,]))[1])
  }
  # Plot the histogram for random correlations
  # Start by saving the plot to a PNG file
  png("correlation_distribution_plot.png", width = 10, height = 6, units = "in", res = 300)

  # Combine both actual and random correlations for checking
  all_correlations <- c(all.pcc, all.pcc.rand)

  # Calculate the proportion of positive correlations
  positive_correlations <- sum(all_correlations > 0, na.rm = TRUE)
  total_correlations <- sum(!is.na(all_correlations))
  proportion_positive <- positive_correlations / total_correlations

  # Check if the proportion of positive correlations is greater than 0.5
  if (proportion_positive > 0.5) {
    cat("The correlation distribution is predominantly positive.\n")

    # Plot histogram for random correlations
    hist(as.numeric(all.pcc.rand[!is.na(all.pcc.rand)]),
         col='grey', xlim=c(-1,1), ylim=c(0,100),
         border='grey', xlab='Pearson Correlation Coefficient', ylab='Frequency', main='Correlation Distribution')

    # Add main title for random distribution
    mtext(text = "Random", side = 3, line = 20)

    # Overlay density plot for random correlations
    par(new=T)
    plot(density(all.pcc.rand[!is.na(all.pcc.rand)]),
         xlim=c(-1,1), col='black', xlab='', ylab='', axes=FALSE, main='')

    # Overlay histogram for actual lncRNA-CNV correlations
    par(new=T)
    hist(as.numeric(all.pcc[!is.na(all.pcc)]),
         col=rgb(255, 0, 0, 80, maxColorValue=255),
         xlim=c(-1,1), ylim=c(0,100),
         border=NA, xlab='', ylab='', main='')

    # Add main title for actual correlation distribution
    mtext(text = "LncRNA-CNV", side = 2, line = 20)

    # Overlay density plot for actual lncRNA-CNV correlations
    par(new=T)
    plot(density(all.pcc[!is.na(all.pcc)]),
         xlim=c(-1,1), col='blue', xlab='', ylab='', axes=FALSE, main='')

    # Calculate means for both distributions
    mean_random <- mean(all.pcc.rand, na.rm = TRUE)
    mean_actual <- mean(all.pcc, na.rm = TRUE)

    # Add vertical lines for mean values
    abline(v = mean_random, col = 'black', lwd = 2, lty = 2) # Mean of random correlations
    abline(v = mean_actual, col = 'blue', lwd = 2, lty = 2) # Mean of actual correlations

    # Add text to indicate which line is which
    legend("topright", legend = c("Random Mean", "Actual Mean"), col = c("black", "blue"), lty = 2, lwd = 2)

    # Perform a t-test comparing actual vs random correlations
    t_test_result <- t.test(all.pcc, all.pcc.rand, alternative = "two.sided", na.rm = TRUE)

    # Extract p-value from t-test result
    p_value <- t_test_result$p.value

    # Add text annotation for p-value below the plot
    text(x = 0.5, y = -0.1, labels = paste("p-value:", round(p_value, 4)),
         col = "black", hjust = "left")

    # Create a correlation data frame
    corr_data <- as.data.frame(all.pcc)
    rownames(corr_data) <- rownames(df2_subset)
    colnames(corr_data) <- "corr"

  } else {
    cat("The correlation distribution is not predominantly positive. Skipping further analysis.\n")
  }

  # Close the PNG device
  dev.off()

  Final_matrix_CNV_plus_corr <- cbind(df2_subset, corr_data)
  New_matrix_with_only_positive_values <- Final_matrix_CNV_plus_corr[Final_matrix_CNV_plus_corr$corr > 0, ]
  New_matrix_with_only_positive_values_1 <- na.omit(New_matrix_with_only_positive_values)

  # Remove correlation column
  New_matrix_with_only_positive_values_1 <- New_matrix_with_only_positive_values_1[, -ncol(New_matrix_with_only_positive_values_1)]
  # Calculate amp/del ratios and frequency of CNVs
  max_value <- max(as.matrix(New_matrix_with_only_positive_values_1))
  min_value <- min(as.matrix(New_matrix_with_only_positive_values_1))
  count_numbers <- function(row) {
    counts <- table(row)
    names(counts) <- paste(names(counts), "count", sep = "_")
    return(counts)
  }
  counts_list <- apply(New_matrix_with_only_positive_values_1, 1, count_numbers)
  unique_values <- unique(unlist(New_matrix_with_only_positive_values_1))
  counts_df <- data.frame(matrix(0, nrow = nrow(New_matrix_with_only_positive_values_1), ncol = length(unique_values)))
  colnames(counts_df) <- unique_values
  for (i in 1:length(counts_list)) {
    counts_df[i, names(counts_list[[i]])] <- counts_list[[i]]
  }
  counts_df <- counts_df[, -(1:5)]
  counts_df[is.na(counts_df)] <- 0

  New_matrix_with_counts <- cbind(New_matrix_with_only_positive_values_1, counts_df)
  New_matrix_with_counts$sum_del <- rowSums(New_matrix_with_counts[, c("-2_count", "-1_count")])
  New_matrix_with_counts$sum_amp <- rowSums(New_matrix_with_counts[, c("1_count", "2_count")])
  New_matrix_with_counts$ratio_amp_del <- New_matrix_with_counts$sum_amp / New_matrix_with_counts$sum_del
  New_matrix_with_counts$CNV_frequency <- rowSums(New_matrix_with_counts[, c("sum_amp", "sum_del")]) / ncol(New_matrix_with_only_positive_values_1) * 100

  # Add the correlation details
  New_matrix_with_counts <- merge(New_matrix_with_counts, corr_data, by = "row.names")

  # Sort the data frame based on correlation and frequency columns
  New_matrix_with_counts <- New_matrix_with_counts[with(New_matrix_with_counts, order(-corr, -CNV_frequency)), ]
  rownames(New_matrix_with_counts) <- New_matrix_with_counts[, 1]
  New_matrix_with_counts <- New_matrix_with_counts[, -1]

  #### at last fetch these lncRNAs expression data from total lncRNA data in previous function ####
  fetch_lnc_names <- rownames(New_matrix_with_counts)
  final_lncRNAs_after_corr_list_expre <- df1.t[fetch_lnc_names, ]
  #write.table(final_lncRNAs_after_corr_list_expre, "final_expression_data.txt", sep = "\t", quote = FALSE, row.names = TRUE)
  ### adding the circos part ####
  new_matrix <- New_matrix_with_counts[, -seq(ncol(New_matrix_with_counts) - 9, ncol(New_matrix_with_counts))]
  ### for deletion ###
  new_matrix$Data <- apply(new_matrix, 1, function(x) {
    sum(x == -1) / length(x)
  })
  new_matrix_del_info <- new_matrix
  ### for amplification ###
  new_matrix <- new_matrix[, -ncol(new_matrix)]
  new_matrix$Data <- apply(new_matrix, 1, function(x) {
    sum(x == 1) / length(x)
  })
  new_matrix_amp_info <- new_matrix
  new_matrix_amp_info_1 <- data.frame(Data = new_matrix_amp_info[, ncol(new_matrix_amp_info)])
  rownames(new_matrix_amp_info_1) <- rownames(new_matrix_amp_info)

  new_matrix_del_info_1 <- data.frame(Data = new_matrix_del_info[, ncol(new_matrix_del_info)])
  rownames(new_matrix_del_info_1) <- rownames(new_matrix_del_info)

  ### Assuming the lncRNA positions already stored from the input_lncRNA function so before running this function youn need to run the input_lncRNA_positions function###
  #lncRNA_positions <- input_lncRNA_positions()
  #rownames(lncRNA_positions) <- lncRNA_positions[,1]
  #lncRNA_positions <- lncRNA_positions[, -1]
  lnc.locs.del <- merge(lncRNA_positions, new_matrix_del_info_1, by = "row.names")
  rownames(lnc.locs.del) <- lnc.locs.del[, 1]
  lnc.locs.del <- lnc.locs.del[, -1]
  lnc.locs.amp <- merge(lncRNA_positions, new_matrix_amp_info_1, by = "row.names")
  rownames(lnc.locs.amp) <- lnc.locs.amp[, 1]
  lnc.locs.amp <- lnc.locs.amp[, -1]
  ### row run the circos code ###
  library(RCircos)
  data(UCSC.HG38.Human.CytoBandIdeogram)
  data(RCircos.Histogram.Data)
  RCircos.Set.Core.Components(cyto.info = UCSC.HG38.Human.CytoBandIdeogram, chr.exclude = NULL, tracks.inside = 5, tracks.outside = 0)
  png("lncRNA_distribution_across_genome_plot.png", width = 8, height = 6, units = "in", res = 300)
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()
  # Assuming lnc.locs.amp and lnc.locs.del are properly formatted data frames

  RCircos.Histogram.Plot(lnc.locs.amp, data.col = 4, track.num = 2, side = 'in', max.value = 1)
  cat("Inner Circle represents the amplification...\n")
  RCircos.Histogram.Plot(lnc.locs.del, data.col = 4, track.num = 4, max.value = 1)
  cat("IThe outer one represents the deletion...\n")

  # Save the plot to a PNG file
  dev.off()
  #write.table(final_lncRNAs_after_corr_list_expre, "final_expression_data.txt", sep = "\t", quote = FALSE, row.names = TRUE)
  Correlation_lncRNAs_and_CNVs_results <- list(
    t_test_result = t_test_result,
    final_matrix = New_matrix_with_counts,
    final_expre_data = final_lncRNAs_after_corr_list_expre
  )
  return(Correlation_lncRNAs_and_CNVs_results)
}

