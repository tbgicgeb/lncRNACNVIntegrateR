#' Visulalization of correlation between lncRNA expression and CNVs profiles
#'
#' @param result_surv list returned by extract_survival_related_significant_lncRNA
#'
#' @param lncRNA_expression_data_surv_related data frame of expression values
#'
#' @param lncRNAs_surv_related_CNV_details data frame of lncRNAs corresponding CNV details
#'
#' @return result_visualization
#'
#' @import ggplot2
#'
#' @import dplyr
#'
#' @import purrr
#'
#' @import ggpubr
#'
#' @export
#'
#' @examples
#' # Example 1: using the results from previous function:
#' # Assuming result_surv is obtained from the extract_survival_related_significant_lncRNA function
#' result_visualization <- visualization_of_lncRNA_and_CNVs_association(result_surv)
#'
#' # Example 2: providing the data frames directly:
#' # Assuming lncRNA_expression_data_surv_related_df and lncRNAs_surv_related_CNV_details_df are available
#' result_visualization <- visualization_of_lncRNA_and_CNVs_association(lncRNA_expression_data_surv_related, lncRNAs_surv_related_CNV_details)
visualization_of_lncRNA_and_CNVs_association <- function(result_surv = NULL, lncRNA_expression_data_surv_related = NULL, lncRNAs_surv_related_CNV_details = NULL) {
  if (!is.null(result_surv)) {

    lncRNA_expression_data_surv_related <- result_surv$final_expr
    lncRNAs_surv_related_CNV_details <- result_surv$data
  } else if (!is.null(lncRNA_expression_data_surv_related) && !is.null(lncRNAs_surv_related_CNV_details)) {
  } else {
    stop("Either result_surv or lncRNA_expression_data_surv_related and lncRNAs_surv_related_CNV_details must be provided.")  }

  # Extracting the expression data
  vst_data <- lncRNA_expression_data_surv_related
  colnames(vst_data) <- gsub("-", ".", colnames(vst_data))
  lncRNA_list <- rownames(lncRNAs_surv_related_CNV_details)

  # Create a directory to save separate expression files
  dir.create("lncRNA_expression_files", showWarnings = FALSE)

  # Loop through the lncRNA list and save separate expression files
  for (lncRNA_name in lncRNA_list) {
    if (lncRNA_name %in% rownames(vst_data)) {
      lncRNA_expression <- vst_data[lncRNA_name, , drop = FALSE]
      file_name <- paste("lncRNA_expression_files", "/", lncRNA_name, "_expression.txt", sep = "")
      write.table(lncRNA_expression, file_name, sep = "\t", row.names = TRUE, col.names = TRUE)
      cat(paste("Saved expression for", lncRNA_name, "to", file_name, "\n"))
    } else {
      cat(paste("lncRNA", lncRNA_name, "not found in the data\n"))
    }
  }

  # Get the current working directory
  current_directory <- getwd()

  # Specify the output directory within the current working directory
  output_directory <- file.path(current_directory, "lncRNA_files")

  # Create the output directory
  dir.create(output_directory, recursive = TRUE)

  # Step 1: Read the CNV data file
  cnv_data <- lncRNAs_surv_related_CNV_details[, grepl("^TCGA", colnames(lncRNAs_surv_related_CNV_details))]

  # Step 2: Extract TCGA sample names
  sample_names <- colnames(cnv_data)
  amplification_samples <- list()
  deletion_samples <- list()
  normal_samples <- list()

  # Step 3: Categorize samples based on CNV values for each lncRNA
  for (lncRNA_name in rownames(cnv_data)) {
    lncRNA_cnv_values <- cnv_data[lncRNA_name, ]
    amplification_samples[[lncRNA_name]] <- sample_names[lncRNA_cnv_values > 0]
    deletion_samples[[lncRNA_name]] <- sample_names[lncRNA_cnv_values < 0]
    normal_samples[[lncRNA_name]] <- sample_names[lncRNA_cnv_values == 0]
  }

  # Step 4: Store the sample names in each category
  categories <- list(
    amplification = amplification_samples,
    deletion = deletion_samples,
    normal = normal_samples
  )
  lncRNA_names <- lncRNA_list
  for (lncRNA_name in lncRNA_names) {
    for (category_name in names(categories)) {
      file_name <- paste(lncRNA_name, category_name, ".txt", sep = "_")
      file_path <- file.path(output_directory, file_name)
      writeLines(categories[[category_name]][[lncRNA_name]], con = file_path)
    }
  }

  cat("Text files have been created and saved in the directory:", output_directory, "\n")

  # Fetch expression for each category
  input_expression_directory <- "lncRNA_expression_files"
  input_cnv_directory <- "lncRNA_files"
  output_data_directory <- "lncRNA_category_data"
  dir.create(output_data_directory, showWarnings = FALSE)

  for (lncRNA_name in lncRNA_names) {
    for (category_name in c("normal", "amplification", "deletion")) {
      expression_file <- file.path(input_expression_directory, paste0(lncRNA_name, "_expression.txt"))
      category_file <- file.path(input_cnv_directory, paste0(lncRNA_name, "_", category_name, "_.txt"))

      if (file.exists(expression_file) && file.exists(category_file)) {
        expression_data <- read.table(expression_file, header = TRUE, sep = "\t", row.names = 1)
        category_samples <- readLines(category_file)
        category_expression <- expression_data[, category_samples, drop = FALSE]
        output_file <- file.path(output_data_directory, paste0(lncRNA_name, "_", category_name, "_data.txt"))
        write.table(category_expression, output_file, sep = "\t", row.names = TRUE, col.names = TRUE)
        cat(paste("Saved category data for", lncRNA_name, "(", category_name, ") to", output_file, "\n"))
      } else {
        cat(paste("Expression or category file not found for", lncRNA_name, "(", category_name, ")\n"))
      }
    }
  }

  cat("Category data files have been created and saved in the directory:", output_data_directory, "\n")

  # Merge normal and amplification, and normal and deletion
  input_data_directory <- "lncRNA_category_data"
  output_merged_directory <- "merged_category_data"
  dir.create(output_merged_directory, showWarnings = FALSE)

  for (lncRNA_name in lncRNA_names) {
    normal_file <- file.path(input_data_directory, paste0(lncRNA_name, "_normal_data.txt"))
    amplification_file <- file.path(input_data_directory, paste0(lncRNA_name, "_amplification_data.txt"))
    deletion_file <- file.path(input_data_directory, paste0(lncRNA_name, "_deletion_data.txt"))

    if (file.exists(normal_file) && file.exists(amplification_file) && file.exists(deletion_file)) {
      normal_data <- read.table(normal_file, header = TRUE, sep = "\t", row.names = 1)
      amplification_data <- read.table(amplification_file, header = TRUE, sep = "\t", row.names = 1)
      deletion_data <- read.table(deletion_file, header = TRUE, sep = "\t", row.names = 1)
      merged_normal_amplification <- cbind(normal_data, amplification_data)
      merged_normal_deletion <- cbind(normal_data, deletion_data)

      merged_normal_amplification_file <- file.path(output_merged_directory, paste0(lncRNA_name, "_normal_amplification.txt"))
      merged_normal_deletion_file <- file.path(output_merged_directory, paste0(lncRNA_name, "_normal_deletion.txt"))

      write.table(merged_normal_amplification, merged_normal_amplification_file, sep = "\t", row.names = TRUE, col.names = TRUE)
      write.table(merged_normal_deletion, merged_normal_deletion_file, sep = "\t", row.names = TRUE, col.names = TRUE)

      cat(paste("Saved merged data for", lncRNA_name, "normal and amplification to", merged_normal_amplification_file, "\n"))
      cat(paste("Saved merged data for", lncRNA_name, "normal and deletion to", merged_normal_deletion_file, "\n"))
    } else {
      cat(paste("Data files not found for", lncRNA_name, "\n"))
    }
  }

  cat("Merged data files have been created and saved in the directory:", output_merged_directory, "\n")

  # Label files generation
  input_data_directory <- "lncRNA_category_data"
  output_label_directory <- "label_files"
  dir.create(output_label_directory, showWarnings = FALSE)

  for (lncRNA_name in lncRNA_names) {
    normal_file <- file.path(input_data_directory, paste0(lncRNA_name, "_normal_data.txt"))
    amplification_file <- file.path(input_data_directory, paste0(lncRNA_name, "_amplification_data.txt"))
    deletion_file <- file.path(input_data_directory, paste0(lncRNA_name, "_deletion_data.txt"))

    if (file.exists(normal_file) && file.exists(amplification_file) && file.exists(deletion_file)) {
      normal_data <- read.table(normal_file, header = TRUE, sep = "\t", row.names = 1)
      amplification_data <- read.table(amplification_file, header = TRUE, sep = "\t", row.names = 1)
      deletion_data <- read.table(deletion_file, header = TRUE, sep = "\t", row.names = 1)

      normal_labels <- rep("N", ncol(normal_data))
      amplification_labels <- rep("A", ncol(amplification_data))
      deletion_labels <- rep("D", ncol(deletion_data))

      label_normal_amplification_file <- file.path(output_label_directory, paste0("label_", lncRNA_name, "_Normal_amplification.txt"))
      label_normal_deletion_file <- file.path(output_label_directory, paste0("label_", lncRNA_name, "_Normal_deletion.txt"))

      label_data_normal_amplification <- data.frame(Label = c(normal_labels, amplification_labels))
      write.table(label_data_normal_amplification, label_normal_amplification_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

      label_data_normal_deletion <- data.frame(Label = c(normal_labels, deletion_labels))
      write.table(label_data_normal_deletion, label_normal_deletion_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

      cat(paste("Saved labels for", lncRNA_name, "Normal_amplification to", label_normal_amplification_file, "\n"))
      cat(paste("Saved labels for", lncRNA_name, "Normal_deletion to", label_normal_deletion_file, "\n"))
    } else {
      cat(paste("Data files not found for", lncRNA_name, "\n"))
    }
  }

  cat("Label files have been created and saved in the directory:", output_label_directory, "\n")

  # Now create the violin plots to show amplification or deletion correlation with expression
  merged_data_dir <- "merged_category_data"
  label_dir <- "label_files"
  output_dir <- "combined_data"

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  data_files <- list.files(merged_data_dir, pattern = "_normal_(amplification|deletion)\\.txt")

  for (data_file in data_files) {
    gene_name <- sub("_(normal_(amplification|deletion))\\.txt$", "", data_file)
    category <- sub(".*_(normal_(amplification|deletion))\\.txt$", "\\2", data_file)

    data <- read.table(file.path(merged_data_dir, data_file), header = TRUE, row.names = 1, sep = "\t")
    transposed_data <- t(data)
    label_file <- paste0("label_", gene_name, "_Normal_", category, ".txt")
    labels <- read.table(file.path(label_dir, label_file), header = TRUE)

    combined_data <- cbind(transposed_data, labels)

    combined_data_file <- file.path(output_dir, paste0(gene_name, "_", category, "_combined_data_plus_labels.txt"))
    write.table(combined_data, combined_data_file, sep = "\t", quote = FALSE)
  }

  cat("Combined data files with labels (transposed) saved in the", output_dir, "directory.\n")

  # Create and save violin plots
  combined_data_dir <- "combined_data"
  output_dir <- "violin_plots"

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  combined_data_files <- list.files(combined_data_dir, pattern = "_combined_data_plus_labels.txt")

  generate_violin_plot_and_test <- function(data_file, output_dir) {
    combined_data <- read.table(file.path(combined_data_dir, data_file), header = TRUE, sep = "\t")
    gene_name <- sub("_(amplification|deletion)_combined_data_plus_labels.txt$", "", data_file)
    category <- sub(".*_(amplification|deletion)_combined_data_plus_labels.txt$", "\\1", data_file)
    plot <- ggplot(combined_data, aes(x = Label, y = get(gene_name), fill = Label)) +
      geom_violin() +
      geom_boxplot(width = 0.2, position = position_dodge(width = 0.75)) +
      ylab(paste(gene_name, " (Expression)")) +
      xlab("Groups") +
      ylim(min(combined_data[[gene_name]]) - 1, max(combined_data[[gene_name]]) + 1) +
      stat_compare_means(method = "kruskal.test")
    plot_file <- file.path(output_dir, paste0(gene_name, "_", category, "_violin_plot.png"))
    # Check if P-value is significant
    kruskal_result <- kruskal.test(get(gene_name) ~ Label, data = combined_data)
    if (kruskal_result$p.value <= 0.05) {
      ggsave(plot_file, plot, width = 12, height = 8)
      cat(paste("Plot saved for", gene_name, "with significant P-value in category", category, "\n"))
      # Print the plot to the R Studio window
      return(list(GeneName = gene_name, Category = category, KruskalTestResult = kruskal_result, Plot = plot))
    } else {
      cat(paste("Plot not saved for", gene_name, "as P-value is not significant in category", category, "\n"))
      return(list(GeneName = gene_name, Category = category, KruskalTestResult = kruskal_result, Plot = NULL))
    }
  }

  # Assume combined_data_files is a list of files and output_dir is the output directory
  significant_plot <- NULL
  for (file in combined_data_files) {
    result <- generate_violin_plot_and_test(file, output_dir)
    gene_name <- sub("_(amplification|deletion)_combined_data_plus_labels.txt$", "", data_file)
    if (result$KruskalTestResult$p.value <= 0.05) {
      significant_plot <- result$Plot
      cat(paste("Plot saved for", gene_name, "with significant P-value in category", category, "\n"))
    }
  }
  # Print the plot for the last significant gene, if any
  if (!is.null(significant_plot)) {
    print(significant_plot)
  }

  kruskal_results_list <- map(combined_data_files, ~generate_violin_plot_and_test(.x, output_dir))

  kruskal_results_df <- do.call(rbind, lapply(kruskal_results_list, function(result) {
    data.frame(
      GeneName = result$GeneName,
      Category = result$Category,
      PValue = result$KruskalTestResult$p.value
    )
  }))

  new_df_significant <- kruskal_results_df %>% filter(PValue <= 0.05)

  return(new_df_significant)
}
