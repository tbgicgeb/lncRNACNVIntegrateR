#' Preprocessing of input files
#'
#' @param expr_df data frame of expression values
#'
#' @param cnv_df data frame of absolute CNV values
#'
#' @param clin_df data frame of clinical values
#'
#' @return preprocessing_and_preparation_result
#'
#' @import DESeq2
#'
#' @import dplyr
#'
#' @import SummarizedExperiment
#'
#' @import R.utils
#'
#' @export
#'
#' @examples
#' preprocessing_and_preparation_result <- preprocessing_and_preparation(expr_df, cnv_df, clin_df)
preprocessing_and_preparation <- function(expr_df, cnv_df, clin_df) {

  setwd(getwd())

  ### Identify if the dataset contains TCGA IDs ###
  # Extract the IDs (row names) from the expression data
  expr_ids <- rownames(expr_df)
  tcga_expr_check <- all(grepl("^TCGA", expr_ids))
  cnv_ids <- rownames(cnv_df)
  tcga_cnv_check <- all(grepl("^TCGA", cnv_ids))
  ## for clinical data ###
  clinical_data <- clin_df
  colnames(clinical_data) <- clinical_data[1, ]
  clinical_data <- clinical_data[-(1:2), ]
  selected_clin <- as.data.frame(clinical_data[, c("bcr_patient_barcode", "vital_status", "days_to_last_followup")])
  rownames(selected_clin) <- selected_clin[, 1]
  selected_clin <- selected_clin[, -1]
  clin_ids <- rownames(selected_clin)
  tcga_clin_check <- all(grepl("^TCGA", clin_ids))

  ### Preprocessing based on ID check ###
  if (tcga_expr_check && tcga_cnv_check && tcga_clin_check) {
    cat("TCGA IDs detected. Proceeding with TCGA-specific data preprocessing...\n")

    ### Processing Expression Data for TCGA ###
    expression_data <- as.data.frame(expr_df)
    shortened_row_names <- substr(rownames(expression_data), 1, 12)
    duplicated_shortened <- duplicated(shortened_row_names)
    expression_data <- expression_data[!duplicated_shortened, ]

    # After removing duplicates, shorten the row names again to match the remaining rows
    shortened_row_names <- substr(rownames(expression_data), 1, 12)
    rownames(expression_data) <- shortened_row_names

    ### Processing CNV Data for TCGA ###
    cnv_data <- as.data.frame(cnv_df)
    rownames(cnv_data) <- substr(rownames(cnv_data), 1, 12)
    rownames(cnv_data) <- gsub("\\.", "-", rownames(cnv_data))

    ### Processing Clinical Data for TCGA ###
    ### already done intially.

  } else {
    cat("Non-TCGA IDs detected. Proceeding with generic data preprocessing...\n")

    ### Generic Processing for Non-TCGA ###
    expression_data <- as.data.frame(expr_df)
    cnv_data <- as.data.frame(cnv_df)
    clinical_data <- clin_df
    selected_clin <- as.data.frame(clinical_data[, c("IDs", "vital_status", "days_to_last_followup")])
    rownames(selected_clin) <- selected_clin[, 1]
    selected_clin <- selected_clin[, -1]
  }

  ### Common processing for both TCGA and non-TCGA datasets ###

  ### Find Common Samples between Expression, CNV, and Clinical Data ###
  common_samples <- intersect(rownames(cnv_data), rownames(expression_data))

  if (length(common_samples) == 0) {
    stop("No common samples found between CNV and expression data.")
  }

  clin_selected <- selected_clin[common_samples, ]
  clin_selected <- na.omit(clin_selected)
  raw_names_clin_data_avail <- rownames(clin_selected)

  ### Select expression and CNV data for final common samples ###
  df_gene_expr_data.t <- expression_data[raw_names_clin_data_avail, ]
  CNV_data_transposed <- cnv_data[raw_names_clin_data_avail, ]

  # Create metadata for DESeq2
  df_gene_expr_data.t2 <- t(df_gene_expr_data.t)
  column2 <- as.data.frame(colnames(df_gene_expr_data.t2))
  condn2 <- as.data.frame(rep("Case", times = nrow(column2)))
  metaData2 <- cbind(column2, condn2)
  colnames(metaData2)[1:2] <- c("Sample", "condition")
  cat("Starting with expression data preprocessing...\n")

  # Create a DESeqDataSet object
  dds <- DESeq2::DESeqDataSetFromMatrix(df_gene_expr_data.t2, metaData2, design = ~1)
  # Run DESeq2
  ddsDE <- DESeq2::DESeq(dds)
  # Apply Variance Stabilizing Transformation
  vst_normal <- DESeq2::varianceStabilizingTransformation(ddsDE, blind = TRUE)
  vst_normal_results <- SummarizedExperiment::assay(vst_normal)
  ## firstly download the gtf file from gencode ####
  # URL to the GTF file (GENCODE version 22)
  cat("Downloading the gtf file from gencode...\n")

  gtf_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.annotation.gtf.gz"

  # Destination file path where you want to save the downloaded file
  current_working_directory <- getwd()

  # Set the destination path to the current working directory
  destination_path <- file.path(current_working_directory, "gencode.v22.annotation.gtf.gz")

  options(timeout = 600)

  # Download the GTF file with a longer timeout
  download.file(url = gtf_url, destfile = destination_path, method = "auto", mode="wb")

  # Uncompress the downloaded file if it's in .gz format
  if (grepl(".gz$", destination_path)) {
    # Using R's built-in gunzip function which works on all platforms (cross-platform)
    gzfile_path <- sub(".gz$", "", destination_path)

    # Uncompress the file using R's readBin() function and write it out
    R.utils::gunzip(destination_path, destname = gzfile_path)

    # Update destination_path to the uncompressed file
    destination_path <- gzfile_path
  }

  # Check if the download was successful and file exists
  if (file.exists(destination_path)) {
    cat("Downloaded GTF file successfully.\n")

    # Read the downloaded GTF file into a data frame
    gtf <- read.table(destination_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

    # Check if the "V9" column exists
    if ("V9" %in% colnames(gtf)) {
      # Split the "V9" column into separate components
      split_names <- strsplit(as.character(gtf$V9), ';')

      # Merge the resulting list into a data frame
      merged_df <- do.call(rbind, split_names)

      # Combine the merged columns with the original data frame
      df_merged <- cbind(gtf, merged_df)

      # Optionally, check the merged data
      head(df_merged[1:5, 1:5])
    } else {
      cat("Column V9 not found in the GTF file. Please check the file format.\n")
    }
  } else {
    cat("Download failed. Please check the URL or your internet connection.\n")
  }
 cat("Extracting lncRNAs and PCGs via mapping from gtf file...\n")
  # Filter and subset lncRNA transcripts
  matched_processed_transcript <- subset(df_merged, grepl("processed_transcript", df_merged$'2', ignore.case = TRUE))
  prime_overlapping_ncrna <- subset(df_merged, grepl("3prime_overlapping_ncrna", df_merged$'2', ignore.case = TRUE))
  sense_intronic <- subset(df_merged, grepl("sense_intronic", df_merged$'2', ignore.case = TRUE))
  sense_overlapping <- subset(df_merged, grepl("sense_overlapping", df_merged$'2', ignore.case = TRUE))
  antisense <- subset(df_merged, grepl("antisense", df_merged$'2', ignore.case = TRUE))
  lncRNA <- subset(df_merged, grepl("lincRNA", df_merged$'2', ignore.case = TRUE))

  # Combine all lncRNA subsets into one data frame
  All_lnc <- rbind(matched_processed_transcript, prime_overlapping_ncrna, sense_intronic, sense_overlapping, antisense, lncRNA)

  # Extract and clean lncRNA names
  All_lnc_names <- as.data.frame(All_lnc$'4')
  All_lnc_names <- gsub("gene_name", "", All_lnc_names$`All_lnc$"4"`)
  All_lnc_names2 <- as.data.frame(All_lnc_names)
  All_lnc_names2_uniq <- unique(All_lnc_names2)
  colnames(All_lnc_names2_uniq)[colnames(All_lnc_names2_uniq) == "All_lnc_names"] <- "Name"

  ### Now fetch for PCGs ###

  Protein_coding_genes <- subset(df_merged, grepl("protein_coding", df_merged$'2', ignore.case = TRUE))
  All_PCG_names <- as.data.frame(Protein_coding_genes$'4')
  All_PCG_names <- gsub("gene_name", "", All_PCG_names$`Protein_coding_genes$"4"`)
  All_PCG_names2 <- as.data.frame(All_PCG_names)
  All_PCG_names2_uniq <- unique(All_PCG_names2)
  colnames(All_PCG_names2_uniq)[colnames(All_PCG_names2_uniq) == "All_PCG_names"] <- "Name"

  # Extract gene names from All_lnc_names2_uniq
  PCG_names <- as.data.frame(trimws(All_PCG_names2_uniq$Name))  # Remove leading/trailing whitespaces


  # Now fetch the lncRNAs expression and CNV data from common samples
  # library(dplyr)

  vst_normal_results <- as.data.frame(t(vst_normal_results))
  gene_names <- as.data.frame(trimws(All_lnc_names2_uniq$Name))  # Remove leading/trailing whitespaces
  column_names <- as.data.frame(colnames(vst_normal_results))


  # Find the intersection of gene names and column names
  intersecting_genes <- intersect(gene_names$`trimws(All_lnc_names2_uniq$Name)`, column_names$`colnames(vst_normal_results)`)
  # Now you can use the select function
  result_df_Expr_lncRNAs <- vst_normal_results %>%
    select(all_of(intersecting_genes))

  # Find the intersection of gene names and column names
  intersecting_genes_PCGs <- intersect(PCG_names$`trimws(All_PCG_names2_uniq$Name)`, column_names$`colnames(vst_normal_results)`)
  # Now you can use the select function
  result_df_Expr_PCGs <- vst_normal_results %>%
    select(all_of(intersecting_genes_PCGs))
  # write.table(result_df_Expr_PCGs, "result_df_Expr_PCGs.txt", sep="\t", quote=FALSE)

  # Extract CNV data
  common_CNV_data <- as.data.frame(CNV_data_transposed)
  column_names_CNV <- as.data.frame(colnames(CNV_data_transposed))
  intersecting_genes_CNV <- intersect(gene_names$`trimws(All_lnc_names2_uniq$Name)`, column_names_CNV$`colnames(CNV_data_transposed)`)
  result_df_CNVs <- common_CNV_data %>% select(all_of(intersecting_genes_CNV))

  # Fetch lncRNAs expression only for those which have CNV call also
  col_new_lncRNA <- as.data.frame(colnames(result_df_Expr_lncRNAs))
  common_genes <- intersect(intersecting_genes_CNV, col_new_lncRNA$`colnames(result_df_Expr_lncRNAs)`)

  # Intersect genes that have expression also from the CNV data
  result_df_Expr_lncRNAs <- result_df_Expr_lncRNAs %>%
    select(all_of(common_genes))
  ### expression final matrix lncRNA ###
  head(result_df_Expr_lncRNAs[1:10,1:10])

  # Also fetch these from CNV as the number is not the same when fetched from expression data
  result_df_Expr_lncRNAs_with_CNV <- result_df_CNVs %>%
    select(all_of(common_genes))
  ### CNVs final matrix ###
  head(result_df_Expr_lncRNAs_with_CNV[1:10,1:10])

  # Return the data frames as a list
  return(list(result_df_Expr_lncRNAs = result_df_Expr_lncRNAs,
              result_df_Expr_lncRNAs_with_CNV = result_df_Expr_lncRNAs_with_CNV,
              clin_selected = clin_selected, PCG_matrix = result_df_Expr_PCGs))
}



