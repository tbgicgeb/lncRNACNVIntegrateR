% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/Correlation_with_PCGs_and_functional_enrichment.R
\name{Correlation_with_PCGs_and_functional_enrichment}
\alias{Correlation_with_PCGs_and_functional_enrichment}
\title{Prognostic lncRNA signatures functional enrichment analysis}
\usage{
Correlation_with_PCGs_and_functional_enrichment(
  preprocessing_and_preparation_result = NULL,
  result_riskscore_model = NULL,
  PCG_expression = NULL,
  lncRNA_expression = NULL
)
}
\arguments{
\item{preprocessing_and_preparation_result}{list returned by preprocessing_and_preparation function}

\item{result_riskscore_model}{resulted dataframe from Risk_score_model_development function}

\item{PCG_expression}{data frame of lncRNAs corresponding CNV details}

\item{lncRNA_expression}{data frame of expression values of identified final lncRNA signtaures}
}
\value{
result_enrichment
}
\description{
Prognostic lncRNA signatures functional enrichment analysis
}
\examples{
# Example 1: using the results from previous functions:
# Assuming preprocessing_and_preparation_result and result_riskscore_model are obtained from previous functions
result_enrichment <- Correlation_with_PCGs_and_functional_enrichment(preprocessing_and_preparation_result, result_riskscore_model)

# Example 2: providing the data frames directly:
# Assuming PCG_expression_df and lncRNA_expression_df are available
result_enrichment <- Correlation_with_PCGs_and_functional_enrichment(PCG_expression, lncRNA_expression)
}
