% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/groupComparisonLiP.R
\name{groupComparisonLiP}
\alias{groupComparisonLiP}
\title{Model LiP and TrP data and make adjustments if needed
Returns list of three modeled datasets}
\usage{
groupComparisonLiP(
  data,
  contrast.matrix = "pairwise",
  fasta.path = NULL,
  log_base = 2,
  use_log_file = TRUE,
  append = FALSE,
  verbose = TRUE,
  log_file_path = NULL,
  base = "MSstatsLiP_log_"
)
}
\arguments{
\item{data}{list of summarized datasets. Can be output of MSstatsLiP
summarization function \code{\link[MSstatsLiP]{dataSummarizationLiP}}. Must
include dataset named "LiP" as minimum.}

\item{contrast.matrix}{comparison between conditions of interests. Default
models full pairwise comparison between all conditions}

\item{fasta.path}{a file path to a fasta file that includes the proteins
listed in the data. Default is NULL. Include this parameter to determine
trypticity of peptides in LiP models.}
}
\value{
list of modeling results. Includes LiP, PROTEIN, and ADJUSTED LiP
data.tables with their corresponding model results.
}
\description{
Takes summarized LiP peptide and TrP protein data from dataSummarizationLiP
If global protein data is unavailable, LiP data only can be passed into the
function. Including protein data allows for adjusting LiP Fold Change by the
change in global protein abundance..
}
\examples{
# Convert and summarize data
fasta_path <- "../inst/extdata/ExampleFastaFile.fasta"

# Convert into MSstatsLiP format
MSstatsLiP_data <- SpectronauttoMSstatsLiPFormat(LiPRawData,
                                                 fasta_path,
                                                 TrPRawData)
# Run summarization without LiP missing value imputation
QuantData <- dataSummarizationLiP(MSstatsLiP_data)

# Test for pairwise comparison
ModelResults <- groupComparisonLiP(QuantData, contrast.matrix = "pairwise",
                                   fasta.path = fasta_path)

# Returns list of three models
names(ModelResults)
head(MSstatsLiP_model$LiP.Model)
head(MSstatsLiP_model$TrP.Model)
head(MSstatsLiP_model$Adjusted.LiP.Model)

}
