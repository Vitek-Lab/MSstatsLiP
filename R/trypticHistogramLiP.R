#' Histogram of Half vs Fully tryptic peptides. Calculates proteotypicity,
#' and then uses calcualtions in histogram.
#'
#' @export
#' @importFrom data.table as.data.table
trypticHistogramLiP <- function(data, fasta) {

  ## TODO: Add checks on input and parameters
  ## TODO: Add logging

  ## Format input data
  lip.data <- data[["LiP"]]
  trp.data <- data[["TrP"]]

  formated_fasta <- tidyFasta(fasta)
  formated_fasta <- as.data.table(formated_fasta)
  #formated_fasta <- formated_fasta[, c("sequence", "uniprot_iso")]

  ## TODO: Replace this with more robust Rcpp function
  regex_protein <- '([^-]+)(?:_[^-]+){1}$'
  lip.data[, ProteinName := factor(str_match(FULL_PEPTIDE, regex_protein)[,2])]

  sequences <- unique(lip.data[, c("PeptideSequence", "ProteinName")])

}
