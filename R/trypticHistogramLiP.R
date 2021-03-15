#' Histogram of Half vs Fully tryptic peptides. Calculates proteotypicity,
#' and then uses calcualtions in histogram.
#'
#' @export
#' @import ggplot2
#' @importFrom data.table as.data.table `:=`
trypticHistogramLiP <- function(data, fasta, x.axis.size = 10,
                                y.axis.size = 10, legend.size = 10) {

  ## TODO: Add checks on input and parameters
  ## TODO: Add logging
  ## TODO: Add address and pdf saving functionality

  ## Format input data
  lip.data <- data[["LiP"]]
  trp.data <- data[["TrP"]]

  format_fasta <- tidyFasta(fasta)
  format_fasta <- as.data.table(format_fasta)
  #formated_fasta <- formated_fasta[, c("sequence", "uniprot_iso")]

  ## TODO: Replace this with more robust Rcpp function
  ## All proteins
  available_proteins <- unique(as.character(trp.data$ProteinName))
  available_proteins <- available_proteins[order(nchar(available_proteins),
                                                 available_proteins,
                                                 decreasing = TRUE)]
  available_ptms <- unique(as.character(lip.data$FULL_PEPTIDE))

  ## Call Rcpp function
  ptm_proteins <- extract_protein_name(available_ptms,
                                       available_proteins)
  global_protein_lookup <- data.table(FULL_PEPTIDE = available_ptms,
                                      ProteinName = ptm_proteins)

  ## Add extracted protein name into model
  lip.data <- merge(lip.data, global_protein_lookup,
                     all.x = TRUE, by = 'FULL_PEPTIDE')

  tryptic.label <- calculateTrypticity(lip.data, format_fasta)

  plot_df <- merge(lip.data, tryptic.label[, c("ProteinName", "PeptideSequence",
                                               "fully_TRI")], all.x = TRUE,
                   by = c("ProteinName", "PeptideSequence"))

  plot_df[, count := .N, by=.(Condition, BioReplicate, fully_TRI)]
  plot_df <- unique(plot_df[, c("Condition", "BioReplicate", "fully_TRI", "count")])

  ggplot(data = plot_df) +
    geom_col(aes(x = BioReplicate, y = count, fill = fully_TRI)) +
    facet_wrap(.~Condition) +
    labs(title = "Proteotrypticity", x = "Replicate", y = "count") +
    theme(
      panel.background = element_rect(fill = 'white', colour = "black"),
      legend.key = element_rect(fill = 'white', colour = 'white'),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = 'gray95'),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = y.axis.size, colour = "black"),
      axis.ticks = element_line(colour = "black"),
      axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
      axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
      title = element_text(size = x.axis.size + 8, vjust = 1.5),
      legend.position = "bottom",
      legend.text = element_text(size = legend.size))


}
