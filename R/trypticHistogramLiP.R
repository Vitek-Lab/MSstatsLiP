#' Histogram of Half vs Fully tryptic peptides. Calculates proteotypicity,
#' and then uses calcualtions in histogram.
#'
#' @export
#' @import ggplot2
#' @importFrom data.table as.data.table `:=`
#' @importFrom grDevices dev.off pdf
#' @importFrom scales percent
#'
#' @param data output of MSstatsLiP converter function. Must include at least
#' ProteinName, PeptideSequence, BioReplicate, and Condition columns
#' @param fasta A string of path to a FASTA file, used to match LiP peptides.
#' @param x.axis.size size of x-axis labeling for plot. Default is 10.
#' @param y.axis.size size of y-axis labeling for plot. Default is 10.
#' @param legend.size size of feature legend for half vs fully tryptic peptides
#' below graph. Default is 7.
#' @param address the name of folder that will store the results. Default folder
#'  is the current working directory. The other assigned folder has to be
#'  existed under the current working directory. An output pdf file is
#'  automatically created with the default name of "TyrpticPlot.pdf". If
#'  address=FALSE, plot will be not saved as pdf file but shown in window..
#' @return plot or pdf
#' @examples
#' # Specify fasta file
#' fasta_path <- "../data/ExampleFastaFile.fasta"
#'
#' MSstatsLiP_data <- SpectronauttoMSstatsLiPFormat(LiPRawData,
#'                                                  fasta_path,
#'                                                  TrPRawData)
#' trypticHistogramLiP(MSstatsLiP_data, fasta_path, address = FALSE)
#'
trypticHistogramLiP <- function(data, fasta, x.axis.size = 10,
                                y.axis.size = 10, legend.size = 10,
                                address = "") {

  ## TODO: Add checks on input and parameters
  ## TODO: Add logging

  ## Format input data
  lip.data <- data[["LiP"]]
  trp.data <- data[["TrP"]]

  format_fasta <- tidyFasta(fasta)
  format_fasta <- as.data.table(format_fasta)

  ## Add tryptic data
  tryptic.label <- calculateTrypticity(lip.data, format_fasta)


  plot_df <- merge(lip.data, tryptic.label[, c("ProteinName", "PeptideSequence",
                                               "fully_TRI")], all.x = TRUE,
                   by = c("ProteinName", "PeptideSequence"))

  plot_df[, count := .N, by=.(Condition, BioReplicate, fully_TRI)]
  plot_df <- unique(plot_df[, c("Condition", "BioReplicate", "fully_TRI",
                                "count")])

  plot_df[, sum := sum(count), by = .(Condition, BioReplicate)]
  plot_df$percent_plot <- plot_df$count / plot_df$sum

  if (address != FALSE) {
    allfiles <- list.files()

    num <- 0
    filenaming <- paste0(address, "TyrpticPlot")
    finalfile <- paste0(address, "TyrpticPlot.pdf")

    while (is.element(finalfile, allfiles)) {
      num <- num + 1
      finalfile <- paste0(paste(filenaming, num, sep = "-"), ".pdf")
    }

    pdf(finalfile, width = width, height = height)
  }

  hist_temp <- ggplot(data = plot_df) +
    geom_col(aes(x = BioReplicate, y = percent_plot, fill = fully_TRI)) +
    facet_wrap(.~Condition, scales = "free") +
    labs(title = "Proteotrypticity", x = "Replicate", y = "Percent") +
    scale_fill_manual(values = c("red", "blue"), labels = c("Half", "Full"),
                      name = "Trypticity") +
    scale_y_continuous(labels = percent) +
    theme(
      panel.background = element_rect(fill = 'white', colour = "black"),
      legend.key = element_rect(fill = 'white', colour = 'white'),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = 'gray95'),
      axis.text.y = element_text(size = y.axis.size, colour = "black"),
      axis.ticks = element_line(colour = "black"),
      axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
      axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
      title = element_text(size = x.axis.size + 8, vjust = 1.5),
      legend.position = "bottom",
      legend.text = element_text(size = legend.size))

  print(hist_temp)

  if (address != FALSE) {
    dev.off()
  }

}
