#' Barcode plot. Shows protein coverge of LiP modified peptides.
#'
#' @export
#' @import ggplot2
#' @importFrom data.table as.data.table `:=`
#' @importFrom stringr str_match
#' @importFrom grDevices dev.off hcl pdf
#'
#' @param data list of data.tables containing LiP and TrP data in MSstatsLiP
#' format. Can be the output of converters such as
#' \code{\link[MSstatsLiP]{SpectronauttoMSstatsLiPFormat}}.
#' @param fasta A string of path to a FASTA file
#' @param which.prot a list of peptides to be visualized.
#' @param width width of the saved file. Default is 10.
#' @param height height of the saved file. Default is 10.
#' @param address the name of folder that will store the results. Default
#' folder is the current working directory. The other assigned folder has to
#' be existed under the current working directory. An output pdf file is
#' automatically created with the default name of "VolcanoPlot.pdf" or
#' "Heatmap.pdf". The command address can help to specify where to store the
#' file as well as how to modify the beginning of the file name. If
#' address=FALSE, plot will be not saved as pdf file but showed in window
#' @return plot or pdf
#' @examples
#' Add example
BarcodePlotLiP <- function(data,
                           fasta,
                           which.prot = "all",
                           width = 12,
                           height = 4,
                           address = ""){

  ## TODO: Add checks on input and parameters
  ## TODO: Add logging
  lip.data <- data[["LiP"]]
  trp.data <- data[["TrP"]]

  formated_fasta <- tidyFasta(fasta)
  formated_fasta <- as.data.table(formated_fasta)
  formated_fasta <- formated_fasta[, c("sequence", "uniprot_iso")]

  ## TODO: Replace this with more robust Rcpp function
  regex_protein <- '([^-]+)(?:_[^-]+){1}$'
  lip.data[, ProteinName := factor(str_match(FULL_PEPTIDE, regex_protein)[,2])]

  ## Test for missing LiP proteins in FASTA file
  lip.proteins <- unique(lip.data[, ProteinName])
  trp.proteins <- unique(trp.data[, ProteinName])
  fasta.proteins <- unique(formated_fasta[, uniprot_iso])

  missing.lip <- setdiff(lip.proteins, fasta.proteins)
  missing.trp <- setdiff(trp.proteins, fasta.proteins)

  if (!identical(missing.lip, character(0))){
    message(paste0("The following proteins are present in the LiP data but ",
                   "not in the FASTA file. Note proteins must be in the FASTA ",
                   "in order to be plotted: ", paste(missing.lip,collapse=", "))
            )
  }

  ## Bring sequence into LiP data
  lip.coverage <- unique(lip.data[, c("PeptideSequence", "ProteinName")])
  lip.coverage$type <- "lip"

  trp.coverage <- unique(trp.data[, c("PeptideSequence", "ProteinName")])
  trp.coverage$type <- "trp"

  coverage.df <- rbindlist(list(lip.coverage, trp.coverage))
  coverage.df <- merge(coverage.df, formated_fasta, all.x = TRUE,
                       by.x = "ProteinName", by.y = "uniprot_iso")

  ## Create PDF to save plots if requested
  if (address != FALSE) {
    allfiles <- list.files()

    num <- 0
    filenaming <- paste0(address, "Barcode_Plot")
    finalfile <- paste0(address, "Barcode_Plot.pdf")

    while (is.element(finalfile, allfiles)) {
      num <- num + 1
      finalfile <- paste0(paste(filenaming, num, sep = "-"), ".pdf")
    }
    pdf(finalfile, width = width, height = height)
  }

  for (protein in which.prot) {

    ## Build coverage data.table
    temp.seq <- formated_fasta[uniprot_iso == protein, sequence]
    coverage.index <- data.table("Index" = 1:nchar(temp.seq), 'Coverage' = "N")

    temp.coverage.df <- coverage.df[ProteinName == protein, ]

    for (idx in 1:nrow(temp.coverage.df)){
      cov_idx <- gregexpr(pattern = temp.coverage.df[idx, PeptideSequence],
                         temp.seq)
      start <- cov_idx[[1]][1] - 1
      end <- start + attr(cov_idx[[1]],'match.length')

      for (j in start:end){
        if (temp.coverage.df[idx, type] == 'lip')
          coverage.index[j, Coverage := 'lip']
        else if (coverage.index[j, Coverage] == "N"){
          coverage.index[j, Coverage := 'trp']
        }
      }
    }

    barcode_plot <- ggplot(data = coverage.index) +
      geom_col(aes(x = Index, y = 10, fill = Coverage)) +
      scale_fill_manual(values=c("#FF0000", "#808080", "#000000")) +
      labs(title = paste0(protein, " Coverage"), x = "Sequence", y = "") +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.background = element_rect(fill = 'white', colour = 'white')) +
      scale_x_continuous(breaks = 1:nchar(temp.seq),
                         labels = strsplit(temp.seq, split = "")[[1]])
    print(barcode_plot)
  }

  if (address != FALSE) {
    dev.off()
  }

}
