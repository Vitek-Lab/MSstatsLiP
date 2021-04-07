#' Barcode plot. Shows protein coverge of LiP modified peptides.
#'
#' @export
#' @import ggplot2
#' @importFrom data.table as.data.table `:=` rbindlist
#' @importFrom stringr str_match
#' @importFrom grDevices dev.off hcl pdf
#'
#' @param data list of data.tables containing LiP and TrP data in MSstatsLiP
#' format. Should be output of modeling function such as
#' \code{\link[MSstatsLiP]{groupComparisonLiP}}.
#' @param fasta A string of path to a FASTA file
#' @param model_type A string of either "Adjusted" or "Unadjusted", indicating
#' whether to plot the adjusted or unadjusted models. Default is "Adjusted".
#' @param which.prot a list of peptides to be visualized. Default is "all" which
#' will plot a separate barcode plot for each protein.
#' @param which.comp a list of comparisons to be visualized. Default is "all"
#' which will plot a separate barcode plot for each condition and protein.
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
#' # Convert and summarize data
#' #fasta_path <- "../inst/extdata/ExampleFastaFile.fasta"
#'
#' # Convert into MSstatsLiP format
#' #MSstatsLiP_data <- SpectronauttoMSstatsLiPFormat(LiPRawData,
#' #                                                  fasta_path,
#' #                                                  TrPRawData)
#' # Run summarization without LiP missing value imputation
#' #QuantData <- dataSummarizationLiP(MSstatsLiP_data)
#'
#' # Test for pairwise comparison
#' #ModelResults <- groupComparisonLiP(QuantData, contrast.matrix = "pairwise",
#' #                                    fasta_path)
#'
#' # Barcode Plot
#' #BarcodePlotLiP(MSstatsLiP_model, fasta_path,
#' #              model_type = "Adjusted")
#'
BarcodePlotLiP <- function(data,
                           fasta,
                           model_type = "Adjusted",
                           which.prot = "all",
                           which.comp = "all",
                           width = 12,
                           height = 4,
                           address = ""){

  ## TODO: Add checks on input and parameters
  ## TODO: Add logging
  if (toupper(model_type) == "ADJUSTED"){
    model.data <- data[["Adjusted.LiP.Model"]]
    model.data <- as.data.table(model.data)
  } else if (toupper(model_type) == "UNADJUSTED") {
    model.data <- data[["LiP.Model"]]
    model.data <- as.data.table(model.data)
  }

  formated_fasta <- tidyFasta(fasta)
  formated_fasta <- as.data.table(formated_fasta)
  formated_fasta <- formated_fasta[, c("sequence", "uniprot_iso")]

  ## Test for missing LiP proteins in FASTA file
  model.proteins <- unique(model.data[, ProteinName])
  fasta.proteins <- unique(formated_fasta[, uniprot_iso])

  missing <- setdiff(model.proteins, fasta.proteins)

  if (!identical(missing, character(0))){
    message(paste0("The following proteins are present in the LiP data but ",
                   "not in the FASTA file. Note proteins must be in the FASTA ",
                   "in order to be plotted: ", paste(missing, collapse=", "))
            )
  }

  ## Bring sequence into LiP data
  sig.coverage <- model.data[adj.pvalue < .05, c("ProteinName",
                                                 "PeptideSequence", "Label")]
  sig.coverage$sig <- TRUE
  insig.coverage <- model.data[adj.pvalue >= .05, c("ProteinName",
                                                    "PeptideSequence", "Label")]
  insig.coverage$sig <- FALSE

  coverage.df <- rbindlist(list(sig.coverage, insig.coverage))
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

  ## Determine which proteins to plot
  if (which.comp == "all"){
    which.comp <- unique(coverage.df[, Label])
  }

  for (c in seq(length(which.comp))){

    cond.coverage.df <- coverage.df[Label == which.comp[[c]], ]

    if (which.prot == "all"){
      which.prot <- unique(cond.coverage.df[, ProteinName])
    }

    for (i in seq(length(which.prot))){

      ## Build coverage data.table
      temp.seq <- formated_fasta[uniprot_iso == which.prot[[i]], sequence]
      coverage.index <- data.table("Index" = 1:nchar(temp.seq),
                                   "Coverage" = "No Coverage")

      temp.coverage.df <- cond.coverage.df[ProteinName == which.prot[[i]], ]

      for (idx in seq(nrow(temp.coverage.df))){
        cov_idx <- gregexpr(pattern = temp.coverage.df[idx, PeptideSequence],
                           temp.seq)
        start <- cov_idx[[1]][1] - 1
        end <- start + attr(cov_idx[[1]],'match.length')

        for (j in start:end){
          if (temp.coverage.df[idx, sig] == TRUE)
            coverage.index[j, Coverage := 'Significant']
          else if (temp.coverage.df[idx, sig] == FALSE){
            coverage.index[j, Coverage := 'Insignificant']
          }
        }
      }

      barcode_plot <- ggplot(data = coverage.index) +
        geom_col(aes(x = Index, y = 10, fill = Coverage)) +
        scale_fill_manual(values = c('Significant' = '#FF0000',
                                     'Insignificant' = '#808080',
                                     'No Coverage' = '#000000')) +
        labs(title = paste0(which.prot[[i]], " Coverage"), x = "Sequence",
             y = "") +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.background = element_rect(fill = 'white', colour = 'white')) +
        scale_x_continuous(breaks = 1:nchar(temp.seq),
                           labels = strsplit(temp.seq, split = "")[[1]])
      print(barcode_plot)
    }
  }

  if (address != FALSE) {
    dev.off()
  }

}
