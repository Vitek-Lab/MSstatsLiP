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
#' which will plot a separate barcode plot for each comparison and protein.
#' @param adj.pvalue.cutoff Default is .05. Alpha value for testing significance
#' of model output.
#' @param FC.cutoff Default is 0. Minimum absolute FC before a comparison will
#' be considered significant.
#' @param FT.only FALSE plots all FT and HT peptides, TRUE plots FT peptides
#' only. Default is FALSE.
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
#' # Specify Fasta path
#' fasta_path <- system.file("extdata", "ExampleFastaFile.fasta", package="MSstatsLiP")
#'
#' # Use model data to create Barcode Plot
#' StructuralBarcodePlotLiP(MSstatsLiP_model, fasta_path,
#'                          model_type = "Adjusted",
#'                          address=FALSE)
#'
StructuralBarcodePlotLiP <- function(data,
                           fasta,
                           model_type = "Adjusted",
                           which.prot = "all",
                           which.comp = "all",
                           adj.pvalue.cutoff = .05,
                           FC.cutoff = 0,
                           FT.only = FALSE,
                           width = 12,
                           height = 4,
                           address = ""){

  .checkBarcodeParams(data, fasta, model_type, which.prot, which.comp,
                      width, height, address)

  fully_TRI <- ProteinName <- uniprot_iso <- Label <- NULL
  PeptideSequence <- sig <- Coverage <- Index <- NULL

  if (toupper(model_type) == "ADJUSTED"){
    model.data <- data[["Adjusted.LiP.Model"]]
    model.data <- as.data.table(model.data)
  } else if (toupper(model_type) == "UNADJUSTED") {
    model.data <- data[["LiP.Model"]]
    model.data <- as.data.table(model.data)
  }

  if (FT.only){
    model.data <- model.data[fully_TRI == TRUE]
  } else {
    model.data <- model.data[fully_TRI == TRUE | NSEMI_TRI == TRUE |
                               CSEMI_TRI == TRUE]
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

  ## Calculate significance
  model.data$sig <-  (model.data$adj.pvalue < adj.pvalue.cutoff &
                        abs(model.data$log2FC) >= FC.cutoff)

  coverage.df <- model.data[, c("ProteinName", "PeptideSequence",
                                "Label", "sig")]

  ## Bring sequence into LiP data
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
      coverage.index <- data.table("Index" = seq_len(nchar(temp.seq)),
                                   "Coverage" = "No Coverage")

      temp.coverage.df <- cond.coverage.df[ProteinName == which.prot[[i]], ]

      for (idx in seq(nrow(temp.coverage.df))){
        cov_idx <- gregexpr(pattern = temp.coverage.df[idx, PeptideSequence],
                           temp.seq)
        start <- cov_idx[[1]][1] - 1
        end <- start + attr(cov_idx[[1]],'match.length')

        for (j in start:end){
          if (temp.coverage.df[idx, sig] == TRUE){
            coverage.index[j, Coverage := 'Significant']
          } else if (temp.coverage.df[idx, sig] == FALSE &
                   coverage.index[j, Coverage] != 'Significant'){
            coverage.index[j, Coverage := 'Insignificant']
          }
        }
      }

      barcode_plot <- ggplot(data = coverage.index) +
        geom_col(aes(x = Index, y = 10, fill = Coverage), width = 1) +
        scale_fill_manual(values = c('Significant' = '#FF0000',
                                     'Insignificant' = '#808080',
                                     'No Coverage' = '#000000')) +
        labs(title = paste0(which.prot[[i]], " Coverage - ", which.comp[[c]]),
             x = "Amino Acid Sequence", y = "") +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.background = element_rect(fill = 'white', colour = 'white'))
      print(barcode_plot)
    }
  }

  if (address != FALSE) {
    dev.off()
  }

}
