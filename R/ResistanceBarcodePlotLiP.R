#' Proteolytic Resistance Barcode plot. Shows accessibility score of different
#' fully tryptic peptides in a protein.
#'
#' @export
#' @import ggplot2
#' @importFrom data.table as.data.table `:=` rbindlist
#' @importFrom stringr str_match str_split
#' @importFrom grDevices dev.off hcl pdf
#'
#' @param data list of data.tables containing LiP and TrP data in MSstatsLiP
#' format. Should be output of summarization function as
#' \code{\link[MSstatsLiP]{dataSummarizationLiP}}.
#' @param fasta A string of path to a FASTA file
#' @param which.prot a list of peptides to be visualized. Default is "all" which
#' will plot a separate barcode plot for each protein.
#' @param which.condition  a list of conditions to be visualized. Default is
#' "all" which will plot all conditions for a single protein in the same barcode
#' plot.
#' @param differential_analysis a boolean indicating if a barcode plot showing 
#' the differential analysis should be plotted. If this is selected you must 
#' have performed differential analysis on the proteoltic data in the 
#' `calculateProteolyticResistance` function. Default is `FALSE`.
#' @param which.comp a list of comparisons to be visualized, if differential
#' analysis is passed to plot_differential variable. Default is "all"
#' which will plot a separate barcode plot for each comparison and protein.
#' @param adj.pvalue.cutoff Default is .05. Alpha value for testing significance
#' of model output.
#' @param FC.cutoff Default is 0. Minimum absolute FC before a comparison will
#' be considered significant.
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
#' fasta_path = system.file("extdata", "ExampleFastaFile.fasta", package="MSstatsLiP")
#'
#' # Use model data to create Barcode Plot
#' #ResistanceBarcodePlotLiP(MSstatsLiP_model, fasta_path)
#'
ResistanceBarcodePlotLiP = function(data,
                           fasta,
                           which.prot = "all",
                           which.condition = "all",
                           differential_analysis = FALSE,
                           which.comp = "all",
                           adj.pvalue.cutoff = .05,
                           FC.cutoff = 0,
                           width = 12,
                           height = 4,
                           address = ""){

  # Load and format Fasta file
  ## Make sure file is loaded into memory
  if (identical(typeof(fasta_file), "character")){
    fasta_file = tidyFasta(fasta_file)
  }
  formated_fasta = as.data.table(fasta_file)
  formated_fasta = formated_fasta[, c("sequence", "uniprot_iso")]

  ## Bring sequence into data
  coverage.df = merge(data$RunLevelData, formated_fasta, all.x = TRUE,
                      by.x = "Protein", by.y = "uniprot_iso")

  ## Create PDF to save plots if requested
  if (address != FALSE) {
    allfiles = list.files()

    num = 0
    filenaming = paste0(address, "ResistanceBarcode_Plot")
    finalfile = paste0(address, "Resistance_Barcode_Plot.pdf")

    while (is.element(finalfile, allfiles)) {
      num = num + 1
      finalfile = paste0(paste(filenaming, num, sep = "-"), ".pdf")
    }
    pdf(finalfile, width = width, height = height)
  }

  ## Determine which proteins to plot
  if (which.condition == "all"){
    which.condition = unique(coverage.df[, GROUP])
  }

  for (c in seq(length(which.condition))){

    cond.coverage.df = coverage.df[GROUP == which.condition[[c]], ]

    if (which.prot == "all"){
      which.prot = unique(cond.coverage.df[, Protein ])
    }

    for (i in seq(length(which.prot))){

      ## Build coverage data.table
      temp.seq = formated_fasta[uniprot_iso == which.prot[[i]], sequence]
      coverage.index = data.table("Index" = seq_len(nchar(temp.seq)),
                                   "Accessibility_ratio" = 0)

      temp.coverage.df = cond.coverage.df[Protein == which.prot[[i]], ]

      for (idx in seq(nrow(temp.coverage.df))){
        cov_idx = gregexpr(pattern = temp.coverage.df[idx, PeptideSequence],
                            temp.seq)
        start = cov_idx[[1]][1] - 1
        end = start + attr(cov_idx[[1]],'match.length')

        for (j in start:end){
          coverage.index[j, Accessibility_ratio := temp.coverage.df[idx, Accessibility_ratio]]
        }
      }
      coverage.index$Accessibility_ratio = ifelse(coverage.index$Accessibility_ratio == 0.,
                                                  NA, coverage.index$Accessibility_ratio)
      barcode_plot = ggplot(data = coverage.index) +
        geom_col(aes(x = Index, y = 10, fill = Accessibility_ratio), width = 1) +
        scale_fill_gradient(low = "yellow", high = "red", limits = c(0,1),
                            name = "Proteolytic Resistance") +
        labs(title = paste0(which.prot[[i]], " Coverage - ", which.condition[[c]]),
             x = "Amino Acid Sequence", y = "") +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.background = element_rect(fill = 'white', colour = 'white'))
      print(barcode_plot)
    }
  }

  #
  if (differential_analysis  == TRUE){
    model.data = data$groupComparison
    # model.data = model.data[is.finite(model.data$log2FC), ]
    model.data$sig = (model.data$adj.pvalue < adj.pvalue.cutoff &
                        abs(model.data$log2FC) >= FC.cutoff & 
                        is.finite(model.data$log2FC))

    coverage.df = model.data[, c("Protein", "PeptideSequence",
                                  "Label", "sig")]
    ## Bring sequence into data
    coverage.df = merge(coverage.df, formated_fasta, all.x = TRUE,
                        by.x = "Protein", by.y = "uniprot_iso")
    if (which.comp == "all"){
      which.comp = unique(coverage.df[, Label])
    }

    for (c in seq(length(which.comp))){

      cond.coverage.df = coverage.df[Label == which.comp[[c]], ]

      if (which.prot == "all"){
        which.prot = unique(cond.coverage.df[, Protein ])
      }

      for (i in seq(length(which.prot))){

        ## Build coverage data.table
        temp.seq = formated_fasta[uniprot_iso == which.prot[[i]], sequence]
        coverage.index = data.table("Index" = seq_len(nchar(temp.seq)),
                                    "Coverage" = "No Coverage")

        temp.coverage.df = cond.coverage.df[Protein == which.prot[[i]], ]

        for (idx in seq(nrow(temp.coverage.df))){
          cov_idx = gregexpr(pattern = temp.coverage.df[idx, PeptideSequence],
                             temp.seq)
          start = cov_idx[[1]][1] - 1
          end = start + attr(cov_idx[[1]],'match.length')

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

  }

  if (address != FALSE) {
    dev.off()
  }

}
