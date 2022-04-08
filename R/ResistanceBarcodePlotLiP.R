#' Proteolytic Resistance Barcode plot. Shows accessibility score of different
#' fully tryptic peptides in a protein.
#'
#' @export
#' @import ggplot2
#' @importFrom data.table as.data.table `:=` rbindlist merge
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
#' @param differential_analysis
#' @param which.comp a list of comparisons to be visualized, if differential
#' analysis is passed to plot_differential variable. Default is "all"
#' which will plot a separate barcode plot for each comparison and protein.
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
#' ResistanceBarcodePlotLiP(MSstatsLiP_model, fasta_path)
#'
ResistanceBarcodePlotLiP = function(data,
                           fasta,
                           which.prot = "all",
                           which.condition = "all",
                           differential_analysis = data.frame(),
                           which.comp = "all",
                           width = 12,
                           height = 4,
                           address = ""){
  # TODO: Figure out if exception should be calculated here or in tryptic function
  # TODO: Add checks for resistance barcode
  # .checkBarcodeParams(data, fasta, model_type, which.prot, which.comp,
  #                     width, height, address)

  fully_TRI = ProteinName = uniprot_iso = Label = NULL
  PeptideSequence = sig = Coverage = Index = NULL

  ## Calculate tryptic
  calc_data = data$LiP$ProteinLevelData[, c("FULL_PEPTIDE", "Protein")]
  calc_data$PeptideSequence = unlist(
    lapply(calc_data$FULL_PEPTIDE, function(x){str_split(x, "_")[[1]][[2]]}))
  setnames(calc_data, c("Protein"), c("ProteinName"))
  tryptic_info = calculateTrypticity(calc_data, fasta)
  tryptic_info$FULL_PEPTIDE = paste(tryptic_info$ProteinName,
                                    tryptic_info$PeptideSequence, sep = "_")

  # Extract correct data.frames from input data
  lip_features = data$LiP$ProteinLevelData[, c("FULL_PEPTIDE", "Protein",
                                               "LogIntensities", "GROUP", "SUBJECT")]
  lip_features = as.data.table(lip_features)
  lip_features = lip_features[, mean(LogIntensities),
                              by = list(FULL_PEPTIDE, Protein, GROUP)]
  trp_features = data$TrP$ProteinLevelData[, c("Protein", "LogIntensities",
                                               "GROUP", "SUBJECT")]
  trp_features = as.data.table(trp_features)
  trp_features = trp_features[, mean(LogIntensities), by = list(Protein, GROUP)]

  # Load and format Fasta file
  formated_fasta = tidyFasta(fasta)
  formated_fasta = as.data.table(formated_fasta)
  formated_fasta = formated_fasta[, c("sequence", "uniprot_iso")]

  ## Test for missing proteins in FASTA file
  model.proteins = unique(c(lip_features$PROTEIN,
                            as.character(trp_features$PROTEIN)))
  fasta.proteins = unique(formated_fasta[, uniprot_iso])

  missing = setdiff(model.proteins, fasta.proteins)

  if (!identical(missing, character(0))){
    message(paste0("The following proteins are present in the data but ",
                   "not in the FASTA file. Note proteins must be in the FASTA ",
                   "in order to be plotted: ", paste(missing, collapse=", "))
    )
  }

  # Calculate Resistance for each peptide
  # Join lip and trp
  joined_data = merge(lip_features, trp_features, by = c("Protein", "GROUP"),
                      all.x = TRUE)

  ## Bring tryptic info into data
  joined_data = merge(joined_data, tryptic_info, all.x = TRUE,
                      by = "FULL_PEPTIDE")
  joined_data = joined_data[fully_TRI == TRUE, ]

  joined_data$FULL_PEPTIDE = unlist(
    lapply(joined_data$FULL_PEPTIDE, function(x){str_split(x, "_")[[1]][[2]]}))
  joined_data$Accessibility_ratio = 2 ^ (joined_data$V1.x - joined_data$V1.y)
  joined_data$Accessibility_ratio = ifelse(joined_data$Accessibility_ratio > 1,
                                           1., joined_data$Accessibility_ratio)

  ## Bring sequence into data
  coverage.df = merge(joined_data, formated_fasta, all.x = TRUE,
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
    which.cond = unique(coverage.df[, GROUP])
  }

  for (c in seq(length(which.cond))){

    cond.coverage.df = coverage.df[GROUP == which.cond[[c]], ]

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
        cov_idx = gregexpr(pattern = temp.coverage.df[idx, FULL_PEPTIDE],
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
        labs(title = paste0(which.prot[[i]], " Coverage - ", which.cond[[c]]),
             x = "Amino Acid Sequence", y = "") +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.background = element_rect(fill = 'white', colour = 'white'))
      print(barcode_plot)
    }
  }

  # if (nrow(plot_differential) > 0){
  #   model_data = differential_analysis$Adjusted.LiP.Model
  #   if (which.comp == "all"){
  #     which.comp = unique(model_data[, Label])
  #   }
  #
  #   for (c in seq(length(which.comp))){
  #
  #     cond.coverage.df = coverage.df[GROUP == which.cond[[c]], ]
  #
  #     if (which.prot == "all"){
  #       which.prot = unique(cond.coverage.df[, Protein ])
  #     }
  #
  #     for (i in seq(length(which.prot))){
  #
  #       ## Build coverage data.table
  #       temp.seq = formated_fasta[uniprot_iso == which.prot[[i]], sequence]
  #       coverage.index = data.table("Index" = seq_len(nchar(temp.seq)),
  #                                   "Accessibility_ratio" = 0)
  #
  #       temp.coverage.df = cond.coverage.df[Protein == which.prot[[i]], ]
  #
  #       for (idx in seq(nrow(temp.coverage.df))){
  #         cov_idx = gregexpr(pattern = temp.coverage.df[idx, FULL_PEPTIDE],
  #                            temp.seq)
  #         start = cov_idx[[1]][1] - 1
  #         end = start + attr(cov_idx[[1]],'match.length')
  #
  #         for (j in start:end){
  #           coverage.index[j, Accessibility_ratio := temp.coverage.df[idx, Accessibility_ratio]]
  #         }
  #       }
  #       coverage.index$Accessibility_ratio = ifelse(coverage.index$Accessibility_ratio == 0.,
  #                                                   NA, coverage.index$Accessibility_ratio)
  #       barcode_plot = ggplot(data = coverage.index) +
  #         geom_col(aes(x = Index, y = 10, fill = Accessibility_ratio), width = 1) +
  #         scale_fill_gradient(low = "yellow", high = "red", limits = c(0,1),
  #                             name = "Proteolytic Resistance") +
  #         labs(title = paste0(which.prot[[i]], " Coverage - ", which.cond[[c]]),
  #              x = "Amino Acid Sequence", y = "") +
  #         theme(axis.text.y = element_blank(),
  #               axis.ticks.y = element_blank(),
  #               panel.background = element_rect(fill = 'white', colour = 'white'))
  #       print(barcode_plot)
  #     }
  #   }
  #
  # }

  if (address != FALSE) {
    dev.off()
  }

}
