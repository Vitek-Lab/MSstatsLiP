#' Visualization for explanatory data analysis
#'
#' To illustrate the quantitative data and quality control of MS runs,
#' dataProcessPlotsTMT takes the quantitative data from MSstatsTMT converter
#' functions as input
#' and generate two types of figures in pdf files as output :
#' (1) profile plot (specify "ProfilePlot" in option type), to identify the
#' potential sources of variation for each protein;
#' (2) quality control plot (specify "QCPlot" in option type), to evaluate the
#' systematic bias between MS runs.
#'
#' @export
#' @import ggplot2
#' @importFrom graphics axis image legend mtext par plot.new title plot
#' @importFrom grDevices dev.off hcl pdf
#' @importFrom dplyr mutate n
#' @importFrom reshape2 dcast
#' @importFrom MSstatsTMT proteinSummarization
#' @importFrom gridExtra grid.arrange
#' @param data.ptm name of the data with PTM sites in protein name, which can be
#'  the output of MSstatsTMT converter functions.
#' @param data.protein name of the data with peptide level, which can be the
#' output of MSstatsTMT converter functions.
#' @param data.ptm.summarization name of the data with ptm sites in
#' protein-level name, which can be the output of
#' the MSstatsTMT \code{\link[MSstatsTMT]{proteinSummarization}} function.
#' @param data.protein.summarization name of the data with protein-level, which
#' can be the output of the
#' MSstatsTMT \code{\link[MSstatsTMT]{proteinSummarization}} function.
#' @param type choice of visualization. "ProfilePlot" represents profile plot of
#'  log intensities across MS runs.
#' "QCPlot" represents box plots of log intensities across channels and MS runs.
#' @param ylimUp upper limit for y-axis in the log scale.
#' FALSE(Default) for Profile Plot and QC Plot uses the upper limit as rounded
#' off maximum of log2(intensities) after normalization + 3..
#' @param ylimDown lower limit for y-axis in the log scale. FALSE(Default) for
#' Profile Plot and QC Plot uses 0..
#' @param x.axis.size size of x-axis labeling for "Run" and "channel in Profile
#' Plot and QC Plot.
#' @param y.axis.size size of y-axis labels. Default is 10.
#' @param text.size size of labels represented each condition at the top of
#' Profile plot and QC plot. Default is 4.
#' @param text.angle angle of labels represented each condition at the top of
#' Profile plot and QC plot. Default is 0.
#' @param legend.size size of legend above Profile plot. Default is 7.
#' @param dot.size.profile size of dots in Profile plot. Default is 2.
#' @param ncol.guide number of columns for legends at the top of plot. Default
#' is 5.
#' @param width width of the saved pdf file. Default is 10.
#' @param height height of the saved pdf file. Default is 10.
#' @param which.Protein Protein list to draw plots. List can be names of
#' Proteins or order numbers of Proteins.
#' Default is "all", which generates all plots for each protein. For QC plot,
#' "allonly" will generate one QC plot with all proteins.
#' @param originalPlot TRUE(default) draws original profile plots, without
#' normalization.
#' @param summaryPlot TRUE(default) draws profile plots with protein
#' summarization for each channel and MS run.
#' @param address the name of folder that will store the results. Default folder
#'  is the current working directory.
#' The other assigned folder has to be existed under the current working
#' directory.
#' An output pdf file is automatically created with the default name of
#' "ProfilePlot.pdf" or "QCplot.pdf".
#' The command address can help to specify where to store the file as well as
#' how to modify the beginning of the file name.
#' If address=FALSE, plot will be not saved as pdf file but showed in window.
#' @return plot or pdf
#' @examples
dataProcessPlotsLiP <- function(data,
                                type = 'PROFILEPLOT',
                                ylimUp = FALSE,
                                ylimDown = FALSE,
                                x.axis.size = 10,
                                y.axis.size = 10,
                                text.size = 4,
                                text.angle = 90,
                                legend.size = 7,
                                dot.size.profile = 2,
                                ncol.guide = 5,
                                width = 10,
                                height = 12,
                                lip.title = "All Peptides",
                                protein.title = "All Proteins",
                                which.Protein = "all",
                                originalPlot = TRUE,
                                summaryPlot = TRUE,
                                address = "") {

  ## Format into PTM format
  Lip.data <- data[["LiP"]]
  Lip.data.Processed <- Lip.data$ProcessedData
  Lip.data.Processed$PROTEIN <- Lip.data.Processed$FULL_PEPTIDE
  Lip.data.Processed <- Lip.data.Processed %>% select(-FULL_PEPTIDE)
  Lip.data.Run <- Lip.data$RunlevelData
  Lip.data.Run$Protein <- Lip.data.Run$FULL_PEPTIDE
  Lip.data.Run <- Lip.data.Run %>% select(-FULL_PEPTIDE)

  TrP.data <- data[["TrP"]]

  format.data <- list(PTM = list(ProcessedData = Lip.data.Processed,
                                 RunlevelData = Lip.data.Run,
                                 SummaryMethod = Lip.data[["SummaryMethod"]],
                                 ModelQC = Lip.data[["ModelQC"]],
                                 PredictBySurvival = Lip.data[["PredictBySurvival"]]),
                      PROTEIN = TrP.data)

  ## Plot
  dataProcessPlotsPTM(format.data, type, ylimUp, ylimDown, x.axis.size,
                      y.axis.size, text.size, text.angle, legend.size,
                      dot.size.profile, ncol.guide, width, height, lip.title,
                      protein.title, which.Protein, originalPlot, summaryPlot,
                      address)

}
