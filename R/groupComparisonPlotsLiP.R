#' Visualization for model-based analysis and summarization
#'
#' @export
#'
groupComparisonPlotsPTM <- function(data = data,
                                    type = type,
                                    sig=0.05,
                                    FCcutoff=FALSE,
                                    logBase.pvalue=10,
                                    ylimUp=FALSE,
                                    ylimDown=FALSE,
                                    xlimUp=FALSE,
                                    x.axis.size=10,
                                    y.axis.size=10,
                                    dot.size=3,
                                    text.size=4,
                                    text.angle=0,
                                    legend.size=13,
                                    ProteinName=TRUE,
                                    colorkey=TRUE,
                                    numProtein=100,
                                    clustering="both",
                                    width=10,
                                    height=10,
                                    which.Comparison="all",
                                    which.Protein="all",
                                    which.Models="all",
                                    address="") {

  ## Format into PTM
  LiP.model <- model.data[['LiP.Model']]
  Trp.model <- model.data[['TrP.Model']]
  Adjusted.model <- model.data[['Adjusted.Model']]

  LiP.model$ProteinName <- LiP.model$FULL_PEPTIDE
  Adjusted.model$ProteinName <- Adjusted.model$FULL_PEPTIDE

  formated.data <- list(PTM.Model = LiP.model,
                        PROTEIN.Model = Trp.model,
                        ADJUSTED.Model = Adjusted.model)

  groupComparisonPlotsPTM(formated.data, type, sig, FCcutoff, logBase.pvalue,
                          ylimUp, ylimDown, xlimUp, x.axis.size, y.axis.size,
                          dot.size, text.size, text.angle, legend.size,
                          ProteinName, colorkey, numProtein, clustering,
                          width, height, which.Comparison, which.Protein,
                          which.Models, address)

}
