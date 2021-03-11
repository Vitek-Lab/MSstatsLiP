#' Model LiP and TrP data and make adjustments if needed
#' Returns list of three modeled datasets
#'
#' Takes summarized LiP peptide and TrP protein data from dataSummarizationLiP
#' If global protein data is unavailable, LiP data only can be passed into the
#' function. Including protein data allows for adjusting LiP Fold Change by the
#' change in global protein abundance..
#'
#' @export
#' @importFrom MSstatsPTM groupComparisonPTM
#' @importFrom data.table as.data.table `:=`
#'
#' @param data list of summarized datasets. Can be output of MSstatsLiP
#' summarization function \code{\link[MSstatsLiP]{dataSummarizationLiP}}.
#' @param contrast.matrix comparison between conditions of interests. Default
#' models full pairwise comparison between all conditions
#'
#' @return list of modeling results. Includes LiP, PROTEIN, and ADJUSTED LiP
#'         data.tables with their corresponding model results.
#' @examples
#'
groupComparisonLiP <- function(data, contrast.matrix = "pairwise"){


  ## Put into format for MSstatsPTM function
  data.LiP <- data[["LiP"]]
  Lip.processed <- as.data.table(data.LiP[["ProcessedData"]])
  Lip.run <- as.data.table(data.LiP[["RunlevelData"]])

  Lip.processed$PROTEIN <- Lip.processed$FULL_PEPTIDE
  Lip.run$Protein <- Lip.run$FULL_PEPTIDE

  data.protein <- data[["TrP"]]
  format.data <- list(PTM = list(ProcessedData = Lip.processed,
                                 RunlevelData = Lip.run,
                                 SummaryMethod = data.LiP[["SummaryMethod"]],
                                 ModelQC = data.LiP[["ModelQC"]],
                                 PredictBySurvival =
                                   data.LiP[["PredictBySurvival"]]),
                      PROTEIN = data.protein)

  ## Model
  model.data <- groupComparisonPTM(format.data, "LabelFree", contrast.matrix)

  ## Format into LiP
  LiP.model <- model.data[['PTM.Model']]
  LiP.model <- as.data.table(LiP.model)
  Trp.model <- model.data[['PROTEIN.Model']]
  Trp.model <- as.data.table(Trp.model)
  Adjusted.model <- model.data[['ADJUSTED.Model']]
  Adjusted.model <- as.data.table(Adjusted.model)

  LiP.model$FULL_PEPTIDE <- LiP.model$Protein
  LiP.model[,Protein:=NULL]

  Adjusted.model$FULL_PEPTIDE <- Adjusted.model$Protein
  Adjusted.model[,Protein:=NULL]

  ## Return models
  MSstats.Models <- list(
    LiP.Model = LiP.model,
    TrP.Model = Trp.model,
    Adjusted.LiP.Model = Adjusted.model
  )

}
