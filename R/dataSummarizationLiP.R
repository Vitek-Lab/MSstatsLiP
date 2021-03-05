#' Summarizes LiP and TrP datasets seperately using methods from MSstats.
#' Returns list of two summarized datasets
#'
#' @export
#' @importFrom MSstatsPTM dataSummarizationPTM
#' @importFrom data.table as.data.table `=:`
dataSummarizationLiP <- function(
  data,
  logTrans = 2,
  normalization = "equalizeMedians",
  normalization.LiP = FALSE,
  nameStandards = NULL,
  nameStandards.LiP = NULL,
  address = "",
  fillIncompleteRows = TRUE,
  featureSubset = "all",
  featureSubset.LiP = "all",
  remove_uninformative_feature_outlier = FALSE,
  remove_uninformative_feature_outlier.LiP = FALSE, n_top_feature = 3,
  n_top_feature.LiP = 3,
  summaryMethod = "TMP", equalFeatureVar = TRUE, censoredInt = "NA",
  cutoffCensored = "minFeature", MBimpute = TRUE, MBimpute.LiP = FALSE,
  remove50missing = FALSE,
  fix_missing = NULL, maxQuantileforCensored = 0.999,clusters = NULL
) {

  ## TODO: Add logging
  ## TODO: Add parameter checking
  #.checkDataProcessParams()

  LiP.dataset <- data[["LiP"]]
  protein.dataset <- data[["TrP"]]

  # Check PTM and PROTEIN data for correct format
  #adj.protein <- .summarizeCheck(data, 'LF')

  LiP.dataset$ProteinName <- LiP.dataset$FULL_PEPTIDE

  format.data <- list(PTM = LiP.dataset, PROTEIN = protein.dataset)

  summarized.data <- dataSummarizationPTM(format.data, logTrans, normalization,
                       normalization.LiP, nameStandards, nameStandards.LiP,
                       address, fillIncompleteRows, featureSubset,
                       featureSubset.LiP, remove_uninformative_feature_outlier,
                       remove_uninformative_feature_outlier.LiP, n_top_feature,
                       n_top_feature.LiP, summaryMethod, equalFeatureVar,
                       censoredInt, cutoffCensored, MBimpute, MBimpute.LiP,
                       remove50missing, fix_missing, maxQuantileforCensored,
                       clusters)

  Lip.summarized <- summarized.data[["PTM"]]
  Lip.processed <- as.data.table(Lip.summarized[["ProcessedData"]])
  Lip.run <- as.data.table(Lip.summarized[["RunlevelData"]])

  ## Naming convention for LiP
  Lip.processed$FULL_PEPTIDE <- Lip.processed$PROTEIN
  Lip.processed[,PROTEIN:=NULL]
  Lip.run$FULL_PEPTIDE <- Lip.run$Protein
  Lip.run[,Protein:=NULL]
  Lip.summarized.format <- list(ProcessedData = Lip.processed,
                                RunlevelData = Lip.run,
                                SummaryMethod = Lip.summarized[["SummaryMethod"]],
                                ModelQC = Lip.summarized[["ModelQC"]],
                                PredictBySurvival = Lip.summarized[["PredictBySurvival"]])

  Trp.summarized <- summarized.data[["PROTEIN"]]

  MSstats.Summarized <- list(
    LiP = Lip.summarized.format,
    TrP = Trp.summarized
  )

  return(MSstats.Summarized)

}
