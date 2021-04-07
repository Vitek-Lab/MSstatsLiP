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
#' @importFrom data.table as.data.table `:=` tstrsplit setnames setDT
#'
#' @param data list of summarized datasets. Can be output of MSstatsLiP
#' summarization function \code{\link[MSstatsLiP]{dataSummarizationLiP}}. Must
#' include dataset named "LiP" as minimum.
#' @param contrast.matrix comparison between conditions of interests. Default
#' models full pairwise comparison between all conditions
#' @param fasta.path a file path to a fasta file that includes the proteins
#' listed in the data. Default is NULL. Include this parameter to determine
#' trypticity of peptides in LiP models.
#' @return list of modeling results. Includes LiP, PROTEIN, and ADJUSTED LiP
#'         data.tables with their corresponding model results.
#' @examples
#' # Convert and summarize data
#' fasta_path <- "../inst/extdata/ExampleFastaFile.fasta"
#'
#' # Convert into MSstatsLiP format
#' MSstatsLiP_data <- SpectronauttoMSstatsLiPFormat(LiPRawData,
#'                                                  fasta_path,
#'                                                  TrPRawData)
#' # Run summarization without LiP missing value imputation
#' QuantData <- dataSummarizationLiP(MSstatsLiP_data)
#'
#' # Test for pairwise comparison
#' ModelResults <- groupComparisonLiP(QuantData, contrast.matrix = "pairwise",
#'                                    fasta.path = fasta_path)
#'
#' # Returns list of three models
#' names(ModelResults)
#' head(MSstatsLiP_model$LiP.Model)
#' head(MSstatsLiP_model$TrP.Model)
#' head(MSstatsLiP_model$Adjusted.LiP.Model)
#'
groupComparisonLiP <- function(data, contrast.matrix = "pairwise",
                               fasta.path = NULL){

  ## TODO: Logging
  ## Put into format for MSstatsPTM function
  .groupComparisonCheck(data)
  data.LiP <- data[["LiP"]]
  Lip.processed <- as.data.table(data.LiP[["ProcessedData"]])
  Lip.run <- as.data.table(data.LiP[["RunlevelData"]])

  ## keep peptide, protein match info
  lookup_table <- unique(Lip.processed[, c("PROTEIN", "FULL_PEPTIDE",
                                           "PEPTIDE")])
  setDT(lookup_table)[, "PEPTIDE" := tstrsplit(PEPTIDE, "_", keep = 1)]

  Lip.processed$PROTEIN <- Lip.processed$FULL_PEPTIDE
  Lip.run$Protein <- Lip.run$FULL_PEPTIDE

  data.protein <- data[["TrP"]]
  format.data <- list(PTM = list(ProcessedData = as.data.frame(Lip.processed), ## TODO: Replace when new MSstats version comes out
                                 RunlevelData = as.data.frame(Lip.run),
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

  LiP.model$FULL_PEPTIDE <- LiP.model$Protein
  LiP.model[,Protein:=NULL]

  LiP.model <- merge(LiP.model, lookup_table, all.x = TRUE,
                     by = "FULL_PEPTIDE")

  setnames(LiP.model, old = c("PROTEIN", "PEPTIDE"),
           new = c("ProteinName", "PeptideSequence"))

  if (!is.null(fasta.path)){

    ## Load fasta
    format_fasta <- tidyFasta(fasta.path)
    format_fasta <- as.data.table(format_fasta)

    ## Get tryptic labels
    tryptic.label <- calculateTrypticity(LiP.model, format_fasta)

    ## Merge back into model
    LiP.model <- merge(LiP.model, tryptic.label, all.x = TRUE,
                            by = c("ProteinName", "PeptideSequence"))
  }

  MSstats.Models <- list(LiP.Model = LiP.model)

  Trp.model <- model.data[['PROTEIN.Model']]
  Trp.model <- as.data.table(Trp.model)
  Adjusted.model <- model.data[['ADJUSTED.Model']]
  Adjusted.model <- as.data.table(Adjusted.model)

  if (nrow(Trp.model) != 0){

    Adjusted.model$FULL_PEPTIDE <- Adjusted.model$Protein
    Adjusted.model[,Protein:=NULL]

    ## Add protein name back in
    Adjusted.model <- merge(Adjusted.model, lookup_table, all.x = TRUE,
                            by = "FULL_PEPTIDE")

    setnames(Adjusted.model, old = c("PROTEIN", "PEPTIDE"),
             new = c("ProteinName", "PeptideSequence"))

    if (!is.null(fasta.path)){

      ## Get tryptic labels
      tryptic.label <- calculateTrypticity(Adjusted.model, format_fasta)

      ## Merge back into model
      Adjusted.model <- merge(Adjusted.model, tryptic.label, all.x = TRUE,
                              by = c("ProteinName", "PeptideSequence"))
      Adjusted.model[,GlobalProtein:=NULL]
    }

    ## Return models
    MSstats.Models <- list(
      LiP.Model = LiP.model,
      TrP.Model = Trp.model,
      Adjusted.LiP.Model = Adjusted.model
    )
  }
  return(MSstats.Models)
}
