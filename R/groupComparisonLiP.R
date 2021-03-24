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
#' summarization function \code{\link[MSstatsLiP]{dataSummarizationLiP}}. Must
#' include dataset named "LiP" as minimum.
#' @param contrast.matrix comparison between conditions of interests. Default
#' models full pairwise comparison between all conditions
#' @param fasta.path a file path to a fasta file that includes the proteins
#' listed in the data
#' @return list of modeling results. Includes LiP, PROTEIN, and ADJUSTED LiP
#'         data.tables with their corresponding model results.
#' @examples
#' Add example
groupComparisonLiP <- function(data, contrast.matrix = "pairwise",
                               fasta.path = NA){

  ## TODO: Logging
  ## Put into format for MSstatsPTM function
  .groupComparisonCheck(data)
  data.LiP <- data[["LiP"]]
  Lip.processed <- as.data.table(data.LiP[["ProcessedData"]])
  Lip.run <- as.data.table(data.LiP[["RunlevelData"]])

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
  Trp.model <- model.data[['PROTEIN.Model']]
  Trp.model <- as.data.table(Trp.model)
  Adjusted.model <- model.data[['ADJUSTED.Model']]
  Adjusted.model <- as.data.table(Adjusted.model)

  LiP.model$FULL_PEPTIDE <- LiP.model$Protein
  LiP.model[,Protein:=NULL]

  Adjusted.model$FULL_PEPTIDE <- Adjusted.model$Protein
  Adjusted.model[,Protein:=NULL]

  if (!is.na(fasta.path)){
    ## Extracted protein name from LiP data
    ## Find unique proteins and peptide combinations
    available_proteins <- unique(as.character(Trp.model$Protein))
    available_proteins <- available_proteins[order(nchar(available_proteins),
                                                   available_proteins,
                                                   decreasing = TRUE)]
    available_ptms <- unique(as.character(Adjusted.model$FULL_PEPTIDE))

    ## Call Rcpp function to extract protein name
    ptm_proteins <- extract_protein_name(available_ptms,
                                         available_proteins)
    global_protein_lookup <- data.table(FULL_PEPTIDE = available_ptms,
                                        ProteinName = ptm_proteins)

    ## Add extracted protein name into dataset
    Adjusted.model <- merge(Adjusted.model, global_protein_lookup,
                            all.x = TRUE, by = 'FULL_PEPTIDE')

    setnames(Adjusted.model, old = "FULL_PEPTIDE", new = "PeptideSequence")
    Adjusted.model$PeptideSequence <- mapply(function(x,y){gsub(paste0(x, "_"),
                                                                "", y)},
                                             Adjusted.model$ProteinName,
                                             Adjusted.model$PeptideSequence)

    ## Load fasta
    format_fasta <- tidyFasta(fasta.path)
    format_fasta <- as.data.table(format_fasta)

    ## Get tryptic labels
    tryptic.label <- calculateTrypticity(Adjusted.model, format_fasta)

    ## Merge back into model
    Adjusted.model <- merge(Adjusted.model, tryptic.label, all.x = TRUE,
                            by = c("ProteinName", "PeptideSequence"))
    Adjusted.model$PeptideSequence <- paste(Adjusted.model$ProteinName,
                                            Adjusted.model$PeptideSequence,
                                            sep = "_")
    setnames(Adjusted.model, old = "PeptideSequence", new = "FULL_PEPTIDE")
    Adjusted.model[,ProteinName:=NULL]
    Adjusted.model[,GlobalProtein:=NULL]
  }

  ## Return models
  MSstats.Models <- list(
    LiP.Model = LiP.model,
    TrP.Model = Trp.model,
    Adjusted.LiP.Model = Adjusted.model
  )
  return(MSstats.Models)
}
