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
#' @param log_base base of the logarithm used in dataProcess.
#' @param use_log_file logical. If TRUE, information about data processing
#' will be saved to a file.
#' @param append logical. If TRUE, information about data processing will be
#' added to an existing log file.
#' @param verbose logical. If TRUE, information about data processing will be
#' printed to the console.
#' @param log_file_path character. Path to a file to which information about
#' data processing will be saved.
#' If not provided, such a file will be created automatically.
#' If `append = TRUE`, has to be a valid path to a file.
#' @param base start of the file name.
#' @return list of modeling results. Includes LiP, PROTEIN, and ADJUSTED LiP
#'         data.tables with their corresponding model results.
#' @examples
#'
#' ## Use output of dataSummarizationLiP function
#' fasta <- system.file("extdata", "ExampleFastaFile.fasta", package="MSstatsLiP")
#'
#' # Test for pairwise comparison
#' MSstatsLiP_model <- groupComparisonLiP(MSstatsLiP_Summarized,
#'                                    contrast.matrix = "pairwise",
#'                                    fasta.path = fasta)
#'
#' # Returns list of three models
#' names(MSstatsLiP_model)
#' head(MSstatsLiP_model$LiP.Model)
#' head(MSstatsLiP_model$TrP.Model)
#' head(MSstatsLiP_model$Adjusted.LiP.Model)
#'
groupComparisonLiP <- function(data, contrast.matrix = "pairwise",
                               fasta.path = NULL,
                               log_base = 2,
                               use_log_file = TRUE,
                               append = FALSE,
                               verbose = TRUE,
                               log_file_path = NULL,
                               base = "MSstatsLiP_log_"){

  ## Start log
  if (is.null(log_file_path) & use_log_file == TRUE){
    time_now <- Sys.time()
    path <- paste0(base, gsub("[ :\\-]", "_", time_now),
                   ".log")
    file.create(path)
  } else {path <- log_file_path}

  ## Put into format for MSstatsPTM function
  .groupComparisonCheck(data)
  data.LiP <- data[["LiP"]]
  Lip.processed <- as.data.table(data.LiP[["FeatureLevelData"]])
  Lip.run <- as.data.table(data.LiP[["ProteinLevelData"]])

  ## keep peptide, protein match info
  lookup_table <- unique(Lip.processed[, c("PROTEIN", "FULL_PEPTIDE",
                                           "PEPTIDE")])
  setDT(lookup_table)[, "PEPTIDE" := tstrsplit(PEPTIDE, "_", keep = 1)]
  lookup_table <- unique(lookup_table)

  Lip.processed$PROTEIN <- as.factor(Lip.processed$FULL_PEPTIDE)
  Lip.run$Protein <- as.factor(Lip.run$FULL_PEPTIDE)

  data.protein <- data[["TrP"]]
  format.data <- list(PTM = list(
    FeatureLevelData = as.data.frame(Lip.processed),
     ProteinLevelData = as.data.frame(Lip.run),
     SummaryMethod = data.LiP[["SummaryMethod"]],
     ModelQC = data.LiP[["ModelQC"]],
     PredictBySurvival = data.LiP[["PredictBySurvival"]]),
                      PROTEIN = data.protein)

  ## Model
  model.data <- groupComparisonPTM(format.data, "LabelFree", contrast.matrix,
                                   FALSE, "BH", log_base, use_log_file, append,
                                   verbose, path, base)
  model.data$ADJUSTED.Model <- model.data$ADJUSTED.Model[!is.na(
    model.data$ADJUSTED.Model$Protein)]

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

    ## Add issue into adjusted model
    Adjusted.model <- merge(Adjusted.model, unique(LiP.model[,c("FULL_PEPTIDE",
                                                                "Label",
                                                                "issue")]),
                            by = c("FULL_PEPTIDE", "Label"), all.x = TRUE)

    ## Return models
    MSstats.Models <- list(
      LiP.Model = LiP.model,
      TrP.Model = Trp.model,
      Adjusted.LiP.Model = Adjusted.model
    )
  } else {
    ## Return models
    MSstats.Models <- list(
      LiP.Model = LiP.model,
      TrP.Model = NULL,
      Adjusted.LiP.Model = NULL
    )
  }

  return(MSstats.Models)
}
