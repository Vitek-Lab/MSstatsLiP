#' Converts raw LiP MS data from Skyline into the format needed for
#' MSstatsLiP.
#'
#' Takes as as input both raw LiP and Trp outputs from Skyline.
#'
#' @export
#' @importFrom data.table as.data.table `:=` .I
#'
#' @param LiP.data name of LiP Spectronaut output, which is long-format.
#' @param TrP.data name of TrP Spectronaut output, which is long-format.
#' @return a `list` of two data.frames in `MSstatsLiP` format
#' @examples
#'
#' ## Output will be in format
#' head(MSstatsLiP_data[["LiP"]])
#' head(MSstatsLiP_data[["TrP"]])

SkylinetoMSstatsLiPFormat <- function(LiP.data,
                                      TrP.data = NULL){

  LiP.data <- as.data.table(LiP.data)

  LiP.data$FULL_PEPTIDE <- paste(LiP.data$ProteinName,
                                 LiP.data$PeptideSequence, sep = "_")
  setnames(LiP.data, 'Replicate.Name', 'Run')

  LiP.data <- LiP.data[LiP.data[, .I[which.max(Intensity)],
                                by = c("PrecursorCharge", "FragmentIon",
                                       "ProductCharge", "IsotopeLabelType",
                                       "Condition", "BioReplicate", "Run",
                                       "FULL_PEPTIDE")]$V1]

  if (!is.null(TrP.data)){
    TrP.data <- as.data.table(TrP.data)
    setnames(TrP.data, 'Replicate.Name', 'Run')

    TrP.data <- TrP.data[TrP.data[, .I[which.max(Intensity)],
                                  by = c("ProteinName", "PeptideSequence",
                                         "PrecursorCharge", "FragmentIon",
                                         "ProductCharge", "IsotopeLabelType",
                                         "Condition", "BioReplicate",
                                         "Run")]$V1]
  }

  return(list(LiP = LiP.data, TrP = TrP.data))

}
