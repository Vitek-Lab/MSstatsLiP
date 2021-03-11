#' Converts raw LiP MS data from Spectronautt into the format needed for
#' MSstatsLiP.
#'
#' Takes as as input both raw LiP and Trp outputs from Spectronautt.
#'
#' @export
#' @importFrom MSstats SpectronauttoMSstatsFormat
#' @importFrom data.table as.data.table
#' @importFrom stringr str_extract str_count str_locate_all
#' @importFrom purrr map_int
#'
#' @param LiP.data name of LiP Spectronaut output, which is long-format.
#' @param Trp.data name of TrP Spectronaut output, which is long-format.
#' @param fasta A string of path to a FASTA file, used to match LiP peptides.
#' @param annotation name of 'annotation.txt' data which includes Condition,
#' BioReplicate, Run. If annotation is already complete in Spectronaut, use
#' annotation=NULL (default). It will use the annotation information from input.
#' @param intensity 'PeakArea'(default) uses not normalized peak area.
#' 'NormalizedPeakArea' uses peak area normalized by Spectronaut
#' @param filter_with_Qvalue TRUE(default) will filter out the intensities that
#' have greater than qvalue_cutoff in EG.Qvalue column. Those intensities will
#' be replaced with zero and will be considered as censored missing values for
#' imputation purpose.
#' @param qvalue_cutoff Cutoff for EG.Qvalue. default is 0.01.
#' @param useUniquePeptide TRUE(default) removes peptides that are assigned for
#' more than one proteins. We assume to use unique peptide for each protein.
#' @param fewMeasurements 'remove' (default) will remove the features that have
#' 1 or 2 measurements across runs.
#' @param removeProtein_with1Feature TRUE will remove the proteins which have
#' only 1 feature, which is the combination of peptide, precursor charge,
#' fragment and charge. FALSE is default.
#' @param summaryforMultipleRows max(default) or sum - when there are multiple
#' measurements for certain feature and certain run, use highest or sum of
#' multiple intensities.
#' @param which.Conditions list of conditions to format into MSstatsPTM format.
#' If "all" all conditions will be used. Default is "all".
#' @examples
#'
SpectronauttoMSstatsLiPFormat <- function(LiP.data,
                                          Trp.data,
                                          fasta,
                                          annotation = NULL,
                                          intensity = 'PeakArea',
                                          filter_with_Qvalue = TRUE,
                                          qvalue_cutoff = 0.01,
                                          useUniquePeptide = TRUE,
                                          fewMeasurements="remove",
                                          removeProtein_with1Feature = FALSE,
                                          summaryforMultipleRows=max,
                                          which.Conditions = 'all'){

  ## TODO: Add in function to check input
  #.check.input.data(LiP.data)
  #.check.input.data(Trp.data)
  ## TODO: Add logging
  ## TODO: Add variable checks

  LiP.data <- as.data.table(LiP.data)
  Trp.data <- as.data.table(Trp.data)

  if (which.Conditions != 'all') {

    LiP.data <- LiP.data[(R.Condition %in% which.Conditions)]
    Trp.data <- Trp.data[(R.Condition %in% which.Conditions)]

  }

  ## MSstats process
  df.lip <- SpectronauttoMSstatsFormat(as.data.frame(LiP.data), annotation, intensity,
                                       filter_with_Qvalue, qvalue_cutoff,
                                       useUniquePeptide, fewMeasurements,
                                       removeProtein_with1Feature,
                                       summaryforMultipleRows)
  df.lip <- as.data.table(df.lip) ## Temp
  df.trp <- SpectronauttoMSstatsFormat(as.data.frame(Trp.data), annotation, intensity,
                                       filter_with_Qvalue, qvalue_cutoff,
                                       useUniquePeptide, fewMeasurements,
                                       removeProtein_with1Feature,
                                       summaryforMultipleRows)
  df.trp <- as.data.table(df.trp) ## Temp

  ## Remove non-unique proteins and modified peptides
  df.trp <- df.trp[which(!grepl(";", df.trp$ProteinName) &
                         !grepl("\\[", df.trp$PeptideSequence) &
                         !grepl("iRT", df.trp$ProteinName)),]
  df.trp$PeptideSequence <- str_extract(df.trp$PeptideSequence,
                                        "([ACDEFGHIKLMNPQRSTVWY]+)")
  df.trp$Intensity <- ifelse(df.trp$Intensity <= 1, NA, df.trp$Intensity)

  df.lip <- df.lip[which(!grepl(";", df.lip$ProteinName) &
                           !grepl("\\[", df.lip$PeptideSequence) &
                           !grepl("iRT", df.lip$ProteinName)),]
  df.lip$PeptideSequence <- str_extract(df.lip$PeptideSequence,
                                        "([ACDEFGHIKLMNPQRSTVWY]+)")
  df.lip$Intensity <- ifelse(df.lip$Intensity <= 1, NA, df.lip$Intensity)

  ## Load and format FASTA file
  formated_fasta <- tidyFasta(fasta)
  formated_fasta <- as.data.table(formated_fasta)

  min_len_peptide <- 6
  df.fasta.lip <- merge(unique(df.lip[, c("ProteinName", "PeptideSequence")]),
        formated_fasta[, c("uniprot_iso", "sequence")], by.x = "ProteinName",
        by.y = "uniprot_iso")
  df.fasta.lip <- df.fasta.lip[
    which(nchar(df.fasta.lip$PeptideSequence) > min_len_peptide & str_count(
      df.fasta.lip$sequence, df.fasta.lip$PeptideSequence) == 1),]

  #Data formatting for MSstatsLiP analysis
  MSstats_LiP <- merge(df.lip, df.fasta.lip[, c("ProteinName", "PeptideSequence")],
                       by = c("ProteinName", "PeptideSequence"))
  MSstats_LiP$FULL_PEPTIDE <- paste(MSstats_LiP$ProteinName,
                                    MSstats_LiP$PeptideSequence, sep = '_')

  MSstats_LiP[,ProteinName:=NULL]
  MSstats_TrP <- merge(df.trp, unique(df.fasta.lip[, "ProteinName"]),
                       by = "ProteinName")

  LipExpt <- list(
    LiP = MSstats_LiP,
    TrP = MSstats_TrP
  )

  return(LipExpt)

}
