#' Converts raw LiP MS data from Spectronautt into the format needed for
#' MSstatsLiP. Needs as input both LiP Trp outputs from Spectronautt.
#'
#' @export
#' @importFrom MSstats SpectronauttoMSstatsFormat
#' @importFrom data.table as.data.table
#' @importFrom stringr str_extract str_count str_locate_all
#' @importFrom purrr map_int
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
  ## TODO: y these matter tho?
  # df.fasta.lip$idx_peptide <- str_locate_all(df.fasta.lip$sequence,
  #                                            df.fasta.lip$PeptideSequence)
  # df.fasta.lip$aa_start <- map_int(df.fasta.lip$idx_peptide, ~.[, "start"])
  # df.fasta.lip$aa_end <- map_int(df.fasta.lip$idx_peptide, ~.[, "end"])

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
