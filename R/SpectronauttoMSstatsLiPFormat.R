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
#' @importFrom checkmate assertChoice assertLogical assertNumeric
#'
#' @param LiP.data name of LiP Spectronaut output, which is long-format.
#' @param fasta A string of path to a FASTA file, used to match LiP peptides.
#' @param Trp.data name of TrP Spectronaut output, which is long-format.
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
#' @param removeFewMeasurements TRUE (default) will remove the features that
#' have 1 or 2 measurements across runs.
#' @param removeProtein_with1Feature TRUE will remove the proteins which have
#' only 1 feature, which is the combination of peptide, precursor charge,
#' fragment and charge. FALSE is default.
#' @param removeNonUniqueProteins TRUE will remove proteins that were not
#' uniquely identified. IE if the protein column contains multiple proteins
#' seperated by ";". TRUE is default
#' @param removeModifications TRUE will remove peptide that contain a
#' modification. Modification must be indicated by "\[". TRUE is default
#' @param removeiRT TRUE will remove proteins that contain iRT. True is default
#' @param summaryforMultipleRows max(default) or sum - when there are multiple
#' measurements for certain feature and certain run, use highest or sum of
#' multiple intensities.
#' @param which.Conditions list of conditions to format into MSstatsLiP format.
#' If "all" all conditions will be used. Default is "all".
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
#' @return a `list` of two `data.frames` in MSstatsLiP format
#' @examples
#' # Output datasets of Spectronaut
#' head(LiPRawData)
#' head(TrPRawData)
#'
#' fasta_path <- system.file("extdata", "ExampleFastaFile.fasta", package="MSstatsLiP")
#'
#' MSstatsLiP_data <- SpectronauttoMSstatsLiPFormat(LiPRawData,
#'                                                  fasta_path,
#'                                                  TrPRawData)
#' head(MSstatsLiP_data[["LiP"]])
#' head(MSstatsLiP_data[["TrP"]])
SpectronauttoMSstatsLiPFormat <- function(LiP.data,
                                          fasta,
                                          Trp.data = NULL,
                                          annotation = NULL,
                                          intensity = 'PeakArea',
                                          filter_with_Qvalue = TRUE,
                                          qvalue_cutoff = 0.01,
                                          useUniquePeptide = TRUE,
                                          removeFewMeasurements = TRUE,
                                          removeProtein_with1Feature = FALSE,
                                          removeNonUniqueProteins = TRUE,
                                          removeModifications = TRUE,
                                          removeiRT = TRUE,
                                          summaryforMultipleRows=max,
                                          which.Conditions = 'all',
                                          use_log_file = TRUE,
                                          append = TRUE,
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

  MSstatsLogsSettings(use_log_file, append,
                      verbose, log_file_path = path)

  getOption("MSstatsLog")("INFO", "Starting parameter and data checks..")

  ## Check variable input
  .checkSpectronautConverterParams(intensity, filter_with_Qvalue,
                                   qvalue_cutoff, useUniquePeptide,
                                   removeProtein_with1Feature,
                                   summaryforMultipleRows)

  ## Ensure format of input data
  LiP.data <- as.data.table(LiP.data)
  if (!is.null(Trp.data)){
    Trp.data <- as.data.table(Trp.data)
  }

  ## Check and filter for available conditions
  if (which.Conditions != 'all') {

    LiP_conditions <- unique(LiP.data$R.Condition)
    if (!is.null(Trp.data)){
      TrP_conditions <- unique(Trp.data$R.Condition)
    } else {
      TrP_conditions <- LiP_conditions
    }
    if((length(setdiff(LiP_conditions, which.Conditions)
               ) == length(LiP_conditions)) |
       (length(setdiff(TrP_conditions, which.Conditions)
               ) == length(TrP_conditions))){
         msg = (
           paste("None of the conditions specified in which.Conditions are",
                  "available in one or both of LiP/TrP datasets. Please",
                  "ensure the conditions listed in which.Conditions appear",
                      "in the input datasets"))
         getOption("MSstatsLog")("ERROR", msg)
         stop(msg)
       }

    LiP.data <- LiP.data[(R.Condition %in% which.Conditions)]
    if (!is.null(Trp.data)){
      Trp.data <- Trp.data[(R.Condition %in% which.Conditions)]
    }

  }

  getOption("MSstatsLog")("INFO", "Formatting LiP data..")
  ## MSstats process
  df.lip <- SpectronauttoMSstatsFormat(LiP.data, annotation, intensity,
                                       filter_with_Qvalue, qvalue_cutoff,
                                       useUniquePeptide, removeFewMeasurements,
                                       removeProtein_with1Feature,
                                       summaryforMultipleRows, use_log_file,
                                       append, verbose, log_file_path = path,
                                       base)
  df.lip <- as.data.table(as.matrix(df.lip))
  if (!is.null(Trp.data)){

    getOption("MSstatsLog")("INFO", "Formatting TrP data..")
    df.trp <- SpectronauttoMSstatsFormat(as.data.frame(Trp.data), annotation,
                                         intensity,
                                         filter_with_Qvalue, qvalue_cutoff,
                                         useUniquePeptide,
                                         removeFewMeasurements,
                                         removeProtein_with1Feature,
                                         summaryforMultipleRows, use_log_file,
                                         append, verbose, log_file_path = path,
                                         base)
    df.trp <- as.data.table(as.matrix(df.trp))
  }

  getOption("MSstatsLog")("INFO", "Cleaning formatted LiP data..")

  ## Remove non-unique proteins and modified peptides if requested
  if (removeNonUniqueProteins){
    df.lip <- df.lip[!grepl(";", df.lip$ProteinName),]
  }
  if (removeModifications){
    df.lip <- df.lip[!grepl("\\[", df.lip$PeptideSequence),]
  }
  if (removeiRT){
    df.lip <- df.lip[!grepl("iRT", df.lip$ProteinName),]
  }

  df.lip$PeptideSequence <- str_extract(df.lip$PeptideSequence,
                                        "([ACDEFGHIKLMNPQRSTVWY]+)")
  df.lip$Intensity <- ifelse(df.lip$Intensity <= 1, NA, df.lip$Intensity)

  getOption("MSstatsLog")("INFO", "Joining with FASTA file..")

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
  MSstats_LiP <- merge(df.lip, df.fasta.lip[, c("ProteinName",
                                                "PeptideSequence")],
                       by = c("ProteinName", "PeptideSequence"))
  MSstats_LiP$FULL_PEPTIDE <- paste(MSstats_LiP$ProteinName,
                                    MSstats_LiP$PeptideSequence, sep = '_')

  if (!is.null(Trp.data)) {
    df.trp <- df.trp[which(!grepl(";", df.trp$ProteinName) &
                             !grepl("\\[", df.trp$PeptideSequence) &
                             !grepl("iRT", df.trp$ProteinName)),]
    df.trp$PeptideSequence <- str_extract(df.trp$PeptideSequence,
                                          "([ACDEFGHIKLMNPQRSTVWY]+)")
    df.trp$Intensity <- ifelse(df.trp$Intensity <= 1, NA, df.trp$Intensity)
    MSstats_TrP <- merge(df.trp, unique(df.fasta.lip[, "ProteinName"]),
                         by = "ProteinName")
  } else {
    MSstats_TrP = NULL
  }

  getOption("MSstatsLog")("INFO", "Converter finished, returning data.")

  LipExpt <- list(
    LiP = MSstats_LiP,
    TrP = MSstats_TrP
  )

  return(LipExpt)

}
