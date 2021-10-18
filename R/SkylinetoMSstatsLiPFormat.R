#' Converts raw LiP MS data from Skyline into the format needed for
#' MSstatsLiP.
#'
#' Takes as as input both raw LiP and Trp outputs from Skyline.
#'
#' @export
#' @importFrom data.table as.data.table `:=` .I
#' @importFrom MSstats SkylinetoMSstatsFormat
#'
#' @param LiP.data name of LiP Skyline output, which is long-format.
#' @param TrP.data name of TrP Skyline output, which is long-format.
#' @param annotation name of 'annotation.txt' data which includes Condition,
#' BioReplicate, Run. If annotation is already complete in Skyline, use
#' annotation=NULL (default). It will use the annotation information from input.
#' @param removeiRT TRUE (default) will remove the proteins or peptides which
#' are labeld 'iRT' in 'StandardType' column. FALSE will keep them.
#' @param filter_with_Qvalue TRUE(default) will filter out the intensities that
#' have greater than qvalue_cutoff in DetectionQValue column. Those intensities
#' will be replaced with zero and will be considered as censored missing values
#' for imputation purpose.
#' @param qvalue_cutoff Cutoff for DetectionQValue. default is 0.01.
#' @param useUniquePeptide TRUE (default) removes peptides that are assigned
#' for more than one proteins. We assume to use unique peptide for each protein.
#' @param removeFewMeasurements TRUE (default) will remove the features that
#' have 1 or 2 measurements across runs.
#' @param removeOxidationMpeptides TRUE will remove the peptides including
#' 'oxidation (M)' in modification. FALSE is default.
#' @param removeProtein_with1Feature TRUE will remove the proteins which have
#' only 1 feature, which is the combination of peptide, precursor charge,
#' fragment and charge. FALSE is default.
#' @param use_log_file logical. If TRUE, information about data processing will
#' be saved to a file.
#' @param append logical. If TRUE, information about data processing will be
#' saved to a file.
#' @param verbose logical. If TRUE, information about data processing wil be
#' printed to the console.
#' @param log_file_path character. Path to a file to which information about
#' data processing will be saved. If not provided, such a file will be created
#' automatically. If 'append = TRUE', has to be a valid path to a file.
#' @return a `list` of two data.frames in `MSstatsLiP` format
#' @examples
#'
#' ## Output will be in format
#' head(MSstatsLiP_data[["LiP"]])
#' head(MSstatsLiP_data[["TrP"]])

SkylinetoMSstatsLiPFormat <- function(LiP.data,
                                      TrP.data = NULL,
                                      annotation = NULL,
                                      removeiRT = TRUE,
                                      filter_with_Qvalue = TRUE,
                                      qvalue_cutoff = 0.01,
                                      useUniquePeptide = TRUE,
                                      removeFewMeasurements = TRUE,
                                      removeOxidationMpeptides = FALSE,
                                      removeProtein_with1Feature = FALSE,
                                      use_log_file = FALSE,
                                      append = FALSE,
                                      verbose = TRUE,
                                      log_file_path = NULL){

  LiP.data <- as.data.table(LiP.data)

  LiP.Skyline <- SkylinetoMSstatsFormat(LiP.data, annotation, removeiRT,
                                        filter_with_Qvalue,
                         qvalue_cutoff, useUniquePeptide,
                         removeFewMeasurements, removeOxidationMpeptides,
                         removeProtein_with1Feature, use_log_file,
                         append, verbose, log_file_path)

  LiP.Skyline$FULL_PEPTIDE <- paste(LiP.Skyline$ProteinName,
                                    LiP.Skyline$PeptideSequence, sep = "_")


  if (!is.null(TrP.data)){
    TrP.data <- as.data.table(TrP.data)

    TrP.Skyline <- SkylinetoMSstatsFormat(TrP.data, annotation, removeiRT,
                                          filter_with_Qvalue,
                                          qvalue_cutoff, useUniquePeptide,
                                          removeFewMeasurements,
                                          removeOxidationMpeptides,
                                          removeProtein_with1Feature,
                                          use_log_file,
                                          append, verbose, log_file_path)

  } else {
    TrP.Skyline <- NULL
  }

  return(list(LiP = LiP.Skyline, TrP = TrP.Skyline))

}
