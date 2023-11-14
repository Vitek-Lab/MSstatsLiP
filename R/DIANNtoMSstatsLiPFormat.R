#' Converts raw LiP MS data from DIA-NN into the format needed for
#' MSstatsLiP.
#'
#' Takes as as input both raw LiP and Trp outputs from DIA-NN
#'
#' @export
#' @importFrom data.table as.data.table `:=` .I
#' @importFrom MSstats DIANNtoMSstatsFormat
#'
#' @param lip_data name of LiP Skyline output, which is long-format.
#' @param trp_data name of TrP Skyline output, which is long-format.
#' @param annotation name of 'annotation.txt' data which includes Condition,
#' BioReplicate, Run. If annotation is already complete in Skyline, use
#' annotation=NULL (default). It will use the annotation information from input.
#' @param global_qvalue_cutoff The global qvalue cutoff. Default is 0.01.
#' @param qvalue_cutoff Cutoff for DetectionQValue. Default is 0.01.
#' @param pg_qvalue_cutoff local qvalue cutoff for protein groups Run should be 
#' the same as filename. Default is .01.
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
DIANNtoMSstatsLiPFormat <- function(lip_data,
                                    trp_data = NULL,
                                    annotation = NULL,
                                    global_qvalue_cutoff=.01,
                                    qvalue_cutoff = 0.01,
                                    pg_qvalue_cutoff = 0.01,
                                    useUniquePeptide = TRUE,
                                    removeFewMeasurements = TRUE,
                                    removeOxidationMpeptides = FALSE,
                                    removeProtein_with1Feature = FALSE,
                                    use_log_file = FALSE,
                                    append = FALSE,
                                    verbose = TRUE,
                                    log_file_path = NULL){
  
  lip_data = as.data.table(lip_data)
  
  lip_msstats = DIANNtoMSstatsFormat(lip_data, annotation, global_qvalue_cutoff,
                                     qvalue_cutoff, pg_qvalue_cutoff, 
                                     useUniquePeptide,
                                     removeFewMeasurements, 
                                     removeOxidationMpeptides,
                                     removeProtein_with1Feature, use_log_file,
                                     append, verbose, log_file_path)
  
  lip_msstats$FULL_PEPTIDE = paste(lip_msstats$ProteinName,
                                   lip_msstats$PeptideSequence, sep = "_")
  
  
  if (!is.null(trp_data)){
    trp_data = as.data.table(trp_data)
    
    trp_msstats = DIANNtoMSstatsFormat(trp_data, annotation, 
                                         global_qvalue_cutoff,
                                         qvalue_cutoff, pg_qvalue_cutoff, 
                                         useUniquePeptide,
                                         removeFewMeasurements, 
                                         removeOxidationMpeptides,
                                         removeProtein_with1Feature, 
                                         use_log_file,
                                         append, verbose, log_file_path)
    
  } else {
    trp_msstats = NULL
  }
  
  return(list(LiP = lip_msstats, TrP = trp_msstats))
  
}
