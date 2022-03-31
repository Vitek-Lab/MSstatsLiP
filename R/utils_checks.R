#' @noRd
.checkSpectronautConverterParams <- function(intensity,
                                             filter_with_Qvalue,
                                             qvalue_cutoff,
                                             useUniquePeptide,
                                             removeProtein_with1Feature,
                                             summaryforMultipleRows){

  ## Check params
  assertChoice(toupper(intensity), c('PEAKAREA', 'NORMALIZEDPEAKAREA'),
               .var.name = 'Intensity')

  assertLogical(filter_with_Qvalue, .var.name = "filter_with_Qvalue")

  assertNumeric(qvalue_cutoff, lower = 0, upper = 1,
                .var.name = "qvalue_cutoff")

  assertLogical(useUniquePeptide, .var.name = "useUniquePeptide")
  assertLogical(removeProtein_with1Feature,
                .var.name = "removeProtein_with1Feature")

  if (!(identical(summaryforMultipleRows, max)|
          identical(summaryforMultipleRows, sum))) {
    msg = (paste('summaryforMultipleRows must be one of max/sum. Please do not',
                 'input the values as strings.'))
    # getOption("MSstatsLog")("ERROR", msg)
    stop(msg)
  }
}

#' @noRd
.summarizeCheck <- function(data) {
  # Check the LiP data
  if (is.null(data[["LiP"]])) {
    msg <- paste('LiP peak list is missing. Input data must be of type list',
                  'with an element named \"LiP\" - stop')
    # getOption("MSstatsLog")("ERROR", msg)
    stop(msg)
  }
  LF.cols <- c("BioReplicate", "Condition", "FragmentIon", "Intensity",
                 "IsotopeLabelType", "PeptideSequence", "PrecursorCharge",
                 "ProductCharge", "FULL_PEPTIDE", "Run")
  provided_cols <- intersect(LF.cols, colnames(data[["LiP"]]))
  if (length(provided_cols) < 10){
    msg <- paste("Missing columns in the LiP input:",
                 paste(setdiff(LF.cols, colnames(data[["LiP"]])),
                       collapse = " "))
    #getOption("MSstatsLog")("ERROR", msg)
    stop(msg)
  }
  if (!is.null(data[["TrP"]])) {
    TrP.cols <- c("BioReplicate", "Condition", "FragmentIon", "Intensity",
                 "IsotopeLabelType", "PeptideSequence", "PrecursorCharge",
                 "ProductCharge", "ProteinName", "Run")
    provided_cols <- intersect(LF.cols, colnames(data[["TrP"]]))

    if (length(provided_cols) < 10){
      msg <- paste("Missing columns in the TrP input:",
                   paste(setdiff(LF.cols, colnames(data[["TrP"]])),
                         collapse = " "))
    }
  }
}

#' @noRd
.groupComparisonCheck <- function(data){
  # Check the LiP data
  if (is.null(data[["LiP"]])) {
    msg <- paste('LiP peak list is missing. Input data must be of type list',
                 'with an element named \"LiP\" - stop')
    # getOption("MSstatsLog")("ERROR", msg)
    stop(msg)
  }

  rawinput <- c("FULL_PEPTIDE", "PeptideSequence", "PrecursorCharge",
                "FragmentIon", "ProductCharge", "IsotopeLabelType",
                "Condition", "BioReplicate", "Run", "Intensity")
  if (length(setdiff(toupper(rawinput),
                     toupper(colnames(data[["LiP"]]$ProcessedData)))) == 0) {

    stop(paste("Please use 'dataProcessLiP' first. Then use output of",
               "dataProcessLiP function as input in groupComparison."))
  }
}

#' @noRd
.checkBarcodeParams <- function(data,
                                fasta,
                                model_type = "Adjusted",
                                which.prot = "all",
                                which.comp = "all",
                                width = 12,
                                height = 4,
                                address){

  if (model_type == "Adjusted") {
    if (is.null(data$Adjusted.LiP.Model)){
      stop(paste("If 'model_type' is Adjusted, the Adjusted.LiP.Model must be",
                  "included as a name in the data list."))
    } else {
      data <- data$Adjusted.LiP.Model
    }
  } else if (model_type == "Unadjusted"){
    if (is.null(data$LiP.Model)){
      stop(paste("If 'model_type' is Undjusted, the LiP.Model must be",
                 "included as a name in the data list."))
    } else {
      data <- data$LiP.Model
    }
  } else {
    stop("model_type must be one of: Adjusted, Unadjusted")
  }

  required_cols <- c("ProteinName", "PeptideSequence", "Label", "adj.pvalue")
  if (length(setdiff(required_cols, colnames(data))) > 0){
    stop(paste("The following columns are missing in the dataset:",
               paste(setdiff(required_cols, colnames(data)), collapse = ", ")))
  }

  if (which.comp[[1]] != "all") {
    if (length(setdiff(which.comp, unique(data$Label))) > 0) {

      stop(paste0("Please check comparison names ",
                  "Result does not have all comparisons listed - ",
                  toString(which.comp)))
    }
  }
  if (which.prot[[1]] != "all") {
    if (length(setdiff(which.prot, unique(data$ProteinName))) > 0) {

      stop(paste0("Please check labels of proteins ",
                  "Result does not have all these proteins. - ",
                  toString(which.prot)))
    }
  }

}

#' @noRd
.locateCheck <- function(peptide, uniprot, fasta, modResidue, modSymbol) {
  if (missing(peptide))
    stop("Input peptide is missing!")
  if (missing(uniprot))
    stop("Input uniprot is missing!")
  if (missing(fasta))
    stop("Input fasta is missing!")
  if (missing(modResidue))
    stop("Input modResidue is missing!")
  if (missing(modSymbol))
    stop("Input modSymbol is missing!")
  if (!is.character(peptide))
    stop("Please provide peptide sequence as character in peptide!")
  if (!is.character(uniprot))
    stop("Please provide Uniprot protein ID as character in uniprot!")
  if (length(peptide) != length(uniprot))
    stop("peptide and uniprot must be of the same length")
  if (!is.data.frame(fasta))
    stop("Please provide the FASTA information in a data frame!")
  if (!all(c("uniprot_iso", "sequence") %in% names(fasta)))
    stop("Uniprot_iso or sequence is missing from FASTA data frame!")
  TRUE
}
