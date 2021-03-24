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
    # processout <- rbind(processout, paste("The required input - LiP data : did",
    #                                       "not process from dataProcess",
    #                                       "function. - stop"))
    stop(paste("Please use 'dataProcessLiP' first. Then use output of",
               "dataProcessLiP function as input in groupComparison."))
  }
}
