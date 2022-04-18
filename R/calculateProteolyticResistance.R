#' Calcutates proteolytic resistance for provided data. Requires input from
#' dataSummarizationLiP function. Can optionally calculate differential
#' analysis using proteolytic resistance. In order for this function to work,
#' Conditions and run numbers must match between the LiP and TrP
#' data.
#'
#' @export
#' @importFrom MSstats groupComparison MSstatsContrastMatrix
#'
#' @param LiP_data name of variable containing LiP data. Must be output of
#' dataSummarizationLiP function.
#' @param fasta_file name of variable containing FASTA data. If FASTA file has
#' not been processed please run the tidyFasta() function on it before inputting
#' into this function. Protein names in file must match those in LiP_data.
#' @param differential_analysis logical indicating whether to run differential
#' analysis. Default is FALSE. Conditions and run numbers must match between
#' the LiP and TrP data.
#' @param contrast.matrix either a string of "pairwise" or a matrix including
#' what comparisons to make in the differential analysis. Only required if
#' differential_analysis=TRUE. Default is "pairwise".
#' @return a `data.frame` including either the summarized Proteolytic data or
#' differential analysis depending on parameter selection.
#' @examples
#' fasta <- tidyFasta(system.file("extdata", "ExampleFastaFile.fasta", package="MSstatsLiP"))
#' #calculateProteolyticResistance(MSstatsLiP_data, fasta)
calculateProteolyticResistance = function(LiP_data,
                                          fasta_file,
                                          differential_analysis = FALSE,
                                          contrast.matrix = "pairwise"){


  # TODO: Add checks
  # .checkProteolyticParams(LiP_data, fasta, differential_analysis)

  fully_TRI = ProteinName = uniprot_iso = Label = NULL
  PeptideSequence = sig = Coverage = Index = NULL

  ## Calculate tryptic
  calc_data = LiP_data$LiP$ProteinLevelData[, c("FULL_PEPTIDE", "Protein")]
  calc_data$PeptideSequence = unlist(
    lapply(calc_data$FULL_PEPTIDE, function(x){str_split(x, "_")[[1]][[2]]}))
  setnames(calc_data, c("Protein"), c("ProteinName"))
  tryptic_info = calculateTrypticity(calc_data, fasta_file)
  tryptic_info$FULL_PEPTIDE = paste(tryptic_info$ProteinName,
                                    tryptic_info$PeptideSequence, sep = "_")

  # Extract correct data.frames from input data
  lip_features = LiP_data$LiP$ProteinLevelData[, c("FULL_PEPTIDE", "Protein",
                                                   "LogIntensities", "GROUP",
                                                   "RUN")]
  lip_features = as.data.table(lip_features)
  trp_features = LiP_data$TrP$ProteinLevelData[, c("Protein", "LogIntensities",
                                                   "GROUP", "RUN")]
  trp_features = as.data.table(trp_features)

  # Load and format Fasta file
  ## Make sure file is loaded into memory
  if (identical(typeof(fasta_file), "character")){
    fasta_file = tidyFasta(fasta_file)
  }
  formated_fasta = as.data.table(fasta_file)
  formated_fasta = formated_fasta[, c("sequence", "uniprot_iso")]

  ## Test for missing proteins in FASTA file
  model.proteins = unique(c(lip_features$PROTEIN,
                            as.character(trp_features$PROTEIN)))
  fasta.proteins = unique(formated_fasta[, uniprot_iso])

  missing = setdiff(model.proteins, fasta.proteins)

  if (!identical(missing, character(0))){
    message(paste0("The following proteins are present in the data but ",
                   "not in the FASTA file. Note proteins must be in the FASTA ",
                   "in order to be plotted: ", paste(missing, collapse=", "))
    )
  }

  # Calculate Resistance for each peptide
  # Join lip and trp
  joined_data = merge(lip_features, trp_features,
                      by = c("Protein", "GROUP", "RUN"), all.x = TRUE)

  ## Bring tryptic info into data
  joined_data = merge(joined_data, tryptic_info, all.x = TRUE,
                      by = "FULL_PEPTIDE")
  joined_data = joined_data[fully_TRI == TRUE, ]

  joined_data$Accessibility_ratio = 2 ^ (joined_data$LogIntensities.x -
                                           joined_data$LogIntensities.y)
  joined_data$Accessibility_ratio = ifelse(joined_data$Accessibility_ratio > 1,
                                           1., joined_data$Accessibility_ratio)
  joined_data = joined_data[, c("FULL_PEPTIDE", "Protein", "PeptideSequence",
                                "GROUP", "RUN", "Accessibility_ratio")]

  MSstatsFormatDF = merge(joined_data, LiP_data$LiP$ProteinLevelData,
                          by = c("FULL_PEPTIDE", "Protein", "GROUP", "RUN"),
                          all.x = TRUE)
  MSstatsFormatDF[,LogIntensities:=NULL]

  if (differential_analysis){

    ## Prep data for group comp function
    da_input = LiP_data$LiP
    da_proteinleveldata = copy(MSstatsFormatDF)
    da_proteinleveldata$LogIntensities = da_proteinleveldata$Accessibility_ratio
    da_proteinleveldata$Protein = da_proteinleveldata$FULL_PEPTIDE
    da_input$ProteinLevelData = da_proteinleveldata

    ## Create pairwise matrix for label free
    if (contrast.matrix[1] == "pairwise"){
      labels = unique(da_input$ProteinLevelData$GROUP)
      contrast.matrix = MSstatsContrastMatrix('pairwise', labels)
    }

    ## TODO: replace group comparison with reference model
    da_analysis = MSstats::groupComparison(contrast.matrix, da_input)
    da_analysis$ComparisonResult = as.data.table(da_analysis$ComparisonResult)
    da_analysis$ComparisonResult = da_analysis$ComparisonResult[
      !is.na(da_analysis$ComparisonResult$Protein),]
    da_analysis$ComparisonResult = merge(da_analysis$ComparisonResult,
                                         unique(joined_data[,c("FULL_PEPTIDE",
                                                          "Protein",
                                                          "PeptideSequence")]),
                                         by.x = 'Protein',
                                         by.y = 'FULL_PEPTIDE', all.x = TRUE)
    da_analysis$ComparisonResult$FULL_PEPTIDE = da_analysis$ComparisonResult$Protein
    da_analysis$ComparisonResult$Protein = da_analysis$ComparisonResult$Protein.y
    da_analysis$ComparisonResult[,Protein.y:=NULL]

    return(list(RunLevelData = MSstatsFormatDF,
                groupComparison = da_analysis$ComparisonResult))
  } else{
    return(list(RunLevelData = MSstatsFormatDF))
  }

}
