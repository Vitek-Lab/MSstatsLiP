#' Calculates level of trypticity
#' @noRd
calculateTrypticity <- function(LiP_data, fasta_file){

  unique_pep <- unique(LiP_data[, c("PeptideSequence",
                                    "ProteinName")])

  ## Add fasta file to get full sequence
  combined_df <- merge(unique_pep, fasta_file, all.x = TRUE,
                       by.x = "ProteinName", by.y = "uniprot_iso")

  ## Identify pre/end/post AA
  combined_df$start <- mapply(function(x,y){gregexpr(x, y)[[1]][1] - 1},
                              combined_df$PeptideSequence, combined_df$sequence)
  combined_df$end <- combined_df$start + nchar(combined_df$PeptideSequence)
  combined_df$preAA <- mapply(function(x,y) {substr(y, x, x)},
                              combined_df$start, combined_df$sequence)
  combined_df$endAA <- mapply(function(x,y) {substr(y, x, x)},
                              combined_df$end, combined_df$sequence)
  combined_df$postAA <- mapply(function(x,y) {substr(y, x, x)},
                               combined_df$end + 1, combined_df$sequence)

  ## Determine trip
  combined_df$fully_TRI <- ifelse(combined_df$preAA %in% c("K", "R") &
                                    combined_df$endAA %in% c("K", "R"),
                                  TRUE, FALSE)
  combined_df$NSEMI_TRI <- ifelse(combined_df$preAA %in% c("K", "R") &
                                    !combined_df$endAA %in% c("K", "R"),
                                  TRUE, FALSE)
  combined_df$CSEMI_TRI <- ifelse(!combined_df$preAA %in% c("K", "R") &
                                    combined_df$endAA %in% c("K", "R"),
                                  TRUE, FALSE)
  combined_df$CTERMINUS <- ifelse(combined_df$postAA == "", TRUE, FALSE)

  return(combined_df[, c("ProteinName", "PeptideSequence", "fully_TRI",
                         "NSEMI_TRI", "CSEMI_TRI", "CTERMINUS")])
}
