
data("summarized_data", package = "MSstatsLiP")
data("model_data", package = "MSstatsLiP")

## Test missing data is handled
expect_error(groupComparisonLiP())
expect_error(groupComparisonLiP(list(LiP = NULL,
                                     TrP = MSstatsLiP_Summarized$TrP)))

expect_silent(groupComparisonLiP(list(LiP = MSstatsLiP_Summarized$LiP,
                                      TrP = NULL)))

## Test expected output - LF
group.comp <- MSstatsLiP::groupComparisonLiP(MSstatsLiP_Summarized)
expect_inherits(group.comp, "list")
expect_inherits(group.comp$LiP.Model, "data.table")
expect_inherits(group.comp$TrP.Model, "data.table")
expect_inherits(group.comp$Adjusted.LiP.Model, "data.table")

expect_equal(colnames(group.comp$LiP.Model), c("FULL_PEPTIDE", "Label", "log2FC",
                                               "SE", "Tvalue", "DF", "pvalue",
                                               "adj.pvalue", "issue",
                                               "MissingPercentage",
                                               "ImputationPercentage", "ProteinName",
                                               "PeptideSequence"))
expect_equal(colnames(group.comp$TrP.Model), c("Protein", "Label", "log2FC",
                                               "SE", "Tvalue", "DF", "pvalue",
                                               "adj.pvalue", "issue",
                                               "MissingPercentage",
                                               "ImputationPercentage"))
expect_equal(colnames(group.comp$Adjusted.LiP.Model), c("FULL_PEPTIDE", "Label",
                                                    "log2FC", "SE", "Tvalue",
                                                    "DF", "pvalue",
                                                    "adj.pvalue",
                                                    "GlobalProtein",
                                                    "ProteinName",
                                                    "PeptideSequence", "issue"))
expect_equal(group.comp$Adjusted.LiP.Model, MSstatsLiP_model$Adjusted.LiP.Model)
