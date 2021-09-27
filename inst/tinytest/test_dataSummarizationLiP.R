data("MSstatsLiP_formatted_data", package = "MSstatsLiP")
data("summarized_data", package = "MSstatsLiP")

## Test missing columns are handled
expect_error(dataSummarizationLiP())
expect_error(dataSummarizationLiP(list(LiP = MSstatsLiP_data$LiP[, -1],
                                       LiP = MSstatsLiP_data$TrP)))
expect_error(dataSummarizationLiP(list(LiP = MSstatsLiP_data$LiP,
                                       TrP = MSstatsLiP_data$TrP[, -1])))
expect_error(dataSummarizationLiP(list(LiP = MSstatsLiP_data$LiP[, -2],
                                       TrP = MSstatsLiP_data$TrP)))
expect_error(dataSummarizationLiP(list(LiP = MSstatsLiP_data$LiP,
                                       TrP = MSstatsLiP_data$TrP[, -2])))
expect_error(dataSummarizationLiP(list(LiP = MSstatsLiP_data$LiP[, -3],
                                       TrP = MSstatsLiP_data$TrP)))
expect_error(dataSummarizationLiP(list(LiP = MSstatsLiP_data$LiP,
                                       TrP = MSstatsLiP_data$TrP[, -3])))
expect_error(dataSummarizationLiP(list(LiP = MSstatsLiP_data$LiP[, -4],
                                       TrP = MSstatsLiP_data$TrP)))
expect_error(dataSummarizationLiP(list(LiP = MSstatsLiP_data$LiP,
                                       TrP = MSstatsLiP_data$TrP[, -4])))
expect_error(dataSummarizationLiP(list(LiP = MSstatsLiP_data$LiP[, -5],
                                       TrP = MSstatsLiP_data$TrP)))
expect_error(dataSummarizationLiP(list(LiP = MSstatsLiP_data$LiP,
                                       TrP = MSstatsLiP_data$TrP[, -5])))
expect_error(dataSummarizationLiP(list(LiP = MSstatsLiP_data$LiP[, -6],
                                       TrP = MSstatsLiP_data$TrP)))
expect_error(dataSummarizationLiP(list(LiP = MSstatsLiP_data$LiP,
                                       TrP = MSstatsLiP_data$TrP[, -6])))
expect_error(dataSummarizationLiP(list(LiP = MSstatsLiP_data$LiP[, -7],
                                       TrP = MSstatsLiP_data$TrP)))
expect_error(dataSummarizationLiP(list(LiP = MSstatsLiP_data$LiP,
                                       TrP = MSstatsLiP_data$TrP[, -7])))
expect_error(dataSummarizationLiP(list(LiP = MSstatsLiP_data$LiP[, -8],
                                       TrP = MSstatsLiP_data$TrP)))
expect_error(dataSummarizationLiP(list(LiP = MSstatsLiP_data$LiP,
                                       TrP = MSstatsLiP_data$TrP[, -8])))
expect_error(dataSummarizationLiP(list(LiP = MSstatsLiP_data$LiP[, -9],
                                       TrP = MSstatsLiP_data$TrP)))
expect_error(dataSummarizationLiP(list(LiP = MSstatsLiP_data$LiP,
                                       TrP = MSstatsLiP_data$TrP[, -9])))

## Test misspecified parameters
expect_error(dataSummarizationLiP(MSstatsLiP_data, logTrans = 20))
expect_error(dataSummarizationLiP(MSstatsLiP_data, logTrans = "test"))
expect_error(dataSummarizationLiP(MSstatsLiP_data, normalization = TRUE))
expect_error(dataSummarizationLiP(MSstatsLiP_data, normalization.LiP = TRUE))
expect_error(dataSummarizationLiP(MSstatsLiP_data, fillIncompleteRows = 20))
expect_error(dataSummarizationLiP(MSstatsLiP_data, summaryMethod = TRUE))
expect_error(dataSummarizationLiP(MSstatsLiP_data, censoredInt = ""))
expect_error(dataSummarizationLiP(MSstatsLiP_data, cutoffCensored = ""))
expect_error(dataSummarizationLiP(MSstatsLiP_data, MBimpute = ""))
expect_error(dataSummarizationLiP(MSstatsLiP_data, MBimpute.LiP = ""))
