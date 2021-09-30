
data("summarized_data", package = "MSstatsLiP")

## Test dataProcessPlotsLiP
expect_error(dataProcessPlotsLiP())

## Incorrect parameters
expect_error(dataProcessPlotsLiP(MSstatsLiP_Summarized, type = "PROFILE"))
expect_error(dataProcessPlotsLiP(MSstatsLiP_Summarized, type = "QCPlo"))

expect_error(dataProcessPlotsLiP(MSstatsLiP_Summarized,
                                 which.Peptide = "test_broken"))
expect_error(dataProcessPlotsLiP(MSstatsLiP_Summarized,
                                 which.Protein = "test_broken"))
expect_error(dataProcessPlotsLiP(MSstatsLiP_Summarized, width = TRUE))
expect_error(dataProcessPlotsLiP(MSstatsLiP_Summarized, height = FALSE))
expect_error(dataProcessPlotsLiP(MSstatsLiP_Summarized, height = "test"))
expect_error(dataProcessPlotsLiP(MSstatsLiP_Summarized, ylimUp = "test"))
expect_error(dataProcessPlotsLiP(MSstatsLiP_Summarized, ylimDown = "test"))

## Normal plotting
expect_silent(dataProcessPlotsLiP(MSstatsLiP_Summarized, address = FALSE))
expect_silent(dataProcessPlotsLiP(MSstatsLiP_Summarized, type = "QCPLOT",
                                  address = FALSE))
