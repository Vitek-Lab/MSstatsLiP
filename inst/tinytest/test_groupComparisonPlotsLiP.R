
data("model_data", package = "MSstatsLiP")

## Test groupComparisonPlotsLiP
expect_error(groupComparisonPlotsLiP())

## Incorrect parameters
expect_error(groupComparisonPlotsLiP(MSstatsLiP_model, type = "Volanoplt"))
expect_error(groupComparisonPlotsLiP(MSstatsLiP_model, type = "Heetmap"))

expect_error(groupComparisonPlotsLiP(MSstatsLiP_model, numProtein  = "yes"))
expect_error(groupComparisonPlotsLiP(MSstatsLiP_model, sig = "test_broken"))
expect_error(groupComparisonPlotsLiP(MSstatsLiP_model, FCcutoff = TRUE))
expect_error(groupComparisonPlotsLiP(MSstatsLiP_model, logBase.pvalue = TRUE))
expect_error(groupComparisonPlotsLiP(MSstatsLiP_model, dot.size = TRUE))

expect_error(groupComparisonPlotsLiP(MSstatsLiP_model, height = FALSE))
expect_error(groupComparisonPlotsLiP(MSstatsLiP_model, height = "test"))
expect_error(groupComparisonPlotsLiP(MSstatsLiP_model, ylimUp = "test"))
expect_error(groupComparisonPlotsLiP(MSstatsLiP_model, ylimDown = "test"))

## Normal plotting
expect_silent(groupComparisonPlotsLiP(MSstatsLiP_model, type = "VolcanoPlot",
                                      address = FALSE))
expect_silent(groupComparisonPlotsLiP(MSstatsLiP_model, type = "VolcanoPlot",
                                      numProtein = 5, address = FALSE))
expect_silent(groupComparisonPlotsLiP(MSstatsLiP_model, type = "Heatmap",
                                      address = FALSE))
expect_silent(groupComparisonPlotsLiP(MSstatsLiP_model, type = "Heatmap",
                                      numProtein = 5, address = FALSE))

