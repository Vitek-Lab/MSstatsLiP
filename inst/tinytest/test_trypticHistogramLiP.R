
data("summarized_data", package = "MSstatsLiP")
data("model_data", package = "MSstatsLiP")
fasta_path = system.file("extdata", "ExampleFastaFile.fasta",
                         package="MSstatsLiP")

## trypticHistogramLiP testing
expect_error(trypticHistogramLiP())

## Parameter error testing
expect_error(trypticHistogramLiP(MSstatsLiP_Summarized, "blah"))
expect_error(trypticHistogramLiP(MSstatsLiP_Summarized, fasta_path,
                                 legened.size = FALSE))
expect_error(trypticHistogramLiP(MSstatsLiP_Summarized, fasta_path,
                                 color_scale = "purple"))

## Normal plotting
expect_silent(trypticHistogramLiP(MSstatsLiP_Summarized, fasta_path,
                                 address = FALSE))

expect_silent(trypticHistogramLiP(MSstatsLiP_Summarized, fasta_path,
                                  color_scale = "grey", address = FALSE))

expect_silent(trypticHistogramLiP(MSstatsLiP_Summarized, fasta_path,
                                  color_scale = "bright", address = FALSE))

## correlationPlotLiP
## correlationPlotLiP testing
expect_error(correlationPlotLiP())

## Parameter error testing
expect_error(correlationPlotLiP(MSstatsLiP_Summarized, method = FALSE))
expect_error(correlationPlotLiP(MSstatsLiP_Summarized, value_columns = FALSE))

## Normal plotting
expect_silent(correlationPlotLiP(MSstatsLiP_Summarized, address = FALSE))

## BarcodePlotLiP
## Test normal plot
expect_warning(StructuralBarcodePlotLiP(MSstatsLiP_model, fasta_path,
                                       address = FALSE))

## Test single protein
expect_warning(StructuralBarcodePlotLiP(MSstatsLiP_model, fasta_path,
                                       which.prot = "P36112"))

## Parameter checking
expect_error(StructuralBarcodePlotLiP(MSstatsLiP_model,
                                      fasta_path,
                            model_type = "test"))
expect_error(StructuralBarcodePlotLiP(MSstatsLiP_model,
                                      fasta_path,
                            which.prot = "test"))
expect_error(StructuralBarcodePlotLiP(MSstatsLiP_model,
                                      fasta_path,
                            which.comp = "test"))
expect_error(StructuralBarcodePlotLiP(MSstatsLiP_model,
                                      fasta_path,
                            FT.only = "test"))

## PCAPlotLiP
## Test normal plot
expect_silent(PCAPlotLiP(MSstatsLiP_Summarized, address = FALSE))

## Test individual plots
expect_silent(PCAPlotLiP(MSstatsLiP_Summarized,
                         which.comparison = c("Ctrl", "Osmo"),
                         address = FALSE))
expect_silent(PCAPlotLiP(MSstatsLiP_Summarized,
                         which.pep = c("P14164_ILQNDLK",
                                       "P17891_ALQLINQDDADIIGGRDR"),
                         address = FALSE))

## Parameter checking
expect_error(PCAPlotLiP(MSstatsLiP_Summarized,
                        bar.plot = "test",
                        address = FALSE))
expect_error(PCAPlotLiP(MSstatsLiP_Summarized,
                        protein.pca = "test",
                        address = FALSE))
expect_error(PCAPlotLiP(MSstatsLiP_Summarized,
                        comparison.pca = "test",
                        address = FALSE))
expect_error(PCAPlotLiP(MSstatsLiP_Summarized,
                        which.pep = "test",
                        address = FALSE))
expect_error(PCAPlotLiP(MSstatsLiP_Summarized,
                        which.comparison = "test",
                        address = FALSE))


