
data("LiPRawData", package = "MSstatsLiP")
data("TrPRawData", package = "MSstatsLiP")

## SpectronauttoMSstatsLiPFormat
## Test normal run
expect_silent(SpectronauttoMSstatsLiPFormat(LiPRawData,
                              "../extdata/ExampleFastaFile.fasta",
                              TrPRawData, use_log_file = FALSE,
                              append = FALSE))

## No Trp Data
expect_silent(SpectronauttoMSstatsLiPFormat(LiPRawData,
                                            "../extdata/ExampleFastaFile.fasta",
                                            use_log_file = FALSE,
                                            append = FALSE))

## Test parameters
expect_error(SpectronauttoMSstatsLiPFormat(LiPRawData,
                                            use_log_file = FALSE,
                                            append = FALSE))
expect_error(SpectronauttoMSstatsLiPFormat(LiPRawData,
                                           "../extdata/ExampleFastaFile.fasta",
                                           TrPRawData,
                                           intensity = TRUE,
                                           use_log_file = FALSE,
                                           append = FALSE))
expect_error(SpectronauttoMSstatsLiPFormat(LiPRawData,
                                           "../extdata/ExampleFastaFile.fasta",
                                           TrPRawData,
                                           filter_with_Qvalue  = "test",
                                           use_log_file = FALSE,
                                           append = FALSE))
expect_error(SpectronauttoMSstatsLiPFormat(LiPRawData,
                                           "../extdata/ExampleFastaFile.fasta",
                                           TrPRawData,
                                           useUniquePeptide  = "test",
                                           use_log_file = FALSE,
                                           append = FALSE))
expect_error(SpectronauttoMSstatsLiPFormat(LiPRawData,
                                           "../extdata/ExampleFastaFile.fasta",
                                           TrPRawData,
                                           removeFewMeasurements = "test",
                                           use_log_file = FALSE,
                                           append = FALSE))
expect_error(SpectronauttoMSstatsLiPFormat(LiPRawData,
                                           "../extdata/ExampleFastaFile.fasta",
                                           TrPRawData,
                                           which.Conditions = "test",
                                           use_log_file = FALSE,
                                           append = FALSE))
expect_error(SpectronauttoMSstatsLiPFormat(LiPRawData,
                                           "../extdata/ExampleFastaFile.fasta",
                                           TrPRawData,
                                           summaryforMultipleRows = "test",
                                           use_log_file = FALSE,
                                           append = FALSE))




