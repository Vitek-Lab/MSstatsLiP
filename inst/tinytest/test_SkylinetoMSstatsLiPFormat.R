
data("SkylineTest", package = "MSstatsLiP")

## SkylinetoMSstatsLiPFormat
## Test normal run
expect_silent(SkylinetoMSstatsLiPFormat(SkylineTest,
                                        TrP.data = SkylineTest,
                                        use_log_file = FALSE,
                                        append = FALSE))

## No Trp Data
expect_silent(SkylinetoMSstatsLiPFormat(SkylineTest,
                                        TrP.data = NULL,
                                        use_log_file = FALSE,
                                        append = FALSE))
## Test parameters
expect_error(SkylinetoMSstatsLiPFormat(LiPRawData,
                                       annotation = "test",
                                       use_log_file = FALSE,
                                       append = FALSE))

expect_error(SkylinetoMSstatsLiPFormat(LiPRawData,
                                       removeiRT = "test",
                                       use_log_file = FALSE,
                                       append = FALSE))

expect_error(SkylinetoMSstatsLiPFormat(LiPRawData,
                                       filter_with_Qvalue = "test",
                                       use_log_file = FALSE,
                                       append = FALSE))

expect_error(SkylinetoMSstatsLiPFormat(LiPRawData,
                                       qvalue_cutoff = "test",
                                       use_log_file = FALSE,
                                       append = FALSE))

expect_error(SkylinetoMSstatsLiPFormat(LiPRawData,
                                       useUniquePeptide = "test",
                                       use_log_file = FALSE,
                                       append = FALSE))

expect_error(SkylinetoMSstatsLiPFormat(LiPRawData,
                                       removeFewMeasurements = "test",
                                       use_log_file = FALSE,
                                       append = FALSE))

expect_error(SkylinetoMSstatsLiPFormat(LiPRawData,
                                       removeOxidationMpeptides = "test",
                                       use_log_file = FALSE,
                                       append = FALSE))

expect_error(SkylinetoMSstatsLiPFormat(LiPRawData,
                                       removeProtein_with1Feature = "test",
                                       use_log_file = FALSE,
                                       append = FALSE))
