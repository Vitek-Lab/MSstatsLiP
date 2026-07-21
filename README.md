# MSstatsLiP

<!-- badges: start -->
[![Bioconductor Release Build](https://bioconductor.org/shields/build/release/bioc/MSstatsLiP.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/MSstatsLiP/)
[![Codecov test coverage](https://codecov.io/gh/Vitek-Lab/MSstatsLiP/branch/devel/graph/badge.svg)](https://codecov.io/gh/Vitek-Lab/MSstatsLiP/branch/devel)
[![Bioconductor Downloads Rank](https://bioconductor.org/shields/downloads/release/MSstatsLiP.svg)](https://bioconductor.org/packages/stats/bioc/MSstatsLiP/)
[![Years in Bioconductor](https://bioconductor.org/shields/years-in-bioc/MSstatsLiP.svg)](https://bioconductor.org/packages/release/bioc/html/MSstatsLiP.html#since)
[![License: Artistic-2.0](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](https://opensource.org/licenses/Artistic-2.0)
<!-- badges: end -->

MSstatsLiP is an R/Bioconductor package for statistical analysis of limited
proteolysis mass spectrometry (LiP-MS) experiments, which detect changes in
protein structure — for example, those caused by compound binding or
conformational shifts — in cellular lysates. The workflow models two datasets in
parallel: the structural LiP peptides (from a limited proteolytic digest) and a
trypsin-only (TrP) control that captures overall protein abundance. By adjusting
LiP-peptide changes for the corresponding protein-level changes, MSstatsLiP
distinguishes genuine structural alterations from differences in protein
expression. The package provides functions for summarization, estimation of LiP
peptide abundance, detection of changes across conditions, and specialized
visualizations such as structural and proteolytic-resistance barcode plots.

MSstatsLiP is part of the [MSstats](https://github.com/Vitek-Lab/MSstats)
family of packages, developed and maintained by the
[Vitek Lab](https://olga-vitek-lab.khoury.northeastern.edu/) at Northeastern
University, in collaboration with the Picotti Lab at ETH Zurich. The package and
its documentation are also available at [msstats.org](http://msstats.org).

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MSstatsLiP")
```

The development version can be installed directly from this repository:

```r
BiocManager::install("Vitek-Lab/MSstatsLiP", ref = "devel")
```

## Quick Start

```r
library(MSstatsLiP)

# Example Spectronaut LiP and TrP (trypsin-only control) reports
data(LiPRawData)
data(TrPRawData)
fasta <- system.file("extdata/ExampleFastaFile.fasta", package = "MSstatsLiP")

# Convert Spectronaut output to MSstatsLiP format
formatted <- SpectronauttoMSstatsLiPFormat(LiPRawData, fasta, TrPRawData,
                                            use_log_file = FALSE)

# Summarize feature intensities to peptide/protein level
summarized <- dataSummarizationLiP(formatted, use_log_file = FALSE)

# Test for structural changes across conditions (LiP, TrP, and adjusted models)
model <- groupComparisonLiP(summarized, fasta = fasta, use_log_file = FALSE)

head(model[["LiP.Model"]])
```

Real experiments typically start by converting a search tool's output into
MSstatsLiP format with one of the `*toMSstatsLiPFormat()` converters below, then
proceed with `dataSummarizationLiP()` and `groupComparisonLiP()` as above. See
the workflow vignette for a complete, worked example.

## Supported Converters

MSstatsLiP does not read raw search-tool output directly. Instead, a converter
translates each tool's LiP and TrP reports into MSstatsLiP format before
`dataSummarizationLiP()` is called:

| Search tool / format | Converter function |
| --- | --- |
| Spectronaut | `SpectronauttoMSstatsLiPFormat()` |
| Skyline | `SkylinetoMSstatsLiPFormat()` |
| DIA-NN | `DIANNtoMSstatsLiPFormat()` |

See the [MSstatsLiP workflow vignette](vignettes/MSstatsLiP_Workflow.Rmd) for the
required input files and options for each converter.

## Documentation

- [MSstatsLiP workflow](vignettes/MSstatsLiP_Workflow.Rmd) — full worked example from raw data to results
- [Proteolytic resistance notebook](vignettes/Proteolytic_resistance_notebook.Rmd) — proteolytic-resistance analysis
- [Official website: msstats.org](http://msstats.org)
- [Bioconductor package page and reference manual](https://bioconductor.org/packages/MSstatsLiP)

## Getting Help / Reporting Bugs

- **Questions about usage, statistical methods, or troubleshooting:** please
  post to the [MSstats Google Group](https://groups.google.com/forum/#!forum/msstats).
  This is monitored by the development team and searchable, so it's the fastest
  way to get help and to see if your question has already been answered.
- **Bug reports and feature requests for this repository:** please open a
  [GitHub issue](https://github.com/Vitek-Lab/MSstatsLiP/issues).

## References

If you use MSstatsLiP, please cite:

1. Malinovska L, Cappelletti V, Kohler D, Piazza I, Tsai TH, Pepelnjak M,
   Stalder P, Dörig C, Sesterhenn F, Elsässer F, Kralickova L, Beaton N,
   Reiter L, de Souza N, Vitek O, Picotti P. **Proteome-wide structural changes
   measured with limited proteolysis-mass spectrometry: an advanced protocol for
   high-throughput applications.** *Nat Protoc*. 2023;18(3):659-682.
   [DOI: 10.1038/s41596-022-00771-x](https://doi.org/10.1038/s41596-022-00771-x)

## Funding

MSstats development has been supported by the Chan Zuckerberg Initiative's
[Essential Open Source Software for Science](https://chanzuckerberg.com/eoss/proposals/).

## License

MSstatsLiP is released under the [Artistic-2.0](https://opensource.org/licenses/Artistic-2.0) license.
