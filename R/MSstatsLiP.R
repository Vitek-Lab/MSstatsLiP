#' MSstatsLiP: A package for identifying and analyzing changes in protein
#' structures caused by compound binding in cellur lysates.
#'
#' A set of tools for detecting differentially abundant LiP peptides in
#' shotgun mass spectrometry-based proteomic experiments. The
#' package includes tools to convert raw data from different spectral processing
#' tools, summarize feature intensities, and fit a linear mixed effects model.
#' If overall protein abundance changes are included, the package will also
#' adjust the LiP peptide fold change for changes in overall protein abundace.
#' Additionally the package includes functionality to plot a variety of data
#' visualizations.
#'
#' @section functions :
#' \itemize{
#'   \item \code{\link{SpectronauttoMSstatsLiPFormat}} : Generates MSstatsLiP
#'   required input format for Spectronaut outputs.
#'   \item \code{\link{trypticHistogramLiP}} : Histogram of Half vs Fully
#'   tryptic peptides. Calculates proteotypicity, and then uses calcualtions in
#'   histogram.
#'   \item \code{\link{correlationPlotLiP}} : Plot run correlation for provided
#'   LiP and TrP experiment.
#'   \item \code{\link{dataSummarizationLiP}} : Summarizes PSM level quantification to
#'   peptide (LiP) and protein level quantification.
#'   \item \code{\link{dataProcessPlotsLiP}} : Visualization for explanatory
#'   data analysis. Specifically gives ability to plot Profile and Quality
#'   Control plots.
#'   \item \code{\link{PCAPlotLiP}} :Visualize PCA analysis for LiP and TrP
#'   datasets. Specifically gives ability to plot explanined varaince per
#'   component, Protein/Peptide PCA, and Condition PCA.
#'   \item \code{\link{groupComparisonLiP}} : Tests for significant changes in
#'   LiP and protein abundance across conditions. Adjusts LiP fold change for
#'   changes in protein abundance.
#'   \item \code{\link{groupComparisonPlotsLiP}} : Visualization for model-based
#'    analysis and summarization.
#'   \item \code{\link{PCAPlotLiP}} : Runs PCA on the summarized data. Can
#'   visualize the PCA analysis in three different plots.
#'   \item \code{\link{BarcodePlotLiP}} :  Shows protein coverage of LiP
#'   modified peptides. Shows significant, insignificant, and missing coverage.
#' }
#'
#' @docType package
#' @name MSstatsLiP
NULL
