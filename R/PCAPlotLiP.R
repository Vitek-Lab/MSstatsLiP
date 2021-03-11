#' Visualize PCA analysis for LiP and TrP datasets.
#'
#' Takes as input LiP and TrP data from summarization function
#' dataSummarizationLiP. Runs PCA on the summarized data. Can visualize the PCA
#' analysis in three different plots:
#' (1) BarPlot (specify "bar.plot=TRUE" in option bar.plot), to plot a bar plot
#' showing the explained variance per PCA component
#' (2) Peptide/Protein PCA (specify "protein.pca = TRUE" in option protein.pca),
#' to create a dot plot with PCA component 1 and 2 on the axis, for different
#' peptides and proteins.
#' (3) Comparison PCA (specify "comparison.pca = TRUE" in option comparison.pca)
#' , to create a arrow plot with PCA component 1 and 2 on the axis, for
#' different comparisons
#'
#' @export
#' @importFrom tidyr pivot_wider
#' @importFrom data.table as.data.table
#' @importFrom factoextra fviz_eig fviz_pca_ind fviz_pca_var fviz_pca_biplot
#' @importFrom gridExtra grid.arrange
#' @importFrom ggpubr ggpar
#'
#' @param data data name of the list with LiP and (optionally) Protein data, which
#' can be the output of the MSstatsLiP.
#' \code{\link[MSstatsLiP]{dataSummarizationLiP}} function.
#' @param center.pca a logical value indicating whether the variables should be
#' shifted to be zero centered. Alternately, a vector of length equal the number
#' of columns of x can be supplied. The value is passed to scale
#' @param scale.pca a logical value indicating whether the variables should be
#' scaled to have unit variance before the analysis takes place. The default is
#' FALSE for consistency with S, but in general scaling is advisable.
#' Alternatively, a vector of length equal the number of columns of x can be
#' supplied. The value is passed to scale.
#' @param n.components an integer of PCA components to be returned. Default is
#' 10.
#' @param bar.plot a logical value indicating if to visualize PCA bar plot
#' @param protein.pca a logical value indicating if to visualize PCA peptide
#' plot
#' @param comparison.pca a logical value indicating if to visualize PCA
#' comparison plot
#' @param which.pep a list of peptides to be visualized. Default is "all". If
#' too many peptides are plotted the names can overlap.
#' @param which.comparison a list of comparisons to be visualized. Default is
#' "all".
#' @param width width of the saved file. Default is 10.
#' @param height height of the saved file. Default is 10.
#' @param address the name of folder that will store the results. Default
#' folder is the current working directory. The other assigned folder has to
#' be existed under the current working directory. An output pdf file is
#' automatically created with the default name of "VolcanoPlot.pdf" or
#' "Heatmap.pdf". The command address can help to specify where to store the
#' file as well as how to modify the beginning of the file name. If
#' address=FALSE, plot will be not saved as pdf file but showed in window
#' @return plot or pdf
#' @examples
#' Add example
PCAPlotLiP <- function(data,
                       center.pca = TRUE,
                       scale.pca = TRUE,
                       n.components = 10,
                       bar.plot = TRUE,
                       protein.pca = TRUE,
                       comparison.pca = FALSE,
                       which.pep = "all",
                       which.comparison = "all",
                       width=10,
                       height=10,
                       address=""){

  ## TODO: Add input checks
  ## TODO: Add logging

  ## Format Dataset
  lip.data <- data[["LiP"]]$RunlevelData
  trp.data <- data[["TrP"]]$RunlevelData

  ## Run PCA
  lip.pca <- calculate.pc(lip.data, center.pca, scale.pca)
  trp.pca <- calculate.pc(trp.data, center.pca, scale.pca)

  ## Create PDF to save plots if requested
  if (address != FALSE) {
    allfiles <- list.files()

    num <- 0
    filenaming <- paste0(address, "PCA_Plot")
    finalfile <- paste0(address, "PCA_Plot.pdf")

    while (is.element(finalfile, allfiles)) {
      num <- num + 1
      finalfile <- paste0(paste(filenaming, num, sep = "-"), ".pdf")
    }
    pdf(finalfile, width = width, height = height)
  }

  ## Bar plot showing variance of each component
  if (bar.plot){
    lip.bar <- pca.component.bar.plot(lip.pca, n.components)
    trp.bar <- pca.component.bar.plot(trp.pca, n.components)
    grid.arrange(lip.bar, trp.bar, ncol=2)
  }

  ## Bar plot showing variance of each component
  if (which.pep != "all") {
    lip.pca$x <- lip.pca$x[rownames(lip.pca$x) %in% which.pep, ]
    ## TODO: Replace this with more robust Rcpp function
    regex_protein <- '([^-]+)(?:_[^-]+){1}$'
    which.prot <- sapply(which.pep,
                         function(x) {str_match(x, regex_protein)[[2]]})
    which.prot <- unique(as.vector(which.prot))
    trp.pca$x <- trp.pca$x[rownames(trp.pca$x) %in% which.prot, ]

  }
  if (protein.pca){
    lip.prot.plot <- pca.component.prot.plot(lip.pca)
    trp.prot.plot <- pca.component.prot.plot(trp.pca)
    grid.arrange(lip.prot.plot, trp.prot.plot, ncol=1)
  }

  if (which.comparison != "all"){
    lip.pca$rotation <- lip.pca$rotation[
      rownames(lip.pca$rotation) %in% which.comparison, ]
    trp.pca$rotation <- trp.pca$rotation[
      rownames(trp.pca$rotation) %in% which.comparison, ]

  }

  if (comparison.pca){
    lip.comp.plot <- pca.component.comparison.plot(lip.pca)
    trp.comp.plot <- pca.component.comparison.plot(trp.pca)
    grid.arrange(lip.comp.plot, trp.comp.plot, ncol=2)
  }

  if (address != FALSE) {
    dev.off()
  }

}

#' Internal function to run Principle Component Analysis on a dataset. Returns
#' results of PCA.
#' @noRd
calculate.pc <- function(data, center.pca, scale.pca){

  ## Create variables
  data$Cond_rep <- paste(data$GROUP_ORIGINAL,
                         data$SUBJECT, sep = '_')

  if ("FULL_PEPTIDE" %in% colnames(data)){
    ## Convert to long format
    long.format <- pivot_wider(data[,c("Cond_rep", "FULL_PEPTIDE",
                                       "LogIntensities")],
                               names_from = Cond_rep,
                               values_from = LogIntensities)

    ## Remove rows that contain NA
    long.format <- long.format[rowSums(is.na(long.format)) == 0, ]
    long.format <- as.data.table(long.format)
    prot.names <- long.format$FULL_PEPTIDE
    long.format[, FULL_PEPTIDE:=NULL]
  } else {
    ## Convert to long format
    long.format <- pivot_wider(data[,c("Cond_rep", "Protein",
                                       "LogIntensities")],
                               names_from = Cond_rep,
                               values_from = LogIntensities)

    ## Remove rows that contain NA
    long.format <- long.format[rowSums(is.na(long.format)) == 0, ]
    long.format <- as.data.table(long.format)
    prot.names <- long.format$Protein
    long.format[, Protein:=NULL]
  }

  ## Run PCA
  pca_analysis <- prcomp(long.format, center = center.pca, scale. = scale.pca)
  rownames(pca_analysis$x) <- prot.names

  return(pca_analysis)
}

#' Bar plot of explained variance per component
#' @noRd
pca.component.bar.plot <- function(data, n.components){

  temp.bar.plot <- ggpar(fviz_eig(
    data,
    choice = "variance",
    geom = "bar",
    barfill = "steelblue",
    barcolor = "steelblue",
    linecolor = "black",
    ncp = n.components,
    addlabels = TRUE,
    hjust = .5,
    main = "Variance Explained by Component",
    xlab = "PC",
    ylab = "Variance",
    ggtheme = theme_minimal())
  )

  return(temp.bar.plot)
}

#' Dot plot of peptides with top two components on the axis
#' @noRd
pca.component.prot.plot <- function(data){

  temp.bar.plot <- ggpar(
    fviz_pca_ind(data,
                 col.ind = "cos2",
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE)
    )

  return(temp.bar.plot)
}

#' Arrow plot of conditions with top two components on the axis
#' @noRd
pca.component.comparison.plot <- function(data){

  temp.bar.plot <- ggpar(
    fviz_pca_var(data,
                 col.var = "contrib",
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE)
  )

  return(temp.bar.plot)
}

