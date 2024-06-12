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
#' @importFrom stats prcomp
#' @importFrom plotly ggplotly style add_trace plot_ly subplot layout
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
#' @param isPlotly Parameter to use Plotly or ggplot2. If set to TRUE, MSstats 
#' will save Plotly plots as HTML files. If set to FALSE MSstats will save ggplot2 plots
#' as PDF files
#' @return plot or pdf
#' @examples
#' # Use output of dataSummarizationLiP function
#'
#' # BarPlot
#' PCAPlotLiP(MSstatsLiP_Summarized, bar.plot = TRUE, protein.pca = FALSE)
#'
#' # Protein/Peptide PCA Plot
#' PCAPlotLiP(MSstatsLiP_Summarized, bar.plot = FALSE, protein.pca = TRUE)
#'
#' # Condition PCA Plot
#' PCAPlotLiP(MSstatsLiP_Summarized, bar.plot = FALSE, protein.pca = FALSE,
#'            comparison.pca = TRUE)
#'
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
                       address="",
                       isPlotly = FALSE){


  FULL_PEPTIDE <- Protein <- NULL

  ## Format Dataset
  lip.data <- data[["LiP"]]$ProteinLevelData
  trp.data <- data[["TrP"]]$ProteinLevelData
  if (which.comparison[[1]] != "all"){
    lip.data <- lip.data[lip.data$GROUP %in% which.comparison, ]
    if (!is.null(trp.data)){
      trp.data <- trp.data[trp.data$GROUP %in% which.comparison, ]
    }
  }

  ## Run PCA
  lip.pca <- calculate.pc(lip.data, center.pca, scale.pca)
  if (!is.null(trp.data)){
    trp.pca <- calculate.pc(trp.data, center.pca, scale.pca)
  }

  ## Create PDF to save plots if requested
  if (!isPlotly && address != FALSE) {
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
    lip.bar <- pca.component.bar.plot(lip.pca, n.components,
                                      "LiP Explained Variance")
    if (!is.null(trp.data)){
      trp.bar <- pca.component.bar.plot(trp.pca, n.components,
                                        "TrP Explained Variance")
      grid.arrange(lip.bar, trp.bar, ncol=2)
      bar_plots <- list(LIP = lip.bar, TRP = trp.bar)
    } else {
      print(lip.bar)
      bar_plots <- list(LIP = lip.bar, TRP = NULL)
    }
  }

  ## PCA plot showing variance of each protein
  if (which.pep[[1]] != "all") {
    lip.pca$x <- lip.pca$x[rownames(lip.pca$x) %in% which.pep, , drop=FALSE]

    ## Extract corresponding protein in TrP dataset
    if (!is.null(trp.data)){
      keep <- lip.data[FULL_PEPTIDE %in% which.pep, c("FULL_PEPTIDE",
                                                      "Protein")]
      which.prot <- unique(keep[, Protein])
      trp.pca$x <- trp.pca$x[rownames(trp.pca$x) %in% which.prot, , drop=FALSE]
    }
  }

  if (protein.pca){
    lip.prot.plot <- pca.component.prot.plot(lip.pca, "LiP Peptide PCA")
    if (!is.null(trp.data)){
      if (nrow(trp.pca$x) > 1){
        trp.prot.plot <- pca.component.prot.plot(trp.pca, "TrP Protein PCA")
        grid.arrange(lip.prot.plot, trp.prot.plot, ncol=1)
        protein_plots <- list(LIP = lip.prot.plot, TRP = trp.prot.plot)
      } else {
        print(lip.prot.plot)
        protein_plots <- list(LIP = lip.prot.plot, TRP = NULL)
      }
    } else {
      print(lip.prot.plot)
      protein_plots <- list(LIP = lip.prot.plot, TRP = NULL)
    }
  }

  if (comparison.pca){
    lip.comp.plot <- pca.component.comparison.plot(lip.pca, "LiP Component PCA")
    if (!is.null(trp.data)){
      trp.comp.plot <- pca.component.comparison.plot(trp.pca,
                                                     "TrP Component PCA")
      grid.arrange(lip.comp.plot, trp.comp.plot, ncol=2)
      comparison_plots <- list(LIP = lip.comp.plot, TRP = trp.comp.plot)
    } else {
      print(lip.comp.plot)
      comparison_plots <- list(LIP = lip.comp.plot, TRP = NULL)
    }
  }

  if (address != FALSE) {
    dev.off()
  }
  
  if(isPlotly) {
    plotly_plots <- list()
    
    # Bar Plots
    if(bar.plot) {
      plotly_plot_bar_lip <- .convertGgplot2Plotly(bar_plots[["LIP"]], width = 1300)
      plotly_plot_bar_lip$x$data[[2]]$textposition='top'
      if (!is.null(bar_plots[["TRP"]])) {
        plotly_plot_bar_trp <- .convertGgplot2Plotly(bar_plots[["TRP"]], width = 1300)
        plotly_plot_bar_trp$x$data[[2]]$textposition='top'
        plotly_plot_combined <- .combineSubPlotsPlotly(plotly_plot_bar_lip, plotly_plot_bar_trp)
        plotly_plots <- c(plotly_plots, list(plotly_plot_combined))
      } else {
        plotly_plots <- c(plotly_plots, list(plotly_plot_bar_lip))
      }
    }
    
    if(protein.pca) {
      plotly_plot_prot_lip <- .convertGgplot2Plotly(protein_plots[["LIP"]], width = 1350)
      # plotly_plot_bar_lip$x$data[[2]]$textposition='top'
      if (!is.null(protein_plots[["TRP"]])) {
        plotly_plot_prot_trp <- .convertGgplot2Plotly(protein_plots[["TRP"]], width = 1350)
        # plotly_plot_bar_trp$x$data[[2]]$textposition='top'
        plotly_plot_combined <- .combineSubPlotsPlotly(plotly_plot_prot_lip, plotly_plot_prot_trp)
        plotly_plots <- c(plotly_plots, list(plotly_plot_combined))
      } else {
        plotly_plots <- c(plotly_plots, list(plotly_plot_prot_lip))
      }
    }
    
    if(comparison.pca) {
      
    }
    
    if(address != FALSE) {
      .savePlotlyPlotHTML(plotly_plots,address,"PCA_Plot" ,width, height)
    }
    plotly_plots <- unlist(plotly_plots, recursive = FALSE)
    plotly_plots
  }

}

#' Internal function to run Principle Component Analysis on a dataset. Returns
#' results of PCA.
#' @noRd
calculate.pc <- function(data, center.pca, scale.pca){

  Cond_rep <- LogIntensities <- FULL_PEPTIDE <- Protein <- NULL

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
pca.component.bar.plot <- function(data, n.components, title){

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
    main = title,
    xlab = "PC",
    ylab = "Variance",
    ggtheme = theme_minimal())
  )
  
  return(temp.bar.plot)
}

#' Dot plot of peptides with top two components on the axis
#' @noRd
pca.component.prot.plot <- function(data, title){
  print("=========")
  print(data)
  print("=========")
  temp.bar.plot <- ggpar(
    fviz_pca_ind(data,
                 col.ind = "cos2",
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE,
                 title = title)
    )

  return(temp.bar.plot)
}

#' Arrow plot of conditions with top two components on the axis
#' @noRd
pca.component.comparison.plot <- function(data, title){

  temp.bar.plot <- ggpar(
    fviz_pca_var(data,
                 geom = c("point", "text"),
                 col.var = "contrib",
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE,
                 title = title)
  )

  return(temp.bar.plot)
}

.combineSubPlotsPlotly = function(plotly_plot_ptm, plotly_plot_protein) {
  title_ptm <- plotly_plot_ptm$x$layout$title$text
  title_protein <- plotly_plot_protein$x$layout$title$text
  
  plotly_plot_ptm <- plotly::layout(plotly_plot_ptm, title = "")
  plotly_plot_protein <- plotly::layout(plotly_plot_protein, title = "")
  
  plotly_plot_combined <- subplot(plotly_plot_ptm, plotly_plot_protein,  margin=0.04,titleX = TRUE, titleY = TRUE)
  plotly_plot_combined <- plotly::layout(plotly_plot_combined,
                                         annotations = list(
                                           list(
                                             x = 0.25,  # Centered horizontally over the first plot
                                             y = 1,  # Above the first plot
                                             text = title_ptm,  # Title for the first plot
                                             showarrow = FALSE,
                                             xref = 'paper',
                                             yref = 'paper',
                                             xanchor = 'center',
                                             yanchor = 'bottom',
                                             font = list(size = 16)
                                           ),
                                           list(
                                             x = 0.75,  # Centered horizontally over the second plot
                                             y = 1,  # Above the second plot
                                             text = title_protein,  # Title for the second plot
                                             showarrow = FALSE,
                                             xref = 'paper',
                                             yref = 'paper',
                                             xanchor = 'center',
                                             yanchor = 'bottom',
                                             font = list(size = 16)
                                           )
                                         ))
  plotly_plot_combined
}

