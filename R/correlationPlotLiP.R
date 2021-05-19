#' Plot run correlation for provided LiP and TrP experiment.
#'
#' @export
#' @import ggplot2
#' @importFrom data.table data.table as.data.table `:=` melt
#' @importFrom grDevices dev.off pdf
#'
#' @param data output of MSstatsLiP converter function. Must include at least
#' ProteinName, Run, and Intensity columns
#' @param method one of "pearson", "kendall", "spearman". Default is pearson.
#' @param value_columns one of "INTENSITY" or "ABUNDANCE". INTENSITY is the raw
#' data, whereas ABUNDANCE is the log transformed INTENSITY column. INTENSITY is
#' default.
#' @param x.axis.size size of axes labels, e.g. name of the comparisons in
#' heatmap, and in comparison plot. Default is 10.
#' @param y.axis.size size of axes labels, e.g. name of targeted proteins in
#' heatmap. Default is 10.
#' @param legend.size size of legend for color at the bottom of volcano plot.
#' Default is 10.
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
#' # add example
correlationPlotLiP <- function(data,
                               method = "pearson",
                               value_columns = "INTENSITY",
                               x.axis.size = 10,
                               y.axis.size = 10,
                               legend.size = 10,
                               width = 10,
                               height = 10,
                               address = ""){


  ##TODO: Add checks and logging
  lip_data <- data$LiP$FeatureLevelData[, c("SUBJECT", value_columns),
                                        with = FALSE]
  runs <- unique(lip_data$SUBJECT)

  ## Create and fill correlation matrix
  cor_mat <- data.table(matrix(0, nrow = length(runs), ncol = length(runs)))
  for (i in seq(runs)){
    for (j in seq(runs)){
      cor_mat[i, j] = cor(lip_data[SUBJECT == runs[i], value_columns, with = FALSE],
                          lip_data[SUBJECT == runs[j], value_columns, with = FALSE],
                          use = "complete.obs", method = method)
    }
  }

  rownames(cor_mat) <- as.character(runs)
  colnames(cor_mat) <- as.character(runs)

  cor_mat <- melt(cor_mat, measure.vars = runs, na.rm = TRUE)
  cor_mat$variable2 <- rep(runs,length(runs))

  min_cor <- min(cor_mat$value)
  max_cor <- max(cor_mat$value)
  mid <- min_cor + ((max_cor - min_cor)/2)

  ## Create PDF to save plots if requested
  if (address != FALSE) {
    allfiles <- list.files()

    num <- 0
    filenaming <- paste0(address, "Correlation_Plot")
    finalfile <- paste0(address, "Correlation_Plot.pdf")

    while (is.element(finalfile, allfiles)) {
      num <- num + 1
      finalfile <- paste0(paste(filenaming, num, sep = "-"), ".pdf")
    }
    pdf(finalfile, width = width, height = height)
  }

  cor_plot <- ggplot(data = cor_mat, aes(variable2, variable, fill = value))+
                    geom_tile(color = "white")+
                scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                    midpoint = mid, limit = c(min_cor, max_cor), space = "Lab",
                    name=paste0(method, "\ncorrelation")) +
                theme_minimal()+
                theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                       size = x.axis.size, hjust = 1),
                axis.text.y = element_text(size = y.axis.size),
                legend.text = element_text(size = legend.size))+
                labs(title = "Run Correlation", x = "", y ="") +
                coord_fixed()
  print(cor_plot)
  if (address != FALSE) {
    dev.off()
  }
}
