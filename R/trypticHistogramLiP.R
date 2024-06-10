#' Histogram of Half vs Fully tryptic peptides. Calculates proteotypicity,
#' and then uses calcualtions in histogram.
#'
#' @export
#' @import ggplot2
#' @importFrom data.table as.data.table `:=` setnames copy .N
#' @importFrom dplyr count
#' @importFrom grDevices dev.off pdf
#' @importFrom scales percent
#' @importFrom plotly ggplotly style add_trace plot_ly subplot layout
#'
#' @param data output of MSstatsLiP converter function. Must include at least
#' ProteinName, PeptideSequence, BioReplicate, and Condition columns
#' @param fasta A string of path to a FASTA file, used to match LiP peptides.
#' @param x.axis.size size of x-axis labeling for plot. Default is 10.
#' @param y.axis.size size of y-axis labeling for plot. Default is 10.
#' @param legend.size size of feature legend for half vs fully tryptic peptides
#' below graph. Default is 7.
#' @param width Width of final pdf to be plotted
#' @param height Height of final pdf to be plotted
#' @param color_scale colors of bar chart. Must be one of "bright" or "grey".
#' Default is "bright".
#' @param address the name of folder that will store the results. Default folder
#'  is the current working directory. The other assigned folder has to be
#'  existed under the current working directory. An output pdf file is
#'  automatically created with the default name of "TyrpticPlot.pdf". If
#'  address=FALSE, plot will be not saved as pdf file but shown in window..
#' @param isPlotly Parameter to use Plotly or ggplot2. If set to TRUE, MSstats 
#' will save Plotly plots as HTML files. If set to FALSE MSstats will save ggplot2 plots
#' as PDF files
#' @return plot or pdf
#' @examples
#' # Use output of summarization function
#' trypticHistogramLiP(MSstatsLiP_Summarized,
#'                     system.file("extdata", "ExampleFastaFile.fasta", package="MSstatsLiP"),
#'                     color_scale = "bright", address = FALSE)
#'
trypticHistogramLiP <- function(data, fasta, x.axis.size = 10,
                                y.axis.size = 10, legend.size = 10,
                                width = 12,
                                height = 4,
                                color_scale = "bright",
                                address = "",
                                isPlotly = FALSE) {

  . <- GROUP <- SUBJECT <- fully_TRI <- percent_plot <- NULL

  ## Format input data
  lip.data <- copy(data[["LiP"]]$FeatureLevelData)
  trp.data <- copy(data[["TrP"]]$FeatureLevelData)

  format_fasta <- tidyFasta(fasta)
  format_fasta <- as.data.table(format_fasta)

  ## Add tryptic data
  setnames(lip.data, c("PROTEIN", "PEPTIDE"), c("ProteinName", "PeptideSequence"))
  lip.data$PeptideSequence <- as.character(lip.data$PeptideSequence)
  lip.data$PeptideSequence <- as.character(lapply(lip.data$PeptideSequence, function(x) {substr(x, 1, nchar(x)-2)}))
  tryptic.label <- calculateTrypticity(lip.data, format_fasta)


  plot_df <- merge(lip.data, tryptic.label[, c("ProteinName", "PeptideSequence",
                                               "fully_TRI")], all.x = TRUE,
                   by = c("ProteinName", "PeptideSequence"))

  plot_df[, count := .N, by=.(GROUP, SUBJECT, fully_TRI)]
  plot_df <- unique(plot_df[, c("GROUP", "SUBJECT", "fully_TRI",
                                "count")])

  plot_df[, sum := sum(count), by = .(GROUP, SUBJECT)]
  plot_df$percent_plot <- plot_df$count / plot_df$sum

  if (!isPlotly && address != FALSE) {
    allfiles <- list.files()

    num <- 0
    filenaming <- paste0(address, "TyrpticPlot")
    finalfile <- paste0(address, "TyrpticPlot.pdf")

    while (is.element(finalfile, allfiles)) {
      num <- num + 1
      finalfile <- paste0(paste(filenaming, num, sep = "-"), ".pdf")
    }

    pdf(finalfile, width = width, height = height)
  }

  if (color_scale == "bright"){
    plot_colors <- c("red", "blue")
  } else if (color_scale == "grey"){
    plot_colors <- c("grey32", "grey")
  }

  hist_temp <- ggplot(data = plot_df) +
    geom_col(aes(x = SUBJECT, y = percent_plot, fill = fully_TRI)) +
    facet_wrap(.~GROUP, scales = "free") +
    labs(title = "Proteotrypticity", x = "Replicate", y = "Percent") +
    scale_fill_manual(values = plot_colors, labels = c("Half tryptic (HT)", 
                                                       "Full tryptic (FT)"),
                      name = "Trypticity") +
    scale_y_continuous(labels = percent) +
    theme(
      panel.background = element_rect(fill = 'white', colour = "black"),
      legend.key = element_rect(fill = 'white', colour = 'white'),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = 'gray95'),
      axis.text.y = element_text(size = y.axis.size, colour = "black"),
      axis.ticks = element_line(colour = "black"),
      axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
      axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
      title = element_text(size = x.axis.size + 8, vjust = 1.5),
      legend.position = "bottom",
      legend.text = element_text(size = legend.size))

  print(hist_temp)
  
  if (address != FALSE) {
    dev.off()
  }
  
  if(isPlotly) {
    plotly_plot <- .convertGgplot2Plotly(hist_temp, width = 1000)
    
    # Fix legend
    for (i in seq_along(plotly_plot$x$data)) {
      current_name <- plotly_plot$x$data[[i]]$name
      if (current_name == "TRUE") {
        plotly_plot$x$data[[i]]$name <- "Full tryptic (FT)"
      } else if (current_name == "FALSE") {
        plotly_plot$x$data[[i]]$name <- "Half tryptic (HT)"
      }
    }
    if(address != FALSE) {
      .savePlotlyPlotHTML(list(plotly_plot),address,"TyrpticPlot" ,width, height)
    }
    plotly_plot
  }

}

#' converter for plots from ggplot to plotly
#' @noRd
.convertGgplot2Plotly = function(plot, tips = "all", width = 1800, height = 600) {
  converted_plot <- ggplotly(plot,tooltip = tips)
  converted_plot <- plotly::layout(
    converted_plot,
    width = width,   # Set the width of the chart in pixels
    height = height,  # Set the height of the chart in pixels
    title = list(
      font = list(
        size = 18
      )
    ),
    legend = list(
      x = 0,     # Set the x position of the legend
      y = -0.25,    # Set the y position of the legend (negative value to move below the plot)
      orientation = "h",  # Horizontal orientation
      font = list(
        size = 12  # Set the font size for legend item labels
      ),
      title = list(
        font = list(
          size = 12  # Set the font size for the legend title
        )
      )
    )
  ) 
  converted_plot
}

.savePlotlyPlotHTML = function(plots, address, file_name, width, height) {
  print("Saving plots as HTML")
  pb <- txtProgressBar(min = 0, max = 4, style = 3)
  
  setTxtProgressBar(pb, 1)
  file_name = getFileName(address, file_name, width, height)
  file_name = paste0(file_name,".html")
  
  setTxtProgressBar(pb, 2)
  doc <- .getPlotlyPlotHTML(plots, width, height)
  
  setTxtProgressBar(pb, 3)
  htmltools::save_html(html = doc, file = file_name) # works but lib same folder
  
  setTxtProgressBar(pb, 4)
  zip(paste0(gsub("\\.html$", "", file_name),".zip"), c(file_name, "lib"))
  unlink(file_name)
  unlink("lib",recursive = T)
  
  close(pb)
}

getFileName = function(name_base, file_name, width, height) {
  all_files = list.files(".")
  if(file_name == 'ProfilePlot'){
    num_same_name = sum(grepl(paste0("^", name_base, file_name, "_[0-9]?"), all_files))
  } else {
    num_same_name = sum(grepl(paste0("^", name_base, file_name, "[0-9]?"), all_files))
  }
  if (num_same_name > 0) {
    file_name = paste(file_name, num_same_name + 1, sep = "_")
  }
  file_path = paste0(name_base, file_name)
  return(file_path)
}

.getPlotlyPlotHTML = function(plots, width, height) {
  doc <- htmltools::tagList(lapply(plots,function(x) htmltools::div(x, style = "float:left;width:100%;")))
  # Set a specific width for each plot
  plot_width <- 800
  plot_height <- 600
  
  # Create a div for each plot with style settings
  divs <- lapply(plots, function(x) {
    htmltools::div(x, style = paste0("width:", plot_width, "px; height:", plot_height, "px; margin: 10px;"))
  })
  
  # Combine the divs into a tagList
  doc <- htmltools::tagList(divs)
  doc
}
