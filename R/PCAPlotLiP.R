#' Plot PCA for LiP and TrP datasets. Takes as input LiP data and outputs plot.
#'
#' @export
#' @importFrom tidyr pivot_wider
PCAPlotLiP <- function(data, x.axis.size=10,
                       y.axis.size=10,
                       dot.size=3,
                       text.size=4,
                       text.angle=0,
                       legend.size=13,
                       numProtein=100,
                       width=10,
                       height=10,
                       address=""){

  ## Format Dataset
  data <- msstats_summarized[["LiP"]]$RunlevelData
  data$Cond_rep <- paste(data$GROUP_ORIGINAL, data$SUBJECT, sep = '_')

  ## Convert to long format
  long_data <- pivot_wider(data[,c("Cond_rep", "FULL_PEPTIDE", "LogIntensities")],
                  names_from = Cond_rep, values_from = LogIntensities)

  ## Remove columns which contain NA
  long_data <- long_data[rowSums(is.na(long_data)) == 0, ]
  long_data <- as.data.table(long_data)
  long_data[, FULL_PEPTIDE:=NULL]
  pca_analysis <- prcomp(long_data, center = TRUE, scale. = TRUE)

  test <- summary(pca_analysis)$importance
  test2 <- t(test[2,])
  plot_df <- data.frame(Component = colnames(test2), Variance = test2[1:length(test2)])
  plot_df$Component <- factor(plot_df$Component, levels = plot_df$Component)
  plot_df[1:10,] %>% ggplot() + geom_col(aes(x = Component, y = Variance), fill = "darkblue") +
    labs(title = "PCA Variance by Component", x = "Component", y = "Variance")

  point_df <- pca_analysis$rotation %>% as.data.frame() %>% mutate(Cond_rep = rownames(pca_analysis$rotation))
  point_df <- point_df %>% merge(data %>% distinct(Cond_rep, GROUP_ORIGINAL, SUBJECT),
                                 all.x = TRUE, by = "Cond_rep")

  point_df %>% ggplot(aes(x = PC1, y = PC2, label = SUBJECT, color = GROUP_ORIGINAL)) +
    geom_point(size = 3) + geom_text(aes(label=SUBJECT),hjust=-.4, vjust=-.25, size = 5,
                                     show.legend = FALSE) +
    scale_color_brewer(palette="Dark2") + labs(title = "PCA Analysis", x = "PC1", y = "PC2",
                                               color = "Condition") +
    theme(
      panel.background = element_rect(fill = 'white', colour = "black"),
      legend.key = element_rect(fill = 'white', colour = 'white'),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = 'gray95'),
      axis.text.x = element_text(size = x.axis.size, colour = "black"),
      axis.text.y = element_text(size = y.axis.size, colour = "black"),
      axis.ticks = element_line(colour = "black"),
      axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
      axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
      title = element_text(size = x.axis.size + 8, vjust = 1.5),
      legend.position = "bottom",
      legend.text = element_text(size = legend.size)) +
    guides(color = guide_legend(order = 1,
                                title = "Condition",
                                label.theme = element_text(
                                  size = legend.size, angle = 0),
                                title.vjust = .55))
}
