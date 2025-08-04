Plot_Size <- function(
    Size_Summary_ByImage, 
    Size_Summary_ByCohort  
){
  Plot <- ggplot(
    data = Size_Summary_ByImage,
    aes(
      fill = SHORT_LABEL
    )
  ) +
    geom_bar(
      data = Size_Summary_ByCohort,
      aes(
        y = SHORT_LABEL, 
        x = 100
      ),
      fill = "#E5E5E5",
      width = 0.5,
      stat = "identity",
      color = "black",
      linewidth = 0.1
    ) +
    geom_bar(
      data = Size_Summary_ByCohort,
      aes(
        x = MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD,
        y = SHORT_LABEL
      ),
      fill = "#B2B2B2",
      width = 0.5,
      stat = "identity",
      color = "black",
      linewidth = 0.1
    ) +
    scale_x_continuous(
      limits = c(-1, 100),
      expand = c(0, 0),
      breaks = c(0, 50, 100),
      position = "top",
    ) +
    labs(
      x = "Percentage"
    ) +
    theme_classic()
  
  return(Plot)
}


Plot_Size_withMeanSEMandpvalues <- function(
    Size_Summary_ByImage, 
    Size_Summary_ByCohort, 
    Plot, 
    p_value_Table  
){
  #p-value Label
  p_value_label <- paste0(
    p_value_Table$COHORT1[1], " vs ", p_value_Table$COHORT2[1], " p_value = ", p_value_Table$p_value[1]
  )
  
  for(index in 2:nrow(p_value_Table)){
    p_value_label <- paste0(
      p_value_label, " \n",
      p_value_Table$COHORT1[index], " vs ", p_value_Table$COHORT2[index], " p_value = ", p_value_Table$p_value[index]
    )
  }
  
  Plot <- Plot +
    geom_jitter(
      data = Size_Summary_ByImage,
      aes(
        y = SHORT_LABEL, 
        x = MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD
      ),
      shape = 21,
      color = "black",
      fill = Size_Summary_ByImage$COLOR,
      size = 3,
      height = 0.1
    ) +
    geom_errorbar(
      data = Size_Summary_ByCohort,
      aes(
        x = MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD,
        y = SHORT_LABEL,
        xmin = Y_MIN,
        xmax = Y_MAX
      ),
      width = 0.1
    ) +
    
    #p-value Annotation
    annotate(
      "text", 
      y = 2, 
      x = 55, 
      label = p_value_label,
      size = 2 
    ) +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 7),
      axis.text = element_text(size = 6),
      legend.position = "none",
      axis.text.y = element_text(hjust = 0.5),
      axis.line.x.top = element_line(color = "black") # Add x-axis line on top
    )
  
  return(Plot)
}

Plot_Size_forPublication <- function(
    Size_Summary_ByImage, 
    Size_Summary_ByCohort, 
    Plot 
){
  Plot <- # Modify plot for figure and save as another PDF
    Plot +
    geom_jitter(
      data = Size_Summary_ByImage,
      aes(
        x = MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD,
        y = SHORT_LABEL
      ),
      shape = 21,
      color = "black",
      fill = Size_Summary_ByImage$COLOR,
      size = 1.5,
      stroke = 0.2,
      height = 0.2  # Control vertical spread
    ) +
    #Adding Mean of 1 Cohort
    geom_segment(
      aes(
        y = 0.95, 
        yend = 1.05,
        x = Size_Summary_ByCohort$MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD[1],
        xend = Size_Summary_ByCohort$MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD[1]
      ),
      size = 0.15,
      color = "black"
    ) +
    # Adding Mean + SEM of 1 Cohort
    geom_segment(
      aes(
        y = 0.95, 
        yend = 1.05,
        x = Size_Summary_ByCohort$Y_MAX[1],
        xend = Size_Summary_ByCohort$Y_MAX[1]
      ),
      size = 0.1,
      color = "black"
    ) +
    # Adding Mean - SEM of 1 Cohort
    geom_segment(
      aes(
        y = 0.95, 
        yend = 1.05,
        x = Size_Summary_ByCohort$Y_MIN[1],
        xend = Size_Summary_ByCohort$Y_MIN[1]
      ),
      size = 0.1,
      color = "black"
    ) +
    # Adding Line connection Fractiom Mean + SEM  to -SEMof 1 Cohort
    geom_segment(
      aes(
        y = 1, 
        yend = 1,
        x = Size_Summary_ByCohort$Y_MAX[1],
        xend = Size_Summary_ByCohort$Y_MIN[1]
      ),
      size = 0.1,
      color = "black"
    ) +
    
    
    #Adding Mean of 2 Cohort
    geom_segment(
      aes(
        y = 1.95, 
        yend = 2.05,
        x = Size_Summary_ByCohort$MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD[2],
        xend = Size_Summary_ByCohort$MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD[2]
      ),
      size = 0.15,
      color = "black"
    ) +
    # Adding Mean + SEM of 2 Cohort
    geom_segment(
      aes(
        y = 1.95, 
        yend = 2.05,
        x = Size_Summary_ByCohort$Y_MAX[2],
        xend = Size_Summary_ByCohort$Y_MAX[2]
      ),
      size = 0.1,
      color = "black"
    ) +
    # Adding Mean - SEM of 2 Cohort
    geom_segment(
      aes(
        y = 1.95, 
        yend = 2.05,
        x = Size_Summary_ByCohort$Y_MIN[2],
        xend = Size_Summary_ByCohort$Y_MIN[2]
      ),
      size = 0.1,
      color = "black"
    ) +
    # Adding Line connection PERCENTAGE Mean + SEM  to -SEM of 2 Cohort
    geom_segment(
      aes(
        y = 2, 
        yend = 2,
        x = Size_Summary_ByCohort$Y_MAX[2],
        xend = Size_Summary_ByCohort$Y_MIN[2]
      ),
      size = 0.1,
      color = "black"
    ) +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none",
      axis.text.y = element_text(size = 6, hjust = 0),
      axis.text.x = element_text(size = 6, hjust = 0),
      panel.background = element_blank(), # Remove the panel background
      plot.background = element_blank(),  # Remove the plot background
      legend.background = element_blank(), # Remove the legend background
      axis.line.x.top = element_line(color = "black"), # Add x-axis line on top
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    )
  
  return(Plot)
}