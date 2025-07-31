Plot_ColocalisationPercentage <- function(Colocalisation_Percentage_byCell, Colocalisation_Percentage_byImage, Colocalisation_Percentage_byCOHORT){
  Plot <- ggplot(
    data = Colocalisation_Percentage_byCOHORT,
    aes(
      y = SHORT_LABEL,
      x = MEAN_COLOCLIZED_SPOT_TEST
    )
  ) +
    geom_violin(
      data = Colocalisation_Percentage_byCell,
      aes(
        x = MEAN_COLOCLIZED_SPOT_TEST,
        y = SHORT_LABEL,
        fill = SHORT_LABEL
      ),
      width = 0.5,
      scale = "width",
      linewidth = 0.2,
      fill = "grey"
    ) +
    geom_errorbar(
      aes(
        x = MEAN_COLOCLIZED_SPOT_TEST,
        xmin = YMIN_COLOCLIZED_SPOT_TEST,
        xmax = YMAX_COLOCLIZED_SPOT_TEST
      ),
      width = 0.1,
      size = 0.1
    ) +
    # Mean of Colocalisation percentage for 1st Cohort
    geom_segment(
      aes(
        y = 0.95, 
        yend = 1.05,
        x = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[1],
        xend = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[1]
      ),
      size = 0.1
    ) +
    # Mean of Colocalisation percentage for 2nd Cohort
    geom_segment(
      aes(
        y = 1.95, 
        yend = 2.05,
        x = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[2],
        xend = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[2]
      ),
      size = 0.1
    ) +
    # Mean of Colocalisation percentage for 3rd Cohort
    geom_segment(
      aes(
        y = 2.95, 
        yend = 3.05,
        x = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[3],
        xend = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[3]
      ),
      size = 0.1
    ) +
    
    # Add Image Means
    geom_jitter(
      data = Colocalisation_Percentage_byImage,
      aes(
        x = MEAN_COLOCLIZED_SPOT_TEST,
        y = SHORT_LABEL
      ),
      shape = 21,
      fill = Colocalisation_Percentage_byImage$JITTER_COLOR,
      size = 1.2,
      stroke = 0.2,
      height = 0.1  # Control vertical spread
    ) +
    
    #Setting X-axis limts
    scale_x_continuous(
      limits = c(LOWER_LIMIT, UPPER_LIMIT), 
      breaks = seq(LOWER_LIMIT_AXIS_TICKS, UPPER_LIMIT, by = X_AXIS_INTERVAL), 
      position = "top",
      expand = c(0,0)
    ) +
    labs(
      x = "Colocalised puncta per Cell (%)"
    ) +
    theme_classic()
  
  return(Plot)
}

Plot_ColocalisationPercentage_pvalues <- function(Plot, p_value_Table){
  
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
  
  ##### Plot with p-value and axis
  Plot <- Plot +
    #p-value Annotation
    annotate(
      "text", 
      y = 2.5, 
      x = 5, 
      label = p_value_label
    ) +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 7),
      axis.text = element_text(size = 6),
      legend.position = "none",
      axis.text.y = element_text(angle = 90, hjust = 0.5),
      axis.line.x.top = element_line(color = "black") # Add x-axis line on top
    )
  
  return(Plot)
}



Plot_ColocalisationPercentage_forPublication <- function(Plot){
  Plot <- Plot +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 6, hjust = 0.5),
      panel.background = element_blank(), # Remove the panel background
      plot.background = element_blank(),  # Remove the plot background
      legend.background = element_blank(), # Remove the legend background
      axis.line.x.top = element_line(color = "black"), # Add x-axis line on top
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    )
  
  return(Plot)
}