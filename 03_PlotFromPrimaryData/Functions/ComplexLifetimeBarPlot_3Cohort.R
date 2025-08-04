Plot_ComplexLifetime <- function(
    Cell_Summary_by_Image, 
    Cell_Summary_by_Cohort  
  ){
  Plot <- ggplot(
    data = Cell_Summary_by_Image,
    aes(
      fill = SHORT_LABEL
    )
  ) +
    geom_bar(
      data = Cell_Summary_by_Cohort,
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
      data = Cell_Summary_by_Cohort,
      aes(
        x = FRACTION,
        y = SHORT_LABEL
      ),
      fill = "#B2B2B2",
      width = 0.5,
      stat = "identity",
      color = "black",
      linewidth = 0.1
    ) +
    scale_x_continuous(
      limits = c(-0.02, 100),
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


Plot_ComplexLifetime_withMeanSEMandpvalues <- function(
    Cell_Summary_by_Image, 
    Cell_Summary_by_Cohort, 
    Plot, 
    p_value_Result_1  
  ){
  
  # Label containing information of Primary Protein Number and Association Events with Complementary Protein
  PrimaryProtien_Number_and_Association_Events <- paste0(
    Cell_Summary_by_Cohort$SHORT_LABEL[1], " - PRIMARY_PROTIEN_TRACK_COUNT_byCOHORT = ", Cell_Summary_by_Cohort$PRIMARY_PROTIEN_TRACK_COUNT_byCOHORT[1], " \n",
    " and PRIMARY_PROTIEN_ASSOCIATION_EVENT_COUNT_byCOHORT = ", Cell_Summary_by_Cohort$PRIMARY_PROTIEN_ASSOCIATION_EVENT_COUNT_byCOHORT[1], " \n", 
    "FRACTION_MEAN = ", signif(Cell_Summary_by_Cohort$FRACTION[1], digits = 3), " +/- SEM = ", signif(Cell_Summary_by_Cohort$SEM[1], digits = 3) 
  )  
  for(index in 2:nrow(Cell_Summary_by_Cohort)){
    PrimaryProtien_Number_and_Association_Events <- paste0(
      PrimaryProtien_Number_and_Association_Events, " \n", " \n", 
      Cell_Summary_by_Cohort$SHORT_LABEL[index], " - PRIMARY_PROTIEN_TRACK_COUNT_byCOHORT = ", Cell_Summary_by_Cohort$PRIMARY_PROTIEN_TRACK_COUNT_byCOHORT[index], " \n",
      " and PRIMARY_PROTIEN_ASSOCIATION_EVENT_COUNT_byCOHORT = ", Cell_Summary_by_Cohort$PRIMARY_PROTIEN_ASSOCIATION_EVENT_COUNT_byCOHORT[index], " \n", 
      "FRACTION_MEAN = ", signif(Cell_Summary_by_Cohort$FRACTION[index], digits = 3), " +/- SEM = ", signif(Cell_Summary_by_Cohort$SEM[index], digits = 3)
    )
  }
  
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
      data = Cell_Summary_by_Image,
      aes(
        y = SHORT_LABEL, 
        x = FRACTION
      ),
      shape = 21,
      color = "black",
      fill = Cell_Summary_by_Image$COLOR,
      size = 3,
      height = 0.1
    ) +
    geom_errorbar(
      data = Cell_Summary_by_Cohort,
      aes(
        x = FRACTION,
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
    
    # Fill in details about number of Primary Protein Number and Association Events
    annotate(
      "text", 
      y = 1, 
      x = 50, 
      label = PrimaryProtien_Number_and_Association_Events,
      size = 2,
      colour = "red"
    )  +
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

Plot_ComplexLifetime_forPublication <- function(
    Cell_Summary_by_Image, 
    Cell_Summary_by_Cohort, 
    Plot 
  ){
  Plot <- # Modify plot for figure and save as another PDF
    Plot +
    geom_jitter(
      data = Cell_Summary_by_Image,
      aes(
        x = FRACTION,
        y = SHORT_LABEL
      ),
      shape = 21,
      color = "black",
      fill = Cell_Summary_by_Image$COLOR,
      size = 1.5,
      stroke = 0.2,
      height = 0.2  # Control vertical spread
    ) +
    #Adding Mean of 1 Cohort
    geom_segment(
      aes(
        y = 0.95, 
        yend = 1.05,
        x = Cell_Summary_by_Cohort$FRACTION[1],
        xend = Cell_Summary_by_Cohort$FRACTION[1]
      ),
      size = 0.15,
      color = "black"
    ) +
    # Adding Mean + SEM of 1 Cohort
    geom_segment(
      aes(
        y = 0.95, 
        yend = 1.05,
        x = Cell_Summary_by_Cohort$Y_MAX[1],
        xend = Cell_Summary_by_Cohort$Y_MAX[1]
      ),
      size = 0.1,
      color = "black"
    ) +
    # Adding Mean - SEM of 1 Cohort
    geom_segment(
      aes(
        y = 0.95, 
        yend = 1.05,
        x = Cell_Summary_by_Cohort$Y_MIN[1],
        xend = Cell_Summary_by_Cohort$Y_MIN[1]
      ),
      size = 0.1,
      color = "black"
    ) +
    # Adding Line connection Fractiom Mean + SEM  to -SEMof 1 Cohort
    geom_segment(
      aes(
        y = 1, 
        yend = 1,
        x = Cell_Summary_by_Cohort$Y_MAX[1],
        xend = Cell_Summary_by_Cohort$Y_MIN[1]
      ),
      size = 0.1,
      color = "black"
    ) +
    
    
    #Adding Mean of 2 Cohort
    geom_segment(
      aes(
        y = 1.95, 
        yend = 2.05,
        x = Cell_Summary_by_Cohort$FRACTION[2],
        xend = Cell_Summary_by_Cohort$FRACTION[2]
      ),
      size = 0.15,
      color = "black"
    ) +
    # Adding Mean + SEM of 2 Cohort
    geom_segment(
      aes(
        y = 1.95, 
        yend = 2.05,
        x = Cell_Summary_by_Cohort$Y_MAX[2],
        xend = Cell_Summary_by_Cohort$Y_MAX[2]
      ),
      size = 0.1,
      color = "black"
    ) +
    # Adding Mean - SEM of 2 Cohort
    geom_segment(
      aes(
        y = 1.95, 
        yend = 2.05,
        x = Cell_Summary_by_Cohort$Y_MIN[2],
        xend = Cell_Summary_by_Cohort$Y_MIN[2]
      ),
      size = 0.1,
      color = "black"
    ) +
    # Adding Line connection Fraction Mean + SEM  to -SEM of 2 Cohort
    geom_segment(
      aes(
        y = 2, 
        yend = 2,
        x = Cell_Summary_by_Cohort$Y_MAX[2],
        xend = Cell_Summary_by_Cohort$Y_MIN[2]
      ),
      size = 0.1,
      color = "black"
    ) +
    
    #Adding Mean of 3 Cohort
    geom_segment(
      aes(
        y = 2.95, 
        yend = 3.05,
        x = Cell_Summary_by_Cohort$FRACTION[3],
        xend = Cell_Summary_by_Cohort$FRACTION[3]
      ),
      size = 0.15,
      color = "black"
    ) +
    # Adding Mean + SEM of 3 Cohort
    geom_segment(
      aes(
        y = 2.95, 
        yend = 3.05,
        x = Cell_Summary_by_Cohort$Y_MAX[3],
        xend = Cell_Summary_by_Cohort$Y_MAX[3]
      ),
      size = 0.1,
      color = "black"
    ) +
    # Adding Mean - SEM of 3 Cohort
    geom_segment(
      aes(
        y = 2.95, 
        yend = 3.05,
        x = Cell_Summary_by_Cohort$Y_MIN[3],
        xend = Cell_Summary_by_Cohort$Y_MIN[3]
      ),
      size = 0.1,
      color = "black"
    ) +
    # Adding Line connection Fraction Mean + SEM  to -SEM of 3 Cohort
    geom_segment(
      aes(
        y = 3, 
        yend = 3,
        x = Cell_Summary_by_Cohort$Y_MAX[3],
        xend = Cell_Summary_by_Cohort$Y_MIN[3]
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

Plot_ComplexLifetime_forPublication_grey_concentration_dosage <- function(Cell_Summary_by_Image, Cell_Summary_by_Cohort, Plot, Y_LABEL){
  Plot <- # Modify plot for figure and save as another PDF
    Plot +
    geom_jitter(
      data = Cell_Summary_by_Image,
      aes(
        x = FRACTION,
        y = SHORT_LABEL
      ),
      shape = 21,
      color = "black",
      fill = Cell_Summary_by_Image$COLOR,
      size = 1.5,
      stroke = 0.2,
      height = 0.2  # Control vertical spread
    ) +
    # Add errorbar
    geom_errorbar(
      data = Cell_Summary_by_Cohort,
      aes(
        x = FRACTION,
        y = SHORT_LABEL,
        xmin = Y_MIN,
        xmax = Y_MAX
      ),
      width = 0.1,
      size = 0.1
    ) +
    
    #Adding y axis label
    labs(
      y = Y_LABEL
    ) +
    
    theme(
      axis.title.y = element_text(size = 7, angle = 90),
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


