Kinase_Domain_Only_Violin_Fx <- function(
    Cell_Summary_by_Cell,
    Cell_Summary_by_Cohort,
    Cell_Summary_by_Image
  ){
  Plot <- ggplot(
    data = Cell_Summary_by_Cell
  ) +
    geom_violin(
      scale = "width",
      linewidth = 0.2,
      aes(
        y = SHORT_LABEL, 
        x = NUMBER_OF_EVENTS, 
        color = SHORT_LABEL
      ),
      width = 0.5,
      color = "grey",
      fill = "grey"
    ) +
    geom_errorbar(
      data = Cell_Summary_by_Cohort,
      aes(
        y = SHORT_LABEL, 
        x = MEAN_NUMBER_OF_EVENTS, 
        xmin = MEAN_MINUS_SEM, 
        xmax = MEAN_PLUS_SEM
      ),
      width = 0.1,
      size = 0.2
    ) + 
    geom_crossbar(
      data = Cell_Summary_by_Cohort,
      aes(
        y = SHORT_LABEL, 
        x = MEAN_NUMBER_OF_EVENTS, 
        xmin = MEAN_NUMBER_OF_EVENTS, 
        xmax = MEAN_NUMBER_OF_EVENTS
      ),
      width = 0.1,
      size = 0.1,
      color = "black"
    ) +
    geom_jitter(
      data = Cell_Summary_by_Image,
      aes(
        y = SHORT_LABEL, 
        x = MEAN_NUMBER_OF_EVENTS
      ),
      color = "black",
      fill = Cell_Summary_by_Image$COLOR,
      shape = 21,
      size = 2,
      stroke = 0.2,
      height = 0.1
    ) +
    labs(
      x = "# Number of Events"
    ) +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(LOWER_LIMIT, UPPER_LIMIT),
      breaks = seq(LOWER_LIMIT_AXIS_TICKS, UPPER_LIMIT, by = X_AXIS_INTERVAL), 
      position = "top"
    ) +
    theme_classic() +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 7),
      axis.text = element_text(size = 6),
      legend.position = "none",
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.background = element_blank(),
      axis.line.x.top = element_line(color = "black"),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    )
  
  return(Plot)
}

Kinase_Domain_Only_Violin_pvalue_Fx <- function(
    Plot, 
    Cell_Number_Label,
    p_value_Result
){
  # Add p-value and cell information to the plot
  Plot_pvalue <- Plot +
    annotate(
      "text",
      y = 1.5,
      x = 1100,
      label = p_value_Result,
      size = 3
    ) +
    annotate(
      "text",
      y = 1.1,
      x = 500,
      label = Cell_Number_Label,
      size = 3
    )
  
  return(Plot_pvalue)
}


Kinase_Domain_Only_Lifetime_Fx <- function(
    Cell_Summary_by_Bin
){
  Plot <- ggplot(
    data = Cell_Summary_by_Bin
  ) +
    geom_bar(
      data = Cell_Summary_by_Bin,
      aes(
        y = BIN_COUNT, 
        x = BIN,
        fill = SHORT_LABEL
      ),
      width = 0.2,
      stat = "identity",
      color = "black",
      linewidth = 0.1
    ) +
    scale_fill_manual(
      values = c("#FED09E", "#CC7B16")
    ) +
    facet_rep_wrap(
      ~SHORT_LABEL,
      scales = "fixed",
      repeat.tick.labels = TRUE,
      nrow = FACET_ROW_NUMBERS
    ) +
    scale_x_continuous(
      limits =  c(LOWER_LIMIT,UPPER_LIMIT), ## Sets x axis limts
      breaks = seq(0.3,UPPER_LIMIT, by = AXIS_BREAK_SEQ),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      breaks = pretty_breaks(n = 3)
    ) +
    labs(
      x = X_LABEL,
      y = Y_LABEL
    ) +
    theme_classic()
  
  return(Plot)
}

Kinase_Domain_Only_Lifetime_Mean_Fx <- function(
    Plot, 
    PLOT_MEAN_POSITIION,
    Label
){
  Plot_Mean <- Plot +
    annotate(
      "text",
      y = PLOT_MEAN_POSITIION,
      x = 0.8,
      label = Label,
      size = 3
    ) +
    theme(
      axis.title.y = element_text(size = 10),  # Style Y-axis title
      axis.text = element_text(size = 10),  # Style axis text
      legend.position = "none" # Remove legend
    )

  return(Plot_Mean)
}

Kinase_Domain_Only_Lifetime_Publication_Fx <- function(
    Plot
){
  Plot <- Plot +
    theme(
      axis.title = element_text(size = 6),  # Style Y-axis title
      axis.text = element_text(size = 6),  # Style axis text
      legend.position = "none", # Remove legend
      strip.text.x = element_blank()
    )
  
  return(Plot)
}