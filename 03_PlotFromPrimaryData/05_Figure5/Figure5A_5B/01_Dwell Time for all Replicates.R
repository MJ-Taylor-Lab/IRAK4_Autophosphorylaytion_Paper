library(pacman)
pacman::p_load(dplyr, tidyr, data.table, ggplot2)

# Setting directory paths
Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots"

Input_Table <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots/05_Figure 5/Figure_5A_5B/01_Dwell Time/01_Dwell Time Figure5A_5B_Cell_Summary_by_Image.csv"


# Creating folder to save data --------------------------------------------
# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Save_Path, "05_Figure 5")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Path, "Figure_5A_5B")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}



Cell_Summary_by_Image <- fread(Input_Table)
# Convert SHORT_LABEL to a factor with specified levels
Cell_Summary_by_Image$SHORT_LABEL <- factor(
  Cell_Summary_by_Image$SHORT_LABEL,
  levels = c("D329A/\nT342E/\nT345E/\nS346E", "D329A", "WT")
)


# p-value caluclation -----------------------------------------------------
p_value_test_confirmation = "Y"
if(p_value_test_confirmation == "Y"){
  # Original list with 2-5 components
  Combination_Table <- unique(Cell_Summary_by_Image$SHORT_LABEL)
  
  # Convert the list to a vector since combn works with vectors
  Combination_Table <- unlist(Combination_Table)
  Combination_Table <- sort(Combination_Table)
  
  # Generate unique combinations of 2 components
  Combination_Table <- combn(Combination_Table, 2, simplify = FALSE)
  
  Combination_Function <- function(index){
    Combination <- Combination_Table[[index]] 
    Combination <- Combination %>% as.data.table()
    colnames(Combination) <- c("SHORT_LABEL")
    Combination <- Combination %>% 
      mutate(
        SHORT_LABEL = as.character(SHORT_LABEL)
      )
    
    p_value_Result <- Cell_Summary_by_Image %>% 
      rowwise() %>% 
      mutate(
        TEST = SHORT_LABEL %in% Combination$SHORT_LABEL
      ) %>% 
      filter(
        TEST == TRUE
      )
    
    p_value_Result_1 <- wilcox.test(
      data = p_value_Result,
      FRACTION ~ SHORT_LABEL
    )$p.value
    p_value_Result_1 <- signif(p_value_Result_1, digits = 3)
    
    Temp <- data.table(
      COHORT1 = Combination$SHORT_LABEL[1],
      COHORT2 = Combination$SHORT_LABEL[2],
      p_value = p_value_Result_1
    )
    
    return(Temp)
  }
  
  p_value_Table <- lapply(1:length(Combination_Table), Combination_Function)
  p_value_Table <- rbindlist(p_value_Table)
  
  rm(Combination_Table)
}



# Summarize the data by cohort --------------------------------------------
Cell_Summary_by_Cohort <- Cell_Summary_by_Image %>%
  group_by(
    SHORT_LABEL, 
    COHORT
  ) %>%
  summarise(
    IMAGE_COUNT = n(),  # Count the number of images
    STANDARD_DEVIATION = sd(FRACTION),  # Calculate standard deviation of fractions
    SEM = STANDARD_DEVIATION / IMAGE_COUNT,  # Calculate standard error of the mean
    FRACTION = mean(FRACTION),  # Calculate mean fraction
    Y_MAX = FRACTION + SEM,  # Calculate upper limit
    Y_MIN = case_when(
      FRACTION - SEM < 0 ~ 0,
      TRUE ~ FRACTION - SEM  # Calculate lower limit
    ),
    PRIMARY_PROTIEN_ASSOCIATION_EVENT_COUNT_byCOHORT = max(PRIMARY_PROTIEN_ASSOCIATION_EVENT_COUNT_byCOHORT),  # Number of Primary Protien association events with Complementary Protienby cohort
    PRIMARY_PROTIEN_TRACK_COUNT_byCOHORT = max(PRIMARY_PROTIEN_TRACK_COUNT_byCOHORT)   # Number of Primary Protien tracks
  )





# Plot Primary--------------------------------------------------------------------
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
  # Mean of Colocalisation percentage for 2nd Cohort
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




# Plot with pvalue and other info -----------------------------------------
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
Plot_Mean_SEM_pvalue <- Plot +
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


# Save the plot as a PDF
Plot_Save_Path_1 <- "01_Dwell Time Bin >30s with p_value labels.pdf"
Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Plot_Save_Path_1)
ggsave(
  Plot_Save_Path,
  plot = Plot_Mean_SEM_pvalue,
  height = 3 * 3,
  width = 5 * 4
)

# Plot for Publication ----------------------------------------------------
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

# Save the plot as a PDF
Plot_Save_Path_1 <- "02_Dwell Time Bin >30s.pdf"
Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Plot_Save_Path_1)
ggsave(
  Plot_Save_Path,
  plot = Plot,
  height = 40,
  width = 40,
  units = "mm"
)


# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()
