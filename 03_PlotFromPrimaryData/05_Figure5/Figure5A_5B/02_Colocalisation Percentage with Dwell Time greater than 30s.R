library(pacman)
pacman::p_load(dplyr, tidyr, data.table, ggplot2)

# Setting directory paths
Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots"

Input_Table <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots/05_Figure 5/Figure_5A_5B/02_Colocalisation/01_1 - Colocalisation Violin Dwell Frames greater than equal to 10 frames or Dwell time greater than or equal to 30 s Figure5A_5B_Colocalisation_Percentage_byCell.csv"


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


# Get tables needed for plot ----------------------------------------------
Colocalisation_Percentage_byCell <- fread(Input_Table)
# Convert SHORT_LABEL to a factor with specified levels
Colocalisation_Percentage_byCell$SHORT_LABEL <- factor(
  Colocalisation_Percentage_byCell$SHORT_LABEL,
  levels = c("D329A/\nT342E/\nT345E/\nS346E", "D329A", "WT")
)


# Colocalisation by Image
Colocalisation_Percentage_byImage <- Colocalisation_Percentage_byCell %>% 
  group_by(
    IMAGE,
    COHORT,
    SHORT_LABEL
  ) %>% 
  summarise(
    MEAN_COLOCLIZED_SPOT_TEST = mean(MEAN_COLOCLIZED_SPOT_TEST)
  ) %>% 
  arrange(
    SHORT_LABEL
  ) %>% 
  as.data.table()


#Calculating Max and Min Number of Cells in replicates to get range
NUMBER_OF_CELLS_PER_REPLICATE_COUNT <- unique(Colocalisation_Percentage_byCell[, .(IMAGE, SHORT_LABEL, CELL)])
NUMBER_OF_CELLS_PER_REPLICATE_COUNT <- NUMBER_OF_CELLS_PER_REPLICATE_COUNT %>% 
  group_by(
    IMAGE,
    SHORT_LABEL
  ) %>% 
  summarise(
    NUMBER_OF_CELLS_PER_REPLICATE = n()
  ) %>% 
  group_by(
    SHORT_LABEL
  ) %>% 
  summarise(
    MIN_NUMBER_OF_CELLS_PER_REPLICATE = min(NUMBER_OF_CELLS_PER_REPLICATE),
    MAX_NUMBER_OF_CELLS_PER_REPLICATE = max(NUMBER_OF_CELLS_PER_REPLICATE),
    TOTAL_NUMBER_OF_CELLS_PER_COHORT = sum(NUMBER_OF_CELLS_PER_REPLICATE)
  )

# COlocalisation Percentage by Cohort
Colocalisation_Percentage_byCOHORT <- Colocalisation_Percentage_byCell %>% 
  arrange(
    SHORT_LABEL,
    IMAGE
  ) %>% 
  transform(
    ID_COHORT = as.numeric(factor(SHORT_LABEL))
  ) %>% 
  group_by(
    ID_COHORT
  ) %>% 
  mutate(
    ID_IMAGE = as.numeric(factor(IMAGE)),
    TOTAL_NUMBER_OF_IMAGES = max(ID_IMAGE)
  ) %>% 
  ungroup() %>% 
  group_by(
    IMAGE,
    COHORT,
    SHORT_LABEL,
    TOTAL_NUMBER_OF_IMAGES
  ) %>% 
  summarise(
    MEAN_COLOCLIZED_SPOT_TEST = mean(MEAN_COLOCLIZED_SPOT_TEST)
  ) %>% 
  ungroup() %>% 
  group_by(
    SHORT_LABEL,
    COHORT
  ) %>%
  summarise(
    TOTAL_NUMBER_OF_IMAGES = mean(TOTAL_NUMBER_OF_IMAGES),
    STANDARD_DEVIATION_OF_MEAN_COLOCLIZED_SPOT_TEST = sd(MEAN_COLOCLIZED_SPOT_TEST),
    STANDARD_ERROR_Of_MEAN_SPOT_TEST = STANDARD_DEVIATION_OF_MEAN_COLOCLIZED_SPOT_TEST/TOTAL_NUMBER_OF_IMAGES,
    MEAN_COLOCLIZED_SPOT_TEST = mean(MEAN_COLOCLIZED_SPOT_TEST),
    YMAX_COLOCLIZED_SPOT_TEST = MEAN_COLOCLIZED_SPOT_TEST + STANDARD_ERROR_Of_MEAN_SPOT_TEST,
    YMIN_COLOCLIZED_SPOT_TEST = case_when(
      MEAN_COLOCLIZED_SPOT_TEST - STANDARD_ERROR_Of_MEAN_SPOT_TEST < 0 ~ 0,
      TRUE ~ MEAN_COLOCLIZED_SPOT_TEST - STANDARD_ERROR_Of_MEAN_SPOT_TEST
    )
  )

Colocalisation_Percentage_byCOHORT <- inner_join(Colocalisation_Percentage_byCOHORT, NUMBER_OF_CELLS_PER_REPLICATE_COUNT, by = "SHORT_LABEL") 

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


Plot_pvalue <-   #p-value Label
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

# Save the reference plot with P-value as a PDF
Plot_Save_Path_1 <- paste0(Dwell_Test_List$SAVE_NAME[Count],".pdf")
Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Plot_Save_Path_1)
ggsave(
  Plot_Save_Path,
  plot = Plot_pvalue,
  height = 3 * 3,
  width = 5 * 4
)



### Plot without stats and axis labels
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

Plot_Save_Path_1 <- paste0(Dwell_Test_List$SAVE_NAME_2[Count],".pdf")
Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Plot_Save_Path_1)
ggsave(
  Plot_Save_Path,
  plot = Plot,
  height = 40,
  width = 32,
  units = "mm"
)
