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




# Plot with pvalue and other info -----------------------------------------
#p-value Label
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
    y = 3, 
    x = 0.55, 
    label = p_value_label,
    size = 2 
  ) +
  
  # Fill in details about number of Primary Protein Number and Association Events
  annotate(
    "text", 
    y = 1, 
    x = 0.5, 
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
