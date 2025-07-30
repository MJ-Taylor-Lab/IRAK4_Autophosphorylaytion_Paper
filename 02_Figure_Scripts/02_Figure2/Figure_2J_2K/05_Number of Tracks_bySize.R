# Load required libraries using pacman (auto-installs if missing)
pacman::p_load(ggbeeswarm, ggforce, gridExtra, writexl)

# Set the subdirectory where the plots will be saved
Plot_Directory_Save_Path_Sub_Dir <- file.path(Plot_Directory_Path, "05_Number of Tracks_bySize")

# Create the directory if it does not already exist
if(!file.exists(Plot_Directory_Save_Path_Sub_Dir)){
  dir.create(Plot_Directory_Save_Path_Sub_Dir)
}

Y_AXIS_LABEL = "# of IRAK2 (per cell)"  # Y-axis label for the plots
YAXIS_LOWER_LIMIT_LESS_THAN_PLOT = 0  # Lower y-axis limit for "<4" size category
YAXIS_UPPER_LIMIT_LESS_THAN_PLOT = 750  # Upper y-axis limit for "<4"
YAXIS_LOWER_LIMIT_GREATER_THAN_PLOT = 0  # Lower y-axis limit for ">=4"
YAXIS_UPPER_LIMIT_GREATER__THAN_PLOT = 100  # Upper y-axis limit for ">=4"



# DATA PROCESSING-----------------------------------------------
# Finding number of IRAK2 spots per cell that are greater or less than size of Threshold of 4
IRAK2_Size_Summary <- Table %>% 
  #Remove rows with missing complementary protien data
  filter(
    !is.na(COMPLEMENTARY_NORMALIZED_INTENSITY_1)
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  # Selecting for IRAK4 tracks that recruit even a little bit of complementary protein
  filter(
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 0.75
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID,
    CELL,
    IMAGE,
    SHORT_LABEL
  ) %>% 
  # Max Size of Complementary protein recruited
  summarise(
    MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1 = max(COMPLEMENTARY_NORMALIZED_INTENSITY_1)
  ) %>% 
  mutate(
    SIZE_CATEGORY = case_when(
      MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1 < 4 ~ "<4",
      TRUE ~ ">=4"
    ),
    IMAGE_CELL = paste0(IMAGE, "_", CELL)
  ) %>% 
  group_by(
    IMAGE_CELL,
    IMAGE,
    SHORT_LABEL,
    SIZE_CATEGORY
  ) %>% 
  summarise(
    NUMBER_OF_TRACKS = n()
  ) %>% 
  group_by(
    IMAGE_CELL
  ) %>% 
  mutate(
    TOTAL_NUMBER_OF_TRACKS = sum(NUMBER_OF_TRACKS)
  ) %>% 
  rowwise() %>% 
  mutate(
    PERCENTAGE_OF_TRACKS_IN_COHORT = NUMBER_OF_TRACKS*100/TOTAL_NUMBER_OF_TRACKS
  )


# Summary statistics by image
IRAK2_Size_Summary_ByImage <- IRAK2_Size_Summary %>% 
  group_by(
    IMAGE,
    SIZE_CATEGORY,
    SHORT_LABEL
  ) %>% 
  summarise(
    MEAN_NUMBER_OF_TRACKS = mean(NUMBER_OF_TRACKS),
    MEAN_NUMBER_OF_TRACKS = as.numeric(MEAN_NUMBER_OF_TRACKS),
    MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT = mean(PERCENTAGE_OF_TRACKS_IN_COHORT)
  ) %>% 
  mutate(
    COLOR = case_when(
      SHORT_LABEL == "IRAK2-WT" ~ "#33AFCC",
      SHORT_LABEL == "IRAK2-DD" ~ "#F99108"
    ),
    COHORT = paste0(SHORT_LABEL, "_", SIZE_CATEGORY)
  ) 

# Summary statistics by cohort
IRAK2_Size_Summary_ByCohort <- IRAK2_Size_Summary_ByImage %>% 
  arrange(
    SHORT_LABEL,
    IMAGE
  ) %>% 
  # Counting number of replicates per cohort
  transform(
    ID_COHORT = as.numeric(factor(SHORT_LABEL))
  ) %>% 
  group_by(
    ID_COHORT
  ) %>% 
  # Assigning highest number as number of replicates per cohort
  mutate(
    ID_IMAGE = as.numeric(factor(IMAGE)),
    TOTAL_NUMBER_OF_IMAGES = max(ID_IMAGE)
  ) %>% 
  ungroup() %>% 
  group_by(
    SIZE_CATEGORY,
    SHORT_LABEL,
    COLOR,
    COHORT
  ) %>% 
  summarise(
    TOTAL_NUMBER_OF_IMAGES = mean(TOTAL_NUMBER_OF_IMAGES),
    #Finding SD, SEM, MEAN, YMAX and YMIN for absolute number of tracks per cell
    STANDARD_DEVIATION_OF_MEAN_OF_MEAN_NUMBER_OF_TRACKS = sd(MEAN_NUMBER_OF_TRACKS),
    STANDARD_ERROR_Of_MEAN_MEAN_OF_MEAN_NUMBER_OF_TRACKS = STANDARD_DEVIATION_OF_MEAN_OF_MEAN_NUMBER_OF_TRACKS/TOTAL_NUMBER_OF_IMAGES,
    MEAN_OF_MEAN_NUMBER_OF_TRACKS = mean(MEAN_NUMBER_OF_TRACKS),
    YMAX_MEAN_OF_MEAN_NUMBER_OF_TRACKS = MEAN_OF_MEAN_NUMBER_OF_TRACKS + STANDARD_ERROR_Of_MEAN_MEAN_OF_MEAN_NUMBER_OF_TRACKS,
    YMIN_MEAN_OF_MEAN_NUMBER_OF_TRACKS = case_when(
      MEAN_OF_MEAN_NUMBER_OF_TRACKS - STANDARD_ERROR_Of_MEAN_MEAN_OF_MEAN_NUMBER_OF_TRACKS < 0 ~ 0,
      TRUE ~ MEAN_OF_MEAN_NUMBER_OF_TRACKS - STANDARD_ERROR_Of_MEAN_MEAN_OF_MEAN_NUMBER_OF_TRACKS
    ),
    #Finding SD, SEM, MEAN, YMAX and YMIN for percentage of tracks per cell
    STANDARD_DEVIATION_OF_MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT = sd(MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT),
    STANDARD_ERROR_Of_MEAN_MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT = STANDARD_DEVIATION_OF_MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT/TOTAL_NUMBER_OF_IMAGES,
    MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT = mean(MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT),
    YMAX_MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT = MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT + STANDARD_ERROR_Of_MEAN_MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT,
    YMIN_MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT = case_when(
      MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT - STANDARD_ERROR_Of_MEAN_MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT < 0 ~ 0,
      TRUE ~ MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT - STANDARD_ERROR_Of_MEAN_MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT
    )
  )



# Subset into <4 and >=4 groups for plotting and stats
# p-Value for less than threshold and seprating tables for plotting-----------------------------------------------------------------
IRAK2_Size_Summary_LessThan <- IRAK2_Size_Summary %>% 
  filter(
    SIZE_CATEGORY == "<4"
  )

IRAK2_Size_Summary_LessThan$SHORT_LABEL <- 
  factor(
    IRAK2_Size_Summary_LessThan$SHORT_LABEL,
    levels = c("IRAK2-WT", "IRAK2-DD")
  )

IRAK2_Size_Summary_ByCohort_LessThan <- IRAK2_Size_Summary_ByCohort %>% 
  filter(
    SIZE_CATEGORY == "<4"
  )

IRAK2_Size_Summary_ByCohort_LessThan$SHORT_LABEL <- 
  factor(
    IRAK2_Size_Summary_ByCohort_LessThan$SHORT_LABEL,
    levels = c("IRAK2-WT", "IRAK2-DD")
  )

IRAK2_Size_Summary_ByImage_LessThan <- IRAK2_Size_Summary_ByImage %>% 
  filter(
    SIZE_CATEGORY == "<4"
  )

p_value_test_confirmation = "Y"
if(p_value_test_confirmation == "Y"){
  # Original list with 2-5 components
  Combination_Table <- unique(IRAK2_Size_Summary_ByImage_LessThan$COHORT)
  
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
    
    p_value_Result <- IRAK2_Size_Summary_ByImage_LessThan %>% 
      rowwise() %>% 
      mutate(
        TEST = COHORT %in% Combination$SHORT_LABEL
      ) %>% 
      filter(
        TEST == TRUE
      )
    
    p_value_Result_1 <- wilcox.test(
      data = p_value_Result,
      MEAN_NUMBER_OF_TRACKS ~ COHORT
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


p_value_Table_Complete <- p_value_Table

# p-Value for greater than threshold-----------------------------------------------------------------
IRAK2_Size_Summary_GreaterThan <- IRAK2_Size_Summary %>% 
  filter(
    SIZE_CATEGORY != "<4"
  )

IRAK2_Size_Summary_GreaterThan$SHORT_LABEL <- 
  factor(
    IRAK2_Size_Summary_GreaterThan$SHORT_LABEL,
    levels = c("IRAK2-WT", "IRAK2-DD")
  )

IRAK2_Size_Summary_GreaterThan$SHORT_LABEL <- 
  factor(
    IRAK2_Size_Summary_GreaterThan$SHORT_LABEL,
    levels = c("IRAK2-WT", "IRAK2-DD")
  )

IRAK2_Size_Summary_ByCohort_GreaterThan <- IRAK2_Size_Summary_ByCohort %>% 
  filter(
    SIZE_CATEGORY != "<4"
  )

IRAK2_Size_Summary_ByImage_GreaterThan <- IRAK2_Size_Summary_ByImage %>% 
  filter(
    SIZE_CATEGORY != "<4"
  ) %>% 
  arrange(
    SHORT_LABEL
  )


p_value_test_confirmation = "Y"
if(p_value_test_confirmation == "Y"){
  # Original list with 2-5 components
  Combination_Table <- unique(IRAK2_Size_Summary_ByImage_GreaterThan$COHORT)
  
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
    
    p_value_Result <- IRAK2_Size_Summary_ByImage_GreaterThan %>% 
      rowwise() %>% 
      mutate(
        TEST = COHORT %in% Combination$SHORT_LABEL
      ) %>% 
      filter(
        TEST == TRUE
      )
    
    p_value_Result_1 <- wilcox.test(
      data = p_value_Result,
      MEAN_NUMBER_OF_TRACKS ~ COHORT
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

# combine results
p_value_Table_Complete <- rbind(p_value_Table_Complete, p_value_Table)


rm(
  IRAK2_Size_Summary,
  IRAK2_Size_Summary_ByCohort,
  IRAK2_Size_Summary_ByImage,
  p_value_Table
)



# # Plot Less Than --------------------------------------------------------------------
# Plot_Less_Than <- ggplot(
#   data = IRAK2_Size_Summary_ByCohort_LessThan
# ) +
#   geom_quasirandom(
#     data = IRAK2_Size_Summary_LessThan,
#     aes(
#       x = SHORT_LABEL,
#       y = NUMBER_OF_TRACKS,
#     ),
#     fill = "grey",
#     size = 0.5,
#     width = 0.2,
#     alpha = 0.5,
#     stroke = 0
#   ) +
#   # Adding Error Bar
#   geom_errorbar(
#     aes(
#       x = SHORT_LABEL,  # Add this line to specify the x aesthetic
#       y = MEAN_OF_MEAN_NUMBER_OF_TRACKS,
#       ymin = YMIN_MEAN_OF_MEAN_NUMBER_OF_TRACKS,
#       ymax = YMAX_MEAN_OF_MEAN_NUMBER_OF_TRACKS
#     ),
#     width = 0.1,
#     size = 0.1,
#     color = "black" 
#   ) +
#   # Adding Mean Line
#   geom_crossbar(
#     aes(
#       x = SHORT_LABEL,
#       y = MEAN_OF_MEAN_NUMBER_OF_TRACKS,
#       ymin = MEAN_OF_MEAN_NUMBER_OF_TRACKS,
#       ymax = MEAN_OF_MEAN_NUMBER_OF_TRACKS
#     ),
#     width = 0.1, 
#     size = 0.1, 
#     color = "black" 
#   ) +
#   # Add Image Means
#   geom_jitter(
#     data = IRAK2_Size_Summary_ByImage_LessThan,
#     aes(
#       y = MEAN_NUMBER_OF_TRACKS,
#       x = SHORT_LABEL
#     ),
#     shape = 21,
#     fill = IRAK2_Size_Summary_ByImage_LessThan$COLOR,
#     size = 1.2,
#     stroke = 0.2,
#     width = 0.1  # Control vertical spread
#   ) +
#   labs(
#     y = Y_AXIS_LABEL
#   ) +
#   scale_y_continuous(
#     limits = c(YAXIS_LOWER_LIMIT_LESS_THAN_PLOT, YAXIS_UPPER_LIMIT_LESS_THAN_PLOT),
#     expand = c(0,0)
#   ) +
#   theme_classic() +
#   theme(
#     legend.position = "none",
#     axis.text = element_text(size = 6, hjust = 0),
#     axis.title.x = element_blank(),
#     panel.background = element_blank(), # Remove the panel background
#     plot.background = element_blank(),  # Remove the plot background
#     legend.background = element_blank(), # Remove the legend background
#     axis.line = element_line(color = "black") # Add x-axis line on top
#   )
# 
# 
# # Plot Greater Than -------------------------------------------------------
# Plot_Greater_Than <- ggplot(
#   data = IRAK2_Size_Summary_ByCohort_GreaterThan
# ) +
#   geom_quasirandom(
#     data = IRAK2_Size_Summary_GreaterThan,
#     aes(
#       x = SHORT_LABEL,
#       y = NUMBER_OF_TRACKS,
#     ),
#     fill = "grey",
#     size = 0.5,
#     width = 0.2,
#     alpha = 0.5,
#     stroke = 0
#   ) +
#   # Adding Error Bar
#   geom_errorbar(
#     aes(
#       x = SHORT_LABEL,  # Add this line to specify the x aesthetic
#       y = MEAN_OF_MEAN_NUMBER_OF_TRACKS,
#       ymin = YMIN_MEAN_OF_MEAN_NUMBER_OF_TRACKS,
#       ymax = YMAX_MEAN_OF_MEAN_NUMBER_OF_TRACKS
#     ),
#     width = 0.1,
#     size = 0.1,
#     color = "black" 
#   ) +
#   # Adding Mean Line
#   geom_crossbar(
#     aes(
#       x = SHORT_LABEL,
#       y = MEAN_OF_MEAN_NUMBER_OF_TRACKS,
#       ymin = MEAN_OF_MEAN_NUMBER_OF_TRACKS,
#       ymax = MEAN_OF_MEAN_NUMBER_OF_TRACKS
#     ),
#     width = 0.1, 
#     size = 0.1, 
#     color = "black" 
#   ) +
#   # Add Image Means
#   geom_jitter(
#     data = IRAK2_Size_Summary_ByImage_GreaterThan,
#     aes(
#       y = MEAN_NUMBER_OF_TRACKS,
#       x = SHORT_LABEL
#     ),
#     shape = 21,
#     fill = IRAK2_Size_Summary_ByImage_GreaterThan$COLOR,
#     size = 1.2,
#     stroke = 0.2,
#     width = 0.1  # Control vertical spread
#   ) +
#   scale_y_continuous(
#     limits = c(YAXIS_LOWER_LIMIT_GREATER_THAN_PLOT, YAXIS_UPPER_LIMIT_GREATER__THAN_PLOT),
#     expand = c(0,0)
#   ) +
#   labs(
#     y = Y_AXIS_LABEL
#   ) +
#   theme_classic() +
#   theme(
#     legend.position = "none",
#     axis.text = element_text(size = 6, hjust = 0),
#     axis.title.x = element_blank(),
#     # axis.title.y = element_blank(),
#     panel.background = element_blank(), # Remove the panel background
#     plot.background = element_blank(),  # Remove the plot background
#     legend.background = element_blank(), # Remove the legend background
#     axis.line = element_line(color = "black") # Add x-axis line on top,
#   )
# 
# 
# 
# 
# 
# # Combine plots -----------------------------------------------------------
# # Arrange the plots side by side
# Plot <- grid.arrange(Plot_Less_Than, Plot_Greater_Than, ncol = 2)
# 
# 
# Plot_Save_Path_1 <- paste0("01_NumberOfComplementaryProtein_vs_Condition_by_Cohort_with_p_value.pdf")
# Plot_Save_Path <- file.path(Plot_Directory_Save_Path_Sub_Dir, Plot_Save_Path_1)
# 
# 
# ggsave(
#   Plot_Save_Path,
#     plot = Plot,
#     height = 40,
#     width = 60,
#     units = "mm"
# )
# 
# 

# Partial clean up --------------------------------------------------------
rm(
  Plot_Save_Path_1,
  Plot_Save_Path,
  Plot,
  Plot_Less_Than,
  Plot_Greater_Than,
  Y_AXIS_LABEL,
  YAXIS_LOWER_LIMIT_LESS_THAN_PLOT,
  YAXIS_UPPER_LIMIT_LESS_THAN_PLOT,
  YAXIS_LOWER_LIMIT_GREATER_THAN_PLOT,
  YAXIS_UPPER_LIMIT_GREATER__THAN_PLOT,
  Combination_Function,
  p_value_test_confirmation
)


# Percentage Plot p value calculation---------------------------------------------------------
p_value_Table_Complete <- p_value_Table_Complete %>% 
  mutate(
    COMPARISON_DATA = "Absolute"
  )

p_value_test_confirmation = "Y"
if(p_value_test_confirmation == "Y"){
  # Original list with 2-5 components
  Combination_Table <- unique(IRAK2_Size_Summary_ByImage_LessThan$COHORT)
  
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
    
    p_value_Result <- IRAK2_Size_Summary_ByImage_LessThan %>% 
      rowwise() %>% 
      mutate(
        TEST = COHORT %in% Combination$SHORT_LABEL
      ) %>% 
      filter(
        TEST == TRUE
      )
    
    p_value_Result_1 <- wilcox.test(
      data = p_value_Result,
      MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT ~ COHORT
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


p_value_Table_Percentage <- p_value_Table

p_value_test_confirmation = "Y"
if(p_value_test_confirmation == "Y"){
  # Original list with 2-5 components
  Combination_Table <- unique(IRAK2_Size_Summary_ByImage_GreaterThan$COHORT)
  
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
    
    p_value_Result <- IRAK2_Size_Summary_ByImage_GreaterThan %>% 
      rowwise() %>% 
      mutate(
        TEST = COHORT %in% Combination$SHORT_LABEL
      ) %>% 
      filter(
        TEST == TRUE
      )
    
    p_value_Result_1 <- wilcox.test(
      data = p_value_Result,
      MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT ~ COHORT
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


p_value_Table_Percentage <- rbind(p_value_Table_Percentage, p_value_Table)
p_value_Table_Percentage <- p_value_Table_Percentage %>% 
  mutate(
    COMPARISON_DATA = "Percentage"
  )

p_value_Table_Complete <- rbind(p_value_Table_Complete, p_value_Table_Percentage)
rm(p_value_Table_Percentage)

# Save "Less Than" data frames to an Excel file
write_xlsx(
  list(
    "p_value" = p_value_Table_Complete,
    "IRAK2_Size_Summary_LessThan" = IRAK2_Size_Summary_LessThan,
    "IRAK2_Size_Summary_ByImage_LessThan" = IRAK2_Size_Summary_ByImage_LessThan,
    "IRAK2_Size_Summary_ByCohort_LessThan" = IRAK2_Size_Summary_ByCohort_LessThan,
    "IRAK2_Size_Summary_GreaterThan" = IRAK2_Size_Summary_GreaterThan,
    "IRAK2_Size_Summary_ByImage_GreaterThan" = IRAK2_Size_Summary_ByImage_GreaterThan,
    "IRAK2_Size_Summary_ByCohort_GreaterThan" = IRAK2_Size_Summary_ByCohort_GreaterThan
  ),
  path = file.path(Plot_Directory_Save_Path_Sub_Dir, "01_NumberOfComplementaryProtein_vs_Condition_by_Cohort_DataAndpvalue.xlsx")
)

# Write csv ---------------------------------------------------------------
Source_Data_Path <- paste0(basename(Plot_Directory_Save_Path_Sub_Dir), " ",  basename(Plot_Script_Directory), ".csv")
Plot_Save_Path <- file.path(Plot_Directory_Save_Path_Sub_Dir, Source_Data_Path)
write.csv(IRAK2_Size_Summary_ByImage_GreaterThan, Plot_Save_Path)


# Percentage Plot ---------------------------------------------------------
Plot <- ggplot(
  data = IRAK2_Size_Summary_ByImage_GreaterThan,
  aes(
    fill = SHORT_LABEL
  )
) +
  geom_bar(
    data = IRAK2_Size_Summary_ByCohort_GreaterThan,
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
    data = IRAK2_Size_Summary_ByCohort_GreaterThan,
    aes(
      x = MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT,
      y = SHORT_LABEL
    ),
    fill = "#B2B2B2",
    width = 0.5,
    stat = "identity",
    color = "black",
    linewidth = 0.1
  ) +
  geom_errorbar(
    data = IRAK2_Size_Summary_ByCohort_GreaterThan,
    aes(
      y = SHORT_LABEL,
      x = MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT,
      xmin = YMIN_MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT,
      xmax = YMAX_MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT
    ),
    width = 0.1,
    size = 0.1
  ) +
  # Adding Mean Line
  geom_crossbar(
    data = IRAK2_Size_Summary_ByCohort_GreaterThan,
    aes(
      y = SHORT_LABEL,
      x = MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT,
      xmin = MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT,
      xmax = MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT
    ),
    width = 0.1, 
    size = 0.1, 
    color = "black" 
  ) +
  # Add Image Means
  geom_jitter(
    data = IRAK2_Size_Summary_ByImage_GreaterThan,
    aes(
      x = MEAN_PERCENTAGE_OF_TRACKS_IN_COHORT,
      y = SHORT_LABEL
    ),
    shape = 21,
    fill = IRAK2_Size_Summary_ByImage_GreaterThan$COLOR,
    size = 1.2,
    stroke = 0.2,
    height = 0.1  # Control vertical spread
  ) + 
  scale_x_continuous(
    limits = c(-0.2, 100),
    expand = c(0, 0),
    position = "top",
  ) +
  labs(
    x = "Fraction"
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 6, hjust = 0),
    panel.background = element_blank(), # Remove the panel background
    plot.background = element_blank(),  # Remove the plot background
    legend.background = element_blank(), # Remove the legend background
    axis.line.x.top = element_line(color = "black"), # Add x-axis line on top
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )


Plot_Save_Path_1 <- paste0("02_PercentageOfComplementaryProtein_vs_Condition_by_Cohort_with_p_value.pdf")
Plot_Save_Path <- file.path(Plot_Directory_Save_Path_Sub_Dir, Plot_Save_Path_1)


ggsave(
  Plot_Save_Path,
  plot = Plot,
  height = 40,
  width = 50,
  units = "mm"
)

# Cleanup --------------------------------------------------------
rm(
  IRAK2_Size_Summary_LessThan,
  IRAK2_Size_Summary_ByImage_LessThan,
  IRAK2_Size_Summary_ByCohort_LessThan,
  IRAK2_Size_Summary_GreaterThan,
  IRAK2_Size_Summary_ByImage_GreaterThan,
  IRAK2_Size_Summary_ByCohort_GreaterThan,
  p_value_Table_Complete,
  Plot_Save_Path_1,
  Plot_Save_Path,
  Plot,
  Plot_Directory_Save_Path_Sub_Dir,
  p_value_test_confirmation,
  p_value_Table
)

