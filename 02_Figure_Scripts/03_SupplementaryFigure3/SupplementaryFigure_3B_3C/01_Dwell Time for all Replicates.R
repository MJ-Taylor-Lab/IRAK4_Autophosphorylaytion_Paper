# Filter and summarize the data by track
Function_script <- file.path(Function_folder, "CellSummarybyImage_ForDwellTimeBins.R")
source(Function_script)
Cell_Summary_by_Image <- Cell_Summary_by_Image(Table)


# color for non illustrator plot
Cell_Summary_by_Image <- Cell_Summary_by_Image %>% 
  mutate(
    COLOR = case_when(
      SHORT_LABEL == "IRAK4-KA" ~ "#F99108",  # Assign colors based on SHORT_LABEL
      SHORT_LABEL == "IRAK4-WT" ~ "#33AFCC"
    )
  ) %>% 
  arrange(
    SHORT_LABEL
  )

# Convert SHORT_LABEL to a factor with specified levels
Cell_Summary_by_Image$SHORT_LABEL <- factor(
  Cell_Summary_by_Image$SHORT_LABEL,
  levels = c("IRAK4-KA", "IRAK4-WT")
)

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

# Summarize the data by cohort
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
    PRIMARY_PROTIEN_TRACK_COUNT_byCOHORT = max(PRIMARY_PROTIEN_TRACK_COUNT_byCOHORT),   # Number of Primary Protien tracks
    MIN_NUMBER_OF_CELLS_PER_REPLICATE = min(MIN_NUMBER_OF_CELLS_PER_REPLICATE),   # Minimum Number of cells in a replicate
    MAX_NUMBER_OF_CELLS_PER_REPLICATE = max(MAX_NUMBER_OF_CELLS_PER_REPLICATE),   # Maximum Number of cells in a replicate
    TOTAL_NUMBER_OF_CELLS_PER_COHORT = max(TOTAL_NUMBER_OF_CELLS_PER_COHORT)
  )


# Write csv ---------------------------------------------------------------
Source_Data_Path <- paste0(basename(Plot_Directory_Save_Path), " ",  basename(Plot_Script_Directory), "_Cell_Summary_by_Image.csv")
Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Source_Data_Path)
write.csv(Cell_Summary_by_Image, Plot_Save_Path)

Source_Data_Path <- paste0(basename(Plot_Directory_Save_Path), " ",  basename(Plot_Script_Directory), "_Cell_Summary_by_Cohort.csv")
Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Source_Data_Path)
write.csv(Cell_Summary_by_Cohort, Plot_Save_Path)

# Plot --------------------------------------------------------------------
# Plot histogram using ggplot2
Function_script <- file.path(Function_folder, "ComplexLifetimeBarPlot_2Cohort.R")
source(Function_script)
Plot <- Plot_ComplexLifetime(Cell_Summary_by_Image, Cell_Summary_by_Cohort)


# Add p-value and other details to the plot
Plot_Mean_SEM_pvalue <- Plot_ComplexLifetime_withMeanSEMandpvalues(
  Cell_Summary_by_Image,
  Cell_Summary_by_Cohort,
  Plot,
  p_value_Table
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


#Plot for Publication
Plot <- Plot_ComplexLifetime_forPublication(
  Cell_Summary_by_Image, 
  Cell_Summary_by_Cohort, 
  Plot
)


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
rm(
  Function_script,
  Cell_Summary_by_Image,
  Cell_Summary_by_Cohort,
  p_value_Result_1,
  Plot,
  Plot_Mean_SEM_pvalue,
  Plot_Save_Path_1,
  Plot_Save_Path,
  Plot_Directory_Save_Path
)

rm(
  Cell_Summary_by_Image,
  p_value_test_confirmation,
  p_value_Table,
  Source_Data_Path,
  Combination_Function,
  Plot_ComplexLifetime,
  Plot_ComplexLifetime_withMeanSEMandpvalues,
  Plot_ComplexLifetime_forPublication,
  Data_Colocalisation_Percentage,
  Data_Dwell_Time,
  Plot_ComplexLifetime_forPublication,
  Plot_ComplexLifetime,
  Plot_ComplexLifetime_withMeanSEMandpvalues
)

gc()
