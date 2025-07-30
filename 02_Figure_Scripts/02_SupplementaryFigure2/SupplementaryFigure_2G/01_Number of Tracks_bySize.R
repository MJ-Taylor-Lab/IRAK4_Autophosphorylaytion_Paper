# DATA PROCESSING-----------------------------------------------
Function_script <- file.path(Function_folder, "CellSummarybyImage_ForSizeofComplementaryProtein.R")
source(Function_script)

Size_Summary_ByImage <- Size_Summary_ByImage(Table, SIZE_THRESHOLD)

# Assigning the right table to the right variable from function
Size_Summary_ByCohort <- Size_Summary_ByImage$Size_Summary_ByCohort
Size_Summary_ByImage <- Size_Summary_ByImage$Size_Summary_ByImage



Size_Summary_ByImage <- Size_Summary_ByImage %>% 
  mutate(
    COLOR = case_when(
      SHORT_LABEL == "IRAK3-DD" ~ "#F99108",  # Assign colors based on SHORT_LABEL
      SHORT_LABEL == "IRAK3-WT" ~ "#33AFCC"
    )
  ) 

# Percentage Plot p value calculation---------------------------------------------------------
p_value_test_confirmation = "Y"
if(p_value_test_confirmation == "Y"){
  # Original list with 2-5 components
  Combination_Table <- unique(Size_Summary_ByImage$SHORT_LABEL)
  
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
    
    p_value_Result <- Size_Summary_ByImage %>% 
      rowwise() %>% 
      mutate(
        TEST = SHORT_LABEL %in% Combination$SHORT_LABEL
      ) %>% 
      filter(
        TEST == TRUE
      )
    
    p_value_Result_1 <- wilcox.test(
      data = p_value_Result,
      MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD ~ SHORT_LABEL
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


# Write csv ---------------------------------------------------------------
Plot_Save_Path <- file.path(Plot_Directory_Path, "00_Size_Summary_ByImage.csv")
write.csv(Size_Summary_ByImage, Plot_Save_Path)

Plot_Save_Path <- file.path(Plot_Directory_Path, "00_Size_Summary_ByCohort.csv")
write.csv(Size_Summary_ByCohort, Plot_Save_Path)

# Plot ---------------------------------------------------------
# Plot histogram using ggplot2
Function_script <- file.path(Function_folder, "SizeBarPlot_2Cohort.R")
source(Function_script)
Plot <- Plot_Size(Size_Summary_ByImage, Size_Summary_ByCohort)


# Add p-value and other details to the plot
Plot_Mean_SEM_pvalue <- Plot_Size_withMeanSEMandpvalues(
  Size_Summary_ByImage,
  Size_Summary_ByCohort,
  Plot,
  p_value_Table
)

# Save the plot as a PDF
Plot_Save_Path_1 <- "01_Percentage of Complementary Protien with size greater than threshold and p_value labels.pdf"
Plot_Save_Path <- file.path(Plot_Directory_Path, Plot_Save_Path_1)
ggsave(
  Plot_Save_Path,
  plot = Plot_Mean_SEM_pvalue,
  height = 3 * 3,
  width = 5 * 4
)


#Plot for Publication
Plot <- Plot_Size_forPublication(
  Size_Summary_ByImage, 
  Size_Summary_ByCohort, 
  Plot
)


Plot_Save_Path_1 <- "02_Percentage of Complementary Protien with size greater than threshold.pdf"
Plot_Save_Path <- file.path(Plot_Directory_Path, Plot_Save_Path_1)
ggsave(
  Plot_Save_Path,
  plot = Plot,
  height = 40,
  width = 40,
  units = "mm"
)

# Cleanup --------------------------------------------------------
rm(
  p_value_Table,
  p_value_test_confirmation,
  Plot,
  Plot_Mean_SEM_pvalue,
  Size_Summary_ByCohort,
  Size_Summary_ByImage,
  Plot_Save_Path_1,
  Plot_Save_Path,
  Combination_Function,
  Plot_Size,
  Plot_Size_forPublication,
  Plot_Size_withMeanSEMandpvalues,
  Function_script
)

