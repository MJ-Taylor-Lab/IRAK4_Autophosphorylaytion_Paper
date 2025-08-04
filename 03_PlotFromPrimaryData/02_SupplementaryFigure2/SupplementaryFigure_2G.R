library(pacman)
pacman::p_load(readxl, dplyr, data.table, ggplot2)

# Setting directory paths
Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots_short"

Primary_Data <- "/Users/u_niranjan/Desktop/Cell_Figures/00_New_Figures/Primary_Data.xlsx"

Function_folder <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/03_PlotFromPrimaryData/Functions"

# Opening Data for relevant plot ----------------------------------------
Size_Summary_ByImage <- read_excel(Primary_Data, sheet = "Supplementary Figure 2G")

# Creating folder to save data --------------------------------------------
# Creating Folder to save Images
Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "02_SupplementaryFigure 2")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

# Creating Folder to save Images
Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "Supplementary Figure 2G")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}


# Convert SHORT_LABEL to a factor with specified levels
Size_Summary_ByImage$SHORT_LABEL <- factor(
  Size_Summary_ByImage$SHORT_LABEL,
  levels = c("IRAK3-DD", "IRAK3-WT")
)

Size_Summary_ByImage <- Size_Summary_ByImage %>% 
  mutate(
    COLOR = case_when(
      SHORT_LABEL == "IRAK3-WT" ~ "#33AFCC",
      SHORT_LABEL == "IRAK3-DD" ~ "#F99108"
    )
  ) 


# Summary statistics by cohort
Size_Summary_ByCohort <- Size_Summary_ByImage %>% 
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
    SHORT_LABEL
  ) %>% 
  summarise(
    TOTAL_NUMBER_OF_IMAGES = mean(TOTAL_NUMBER_OF_IMAGES),
    #Finding SD, SEM, MEAN, YMAX and YMIN
    SD = sd(MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD),
    SEM = SD/TOTAL_NUMBER_OF_IMAGES,
    MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD = mean(MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD),
    Y_MAX = MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD + SEM,
    Y_MIN = case_when(
      MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD - SEM < 0 ~ 0,
      TRUE ~ MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD - SEM
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
Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Plot_Save_Path_1)
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
