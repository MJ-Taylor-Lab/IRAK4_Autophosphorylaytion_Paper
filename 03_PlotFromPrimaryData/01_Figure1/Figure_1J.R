library(pacman)
pacman::p_load(readxl, dplyr, data.table, ggplot2)

# Setting directory paths
Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots_short"

Primary_Data <- "/Users/u_niranjan/Desktop/Cell_Figures/00_New_Figures/Primary_Data.xlsx"

Function_folder <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/03_PlotFromPrimaryData/Functions"

#Variables needed for Plot
LOWER_LIMIT = -0.5 # The Plot x-axis Lower limit
UPPER_LIMIT = 36 # The Plot a-axis Upper limit
LOWER_LIMIT_AXIS_TICKS = 0 # Axis Lower Limit
X_AXIS_INTERVAL = 15 # Axis tick mark Interval 

# Opening Data for relevant plot ----------------------------------------
Colocalisation_Percentage_byCell <- read_excel(Primary_Data, sheet = "Figure 1J")

# Creating folder to save data --------------------------------------------
# Creating Folder to save Images
Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "01_Figure 1")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

# Creating Folder to save Images
Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "Figure_1J")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

# Convert SHORT_LABEL to a factor with specified levels
Colocalisation_Percentage_byCell$SHORT_LABEL <- factor(
  Colocalisation_Percentage_byCell$SHORT_LABEL,
  levels = c("IRAK4WT\r\n+\r\nIRAK1DD", "WT")
)

Function_script <- file.path(Function_folder, "CellSummary_ForColocalisation.R")
source(Function_script)
Colocalisation_Percentage_byCell <- Colocalisation_Percentage_byCellFx(Colocalisation_Percentage_byCell)

# Assigning the right table to the right variable from function
Colocalisation_Percentage_byImage <- Colocalisation_Percentage_byCell$Colocalisation_Percentage_byImage
Colocalisation_Percentage_byCOHORT <- Colocalisation_Percentage_byCell$Colocalisation_Percentage_byCOHORT
Colocalisation_Percentage_byCell <- Colocalisation_Percentage_byCell$Colocalisation_Percentage_byCell

Colocalisation_Percentage_byImage <- Colocalisation_Percentage_byImage %>% 
  mutate(
    JITTER_COLOR = case_when(
      SHORT_LABEL == "WT" ~ "#33AFCC",
      SHORT_LABEL == "IRAK4WT\r\n+\r\nIRAK1DD" ~ "#F99108"
    )
  )

# t-test ------------------------------------------------------------------
p_value_test_confirmation = "Y"
if(p_value_test_confirmation == "Y"){
  # Original list with 2-5 components
  Combination_Table <- unique(Colocalisation_Percentage_byImage$SHORT_LABEL)
  
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
    
    p_value_Result <- Colocalisation_Percentage_byImage %>% 
      rowwise() %>% 
      mutate(
        TEST = SHORT_LABEL %in% Combination$SHORT_LABEL
      ) %>% 
      filter(
        TEST == TRUE
      )
    
    p_value_Result_1 <- wilcox.test(
      data = p_value_Result,
      MEAN_COLOCLIZED_SPOT_TEST ~ SHORT_LABEL
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

# GGPLOT VIOLIN -----------------------------------------------------------
# Generate a Violin Baseplot
Function_script <- file.path(Function_folder, "ColocalisationPercentage_3Cohort.R")
source(Function_script)
Plot <- Plot_ColocalisationPercentage(
  Colocalisation_Percentage_byCell, 
  Colocalisation_Percentage_byImage, 
  Colocalisation_Percentage_byCOHORT
)


Plot_pvalue <- Plot_ColocalisationPercentage_pvalues(
  Plot, 
  p_value_Table
)

# Save the reference plot with P-value as a PDF
Plot_Save_Path <- file.path(Plot_Directory_Save_Path, "01_Colocalisation Percentage Violin Stats.pdf")
ggsave(
  Plot_Save_Path,
  plot = Plot_pvalue,
  height = 3 * 3,
  width = 5 * 4
)



### Plot without stats and axis labels
Plot <- Plot_ColocalisationPercentage_forPublication(Plot)

Plot_Save_Path <- file.path(Plot_Directory_Save_Path, "02_Colocalisation Percentage Violin.pdf")
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