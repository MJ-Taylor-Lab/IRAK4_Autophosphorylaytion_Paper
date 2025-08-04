library(pacman)
pacman::p_load(readxl, dplyr, data.table, ggplot2)

# Setting directory paths
Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots_short"

Primary_Data <- "/Users/u_niranjan/Desktop/Cell_Figures/00_New_Figures/Primary_Data.xlsx"

Function_folder <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/03_PlotFromPrimaryData/Functions"

#Variables needed for Plot
LOWER_LIMIT = -0.1 # The Plot x-axis Lower limit
UPPER_LIMIT = 151 # The Plot a-axis Upper limit
LOWER_LIMIT_AXIS_TICKS = 0 # Axis Lower Limit
X_AXIS_INTERVAL = 75 # Axis tick mark Interval 

# Opening Data for relevant plot ----------------------------------------
Cell_Summary_by_Cell <- read_excel(Primary_Data, sheet = "Figure 4J")

# Creating folder to save data --------------------------------------------
# Creating Folder to save Images
Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "04_Figure 4")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

# Creating Folder to save Images
Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "Figure_4J")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

Cell_Summary_by_Cell <- Cell_Summary_by_Cell %>% 
  mutate(
    SHORT_LABEL = case_when(
      SHORT_LABEL == "IRAK4 Kinase\r\nInhibitor 20 ¬µM" ~ "IRAK4 Kinase\nInhibitor\n20 µM",
      TRUE ~ SHORT_LABEL
    )
  )

# Convert SHORT_LABEL to a factor with specified levels
Cell_Summary_by_Cell$SHORT_LABEL <- factor(
  Cell_Summary_by_Cell$SHORT_LABEL,
  levels = c("IRAK4 Kinase\nInhibitor\n20 µM", "DMSO")
)


# Data Processing ---------------------------------------------------------
Function_script <- file.path(Function_folder, "CellSummary_ForKinaseDomainOnly_ViolinPlot.R")
source(Function_script)
Cell_Summary_by_Cell <- Number_of_Association_Events_Data_Processing(Cell_Summary_by_Cell)

# Assigning the right table to the right variable from function
Cell_Summary_by_Image <- Cell_Summary_by_Cell$Cell_Summary_by_Image
Cell_Summary_by_Cohort <- Cell_Summary_by_Cell$Cell_Summary_by_Cohort
Cell_Number_Label <- Cell_Summary_by_Cell$Cell_Number_Label
p_value_Result <- Cell_Summary_by_Cell$p_value_Result_1
Cell_Summary_by_Cell <- Cell_Summary_by_Cell$Cell_Summary_by_Cell

Cell_Summary_by_Image <- Cell_Summary_by_Image %>% 
  mutate(
    COLOR = case_when(
      SHORT_LABEL == "DMSO" ~ "#FED09E",
      SHORT_LABEL == "IRAK4 Kinase\nInhibitor\n20 µM" ~ "#CC7B16"
    )
  )

# GGPLOT VIOLIN -----------------------------------------------------------
# Generate a Violin Baseplot
Function_script <- file.path(Function_folder, "KinaseDomainOnly_Plot.R")
source(Function_script)
Plot <- Kinase_Domain_Only_Violin_Fx(
  Cell_Summary_by_Cell,
  Cell_Summary_by_Cohort,
  Cell_Summary_by_Image
)


Plot_pvalue <- Kinase_Domain_Only_Violin_pvalue_Fx(
  Plot, 
  Cell_Number_Label,
  p_value_Result
)

# Save the reference plot with P-value as a PDF
Plot_Save_Path <- file.path(Plot_Directory_Save_Path, "01_Number of Association Events Violin Stats.pdf")
ggsave(
  Plot_Save_Path,
  plot = Plot_pvalue,
  height = 3 * 3,
  width = 5 * 4
)

Plot_Save_Path <- file.path(Plot_Directory_Save_Path, "02_Number of Association Events Violin.pdf")
ggsave(
  Plot_Save_Path, 
  plot = Plot, 
  height = 50, 
  width = 50, 
  units = "mm"
)

# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()