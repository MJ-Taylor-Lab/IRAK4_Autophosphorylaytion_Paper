library(pacman)
pacman::p_load(readxl, dplyr, data.table, ggplot2, lemon, scales)

# Setting directory paths
Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots_short"

Primary_Data <- "/Users/u_niranjan/Desktop/Cell_Figures/00_New_Figures/Primary_Data.xlsx"

Function_folder <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/03_PlotFromPrimaryData/Functions"

#Variables needed for Plot
LOWER_LIMIT = 0.1           #Axis Lower Limit
UPPER_LIMIT = 1.6       #Axis Upper Limit
AXIS_BREAK_SEQ = 0.6        #X axis tick marks
FACET_ROW_NUMBERS = 2     #Number of facet rows
X_LABEL = "Lifetime (s)"
Y_LABEL = "# of Events"

# Opening Data for relevant plot ----------------------------------------
Cell_Summary_by_Track <- read_excel(Primary_Data, sheet = "Figure 4H")

# Creating folder to save data --------------------------------------------
# Creating Folder to save Images
Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "04_Figure 4")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

# Creating Folder to save Images
Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "Figure_4H")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

Cell_Summary_by_Track <- Cell_Summary_by_Track %>% 
  mutate(
    SHORT_LABEL = case_when(
      SHORT_LABEL == "IRAK4 Kinase\r\nInhibitor 20 ¬µM" ~ "IRAK4 Kinase\nInhibitor\n20 µM",
      TRUE ~ SHORT_LABEL
    )
  )

# Data Processing ---------------------------------------------------------
Function_script <- file.path(Function_folder, "CellSummary_ForKinaseDomainOnly_ViolinPlot.R")
source(Function_script)
Cell_Summary_by_Track <- Lifetime_Events_Data_Processing(Cell_Summary_by_Track)

# Assigning the right table to the right variable from function
Cell_Summary_by_Bin <- Cell_Summary_by_Track$Cell_Summary_by_Bin
Label <- Cell_Summary_by_Track$Label
Cell_Summary_by_Track <- Cell_Summary_by_Track$Cell_Summary_by_Track

# Convert SHORT_LABEL to a factor with specified levels
Cell_Summary_by_Bin$SHORT_LABEL <- factor(
  Cell_Summary_by_Bin$SHORT_LABEL,
  levels = c("DMSO", "IRAK4 Kinase\nInhibitor\n20 µM")
)

# GGPLOT VIOLIN -----------------------------------------------------------
# Generate a Violin Baseplot
Function_script <- file.path(Function_folder, "KinaseDomainOnly_Plot.R")
source(Function_script)
Plot <- Kinase_Domain_Only_Lifetime_Fx(Cell_Summary_by_Bin)


PLOT_MEAN_POSITIION = max(Cell_Summary_by_Bin$BIN_COUNT)/2 #Place the position of means at around half the y axis limit
Plot_Mean <- Kinase_Domain_Only_Lifetime_Mean_Fx(
  Plot, 
  PLOT_MEAN_POSITIION,
  Label
)

# Save the reference plot with P-value as a PDF
Plot_Save_Path <- file.path(Plot_Directory_Save_Path, "01_Lifetime Events Violin Stats.pdf")
ggsave(
  Plot_Save_Path,
  plot = Plot_Mean,
  height = 3 * 3,
  width = 5 * 4
)

Plot <- Kinase_Domain_Only_Lifetime_Publication_Fx(Plot)

Plot_Save_Path <- file.path(Plot_Directory_Save_Path, "02_Lifetime Events Violin.pdf")
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