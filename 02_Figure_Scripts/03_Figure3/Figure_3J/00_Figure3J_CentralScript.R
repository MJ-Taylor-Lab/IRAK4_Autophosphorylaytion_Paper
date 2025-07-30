library(pacman)
pacman::p_load(dplyr, tidyr, data.table, ggplot2, drc)

# Loading data
Table1 <- fread("/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/23_MyD88-GFP_IRAK1-mScarlet_DMSO_Analysis.csv.gz")
Table2 <- fread("/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/27_MyD88-GFP_IRAK1-mScarlet_Zimlovisertib20uM_Analysis.csv.gz")

# Setting directory paths
Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots"
Plot_Script_Directory <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/02_Figure_Scripts/03_Figure3/Figure_3J"
# Sourcing Function folder
Function_folder <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/02_Figure_Scripts/00_Functions"

# Merging tables
Table <- rbind(
  Table1,
  Table2
)

# Removing unnecessary objects
rm(
  Table1,
  Table2
)

# Adjusting ORDER_NUMBER and filtering data
Table <- Table %>% 
  mutate(
    SHORT_LABEL = case_when(
      SHORT_LABEL == "DMSO" ~ "DMSO",
      SHORT_LABEL == "Kinase Inhibitor 20 uM" ~ "KI"
    )
  ) %>% 
  as.data.table()

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Save_Path, "03_Figure 3")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Path, "Figure_3J")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}


# Primary Protein
Table <- Table %>% 
  filter(
    PROTEIN == "MyD88",                   # Only include rows where the protein is MyD88
    MAX_NORMALIZED_INTENSITY >= 1.5,      # Filter rows with maximum normalized intensity >= 1.5
    NORMALIZED_INTENSITY >= 0.75          # Filter rows with normalized intensity >= 0.75
  )


COLOCALISATION_THREHOLD = 1.5        # Setting Threshold to consider Protien of Interest is Colocalised with Complementary Protien

# Coloclaisation Temp -----------------------------------------------------
Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "01_Colocalisation")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

ScriptList <- c(
  "01_Colocalisation Percentage with Dwell Time greater than 30s means DMSO vs 20uM.R" # Max Normalized Intensity MyD88 Table Creation.R
)


LOWER_LIMIT = 0 # The Plot y-axis Lower limit
UPPER_LIMIT = 40 # The Plot y-axis Upper limit
LOWER_LIMIT_AXIS_TICKS = 0 # Axis Lower Limit
X_AXIS_INTERVAL = 20 # Axis tick mark Interval 

# Command
tryCatch({
  print(paste("::::::::::::::::::::", 1, "::::::::::::::::::::"))
  setwd(Plot_Script_Directory)
  source(ScriptList[1], local = T)
}, error = function(e) {print(paste("Error loading", ScriptList[1]))})


rm(
  LOWER_LIMIT,
  UPPER_LIMIT,
  LOWER_LIMIT_AXIS_TICKS,
  X_AXIS_INTERVAL,
  ScriptList
)


# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()
