library(pacman)
pacman::p_load(dplyr, tidyr, data.table, ggplot2)

# Loading data
Table1 <- fread("/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/03_3xKO+MyD88-3xFLAG-IRES-Puro+mIRAK4WT-GFP+mIRAK1WT-mScarlet_Analysis.csv.gz")
Table2 <- fread("/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/04_3xKO+MyD88-3xFLAG-IRES-Puro+mIRAK4DD-GFP+mIRAK1WT-mScarlet_Analysis.csv.gz")

# Setting directory paths
Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots"
Plot_Script_Directory <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/02_Figure_Scripts/01_Figure1/Figure_1F_1G"
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
    PROTEIN = case_when(
      PROTEIN == "IRAK4WT" | PROTEIN == "IRAK4DD" ~ "IRAK4",
      PROTEIN == "IRAK1WT" ~ PROTEIN 
    ),
    SHORT_LABEL = case_when(
      SHORT_LABEL == "IRAK4 WT + IRAK1 WT" ~ "WT",
      SHORT_LABEL == "IRAK4 DD + IRAK1 WT" ~ "IRAK4DD\n+\nIRAK1WT",
      TRUE ~ SHORT_LABEL
    )
  ) %>% 
  as.data.table()

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Save_Path, "01_Figure 1")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Path, "Figure_1F_1G")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}

# Primary Protein
Table <- Table %>% 
  filter(
    PROTEIN == "IRAK4",                   # Only include rows where the protein is MyD88
    MAX_NORMALIZED_INTENSITY >= 1.5,      # Filter rows with maximum normalized intensity >= 1.5
    NORMALIZED_INTENSITY >= 0.75          # Filter rows with normalized intensity >= 0.75
  )


COLOCALISATION_THREHOLD = 1.5          # Setting Threshold to consider Protien of Interest is Colocalised with Complementary Protien


# Dwell Time IRAK4 WT vs IRAK4 DD------------------------------------------------
Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "01_Dwell Time")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

ScriptList <- c(
  "01_Dwell Time for all Replicates.R" # Dwell Time
)

# Command
for(x in 1:length(ScriptList)){
  tryCatch({
    print(paste("::::::::::::::::::::", x, "::::::::::::::::::::"))
    setwd(Plot_Script_Directory)
    source(ScriptList[x], local = T)
  }, error = function(e) {print(paste("Error loading", ScriptList[x]))})
}

rm(
  x,
  ScriptList
)


# Colocalisation IRAK4 WT vs IRAK4 DD----------------------------------------------------------
Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "02_Colocalisation")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

ScriptList <- c(
  "02_Colocalisation Percentage with Dwell Time greater than 30s.R" # Max Normalized Intensity MyD88 Table Creation.R
)

LOWER_LIMIT = -0.5 # The Plot x-axis Lower limit
UPPER_LIMIT = 36 # The Plot a-axis Upper limit
LOWER_LIMIT_AXIS_TICKS = 0 # Axis Lower Limit
X_AXIS_INTERVAL = 15 # Axis tick mark Interval 

# Command
for(x in 1:length(ScriptList)){
  tryCatch({
    print(paste("::::::::::::::::::::", x, "::::::::::::::::::::"))
    setwd(Plot_Script_Directory)
    source(ScriptList[x], local = T)
  }, error = function(e) {print(paste("Error loading", ScriptList[x]))})
}

rm(
  LOWER_LIMIT,
  UPPER_LIMIT,
  LOWER_LIMIT_AXIS_TICKS,
  X_AXIS_INTERVAL,
  x,
  ScriptList
)


# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()
