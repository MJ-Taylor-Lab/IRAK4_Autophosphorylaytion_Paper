library(pacman)
pacman::p_load(dplyr, tidyr, data.table, ggplot2)

# Loading data
Table1 <- fread("/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/10_3xKO+MyD88-3xFLAG-IRES-Puro+mIRAK4WT-GFP+mIRAK2WT-mScarlet_Analysis.csv.gz")
Table2 <- fread("/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/11_3xKO+MyD88-3xFLAG-IRES-Puro+mIRAK4WT-GFP+mIRAK2DD-mScarlet_Analysis.csv.gz")

# Setting directory paths
Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots"
Plot_Script_Directory <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/02_Figure_Scripts/02_Figure2/Figure_2J_2K"
Plot_Script_Directory_Supplementary <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/02_Figure_Scripts/02_SupplementaryFigure2/SupplementaryFigure_2F"

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
      SHORT_LABEL == "IRAK4 WT + IRAK2 WT" ~ "IRAK2-WT",
      SHORT_LABEL == "IRAK4 WT + IRAK2 DD" ~ "IRAK2-DD"
    ),
    PROTEIN_GENERIC = case_when(
      PROTEIN == "IRAK2WT" ~ "IRAK2",
      PROTEIN == "IRAK2DD" ~ "IRAK2",
      PROTEIN == "IRAK4WT" ~ "IRAK4",
      TRUE ~ PROTEIN
    )
  ) %>% 
  as.data.table()

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Save_Path, "02_Figure 2")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Path, "Figure_2J_2K")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}


# Primary Protein
Table <- Table %>% 
  filter(
    PROTEIN_GENERIC == "IRAK4",                   # Only include rows where the protein is MyD88
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
    print(paste("::::::::::::::::::::", 1, "::::::::::::::::::::"))
    setwd(Plot_Script_Directory)
    source(ScriptList[1], local = T)
  }, error = function(e) {print(paste("Error loading", ScriptList[1]))})
}

rm(ScriptList)


# Colocalisation IRAK4 WT vs IRAK4 DD----------------------------------------------------------
Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "02_Colocalisation")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

ScriptList <- c(
  "02_Colocalisation Percentage with Dwell Time greater than 30s.R" # Max Normalized Intensity MyD88 Table Creation.R
)

LOWER_LIMIT = -1 # The Plot x-axis Lower limit
UPPER_LIMIT = 46 # The Plot a-axis Upper limit
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
  ScriptList
)

# Size of Complementary Protien -------------------------------------------
Plot_Directory_Path <- dirname(Plot_Directory_Path)
Plot_Directory_Path <- dirname(Plot_Directory_Path)

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Path, "02_Supplementary Figure 2")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Path, "Supplementary_Figure_2F")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}

ScriptList <- c(
  "01_Number of Tracks_bySize.R" # Dwell Time
)

SIZE_THRESHOLD = 4

# Command
for(x in 1:length(ScriptList)){
  tryCatch({
    print(paste("::::::::::::::::::::", 1, "::::::::::::::::::::"))
    setwd(Plot_Script_Directory_Supplementary)
    source(ScriptList[1], local = T)
  }, error = function(e) {print(paste("Error loading", ScriptList[1]))})
}

rm(
  ScriptList,
  SIZE_THRESHOLD,
  x
)

# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()
