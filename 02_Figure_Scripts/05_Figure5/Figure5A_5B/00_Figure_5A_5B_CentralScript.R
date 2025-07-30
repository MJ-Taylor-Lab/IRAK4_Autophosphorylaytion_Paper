library(pacman)
pacman::p_load(dplyr, tidyr, data.table, ggplot2)

# Loading data
Table1 <- fread("/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/43_3xKO+hMyD88-3xFLAG-IRES-TagBFP2+hIRAK4WT-GFP+hIRAK1WT-mScarlet_Analysis.csv.gz")
Table2 <- fread("/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/44_3xKO+hMyD88-3xFLAG-IRES-TagBFP2+hIRAK4D329A-GFP+hIRAK1WT-mScarlet_Analysis.csv.gz")
Table3 <- fread("/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/45_3xKO+hMyD88-3xFLAG-IRES-TagBFP2+hIRAK4D329A-Phosphomimetic-GFP+hIRAK1WT-mScarlet_Analysis.csv.gz")

# Setting directory paths
Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots"
Plot_Script_Directory <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/02_Figure_Scripts/05_Figure5/Figure5A_5B"
# Sourcing Function folder
Function_folder <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/02_Figure_Scripts/00_Functions"

# Merging tables
Table <- rbind(
  Table1,
  Table2,
  Table3
)

# Removing unnecessary objects
rm(
  Table1,
  Table2,
  Table3
)

# Adjusting ORDER_NUMBER and filtering data
Table <- Table %>% 
  # filter(
  #   IMAGE != "20250401 plate2_well3G_40nM_cl621_hIRAK4D329A_hIRAK1WT_"
  # ) %>% 
  mutate(
    SHORT_LABEL = case_when(
      SHORT_LABEL == "hIRAK4-WT" ~ "WT",
      SHORT_LABEL == "hIRAK4-D329A" ~ "D329A",
      SHORT_LABEL == "hIRAK4-D329A-Phosphomimetic" ~ "D329A/\nT342E/\nT345E/\nS346E"
    ),
    PROTEIN_GENERIC = case_when(
      PROTEIN == "hIRAK4WT" ~ "IRAK4",
      PROTEIN == "hIRAK4-D329A" ~ "IRAK4",
      PROTEIN == "hIRAK4-D329A-Phosphomimetic" ~ "IRAK4",
      TRUE ~ PROTEIN
    )
  ) %>% 
  as.data.table()

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Save_Path, "05_Figure 5")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Path, "Figure_5A_5B")
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

# Dwell Time hIRAK4-WT vs hIRAK4-D329A vs hIRAK4-D329A-PM ------------------------------------------------
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

rm(
  ScriptList,
  x
)


# Colocalisation IRAK4 WT vs IRAK4 DD----------------------------------------------------------
Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "02_Colocalisation")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

ScriptList <- c(
  "02_Colocalisation Percentage with Dwell Time greater than 30s.R" # Max Normalized Intensity MyD88 Table Creation.R
)

LOWER_LIMIT = -1 # The Plot x-axis Lower limit
UPPER_LIMIT = 71 # The Plot a-axis Upper limit
LOWER_LIMIT_AXIS_TICKS = 0 # Axis Lower Limit
X_AXIS_INTERVAL = 30 # Axis tick mark Interval 

# Command
for(x in 1:length(ScriptList)){
  tryCatch({
    print(paste("::::::::::::::::::::", 1, "::::::::::::::::::::"))
    setwd(Plot_Script_Directory)
    source(ScriptList[1], local = T)
  }, error = function(e) {print(paste("Error loading", ScriptList[1]))})
}

rm(
  LOWER_LIMIT,
  UPPER_LIMIT,
  LOWER_LIMIT_AXIS_TICKS,
  X_AXIS_INTERVAL,
  ScriptList,
  x
)

# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()
