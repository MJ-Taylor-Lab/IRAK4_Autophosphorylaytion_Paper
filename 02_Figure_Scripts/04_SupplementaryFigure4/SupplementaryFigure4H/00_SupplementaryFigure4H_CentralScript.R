library(pacman)
pacman::p_load(dplyr, tidyr, data.table, ggplot2)

# Loading data
Table1 <- fread("/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/40_MyD88-GFP+IRAK4-KO+IRAK4WT-mScarlet_Analysis.csv.gz")
Table2 <- fread("/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/41_MyD88-GFP+IRAK4-KO+IRAK4KinaseDomain-mScarlet_Analysis.csv.gz")
Table3 <- fread("/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/42_MyD88-GFP+IRAK4-KO+IRAK4R12C-mScarlet_Analysis.csv.gz")

# Setting directory paths
Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots"
Plot_Script_Directory <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/02_Figure_Scripts/04_SupplementaryFigure4/SupplementaryFigure4H"
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
  mutate(
    SHORT_LABEL = case_when(
      SHORT_LABEL == "IRAK4-WT" ~ "WT",
      SHORT_LABEL == "IRAK4-KinaseDomain" ~ "KD",
      SHORT_LABEL == "IRAK4-R12C" ~ "R12C"
    ),
    PROTEIN_GENERIC = case_when(
      PROTEIN == "MyD88" ~ "MyD88",
      TRUE ~ "IRAK4"
    )
  ) %>% 
  as.data.table()

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Save_Path, "04_SupplementaryFigure4")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Path, "SupplementaryFigure4H")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}

# Primary Protein
Table <- Table %>% 
  filter(
    PROTEIN_GENERIC == "MyD88",                   # Only include rows where the protein is MyD88
    MAX_NORMALIZED_INTENSITY >= 1.5,      # Filter rows with maximum normalized intensity >= 1.5
    NORMALIZED_INTENSITY >= 0.75          # Filter rows with normalized intensity >= 0.75
  )


COLOCALISATION_THREHOLD = 1.5          # Setting Threshold to consider Protien of Interest is Colocalised with Complementary Protien


# Colocalisation IRAK4 WT vs IRAK4 R12C vs IRAK4 Kinase Domain Only----------------------------------------------------------
ScriptList <- c(
  "01_Colocalisation Percentage with Dwell Time greater than 30s.R" # Max Normalized Intensity MyD88 Table Creation.R
)

LOWER_LIMIT = -1 # The Plot x-axis Lower limit
UPPER_LIMIT = 12 # The Plot a-axis Upper limit
LOWER_LIMIT_AXIS_TICKS = 0 # Axis Lower Limit
X_AXIS_INTERVAL = 5 # Axis tick mark Interval 

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
  ScriptList
)

# Clean up ----------------------------------------------------------------
rm(list = ls())
gc()
