library(pacman)
# `p_load` checks if a package is installed and loads it; if not, it installs and then loads it.
# This makes the script more portable as it handles dependencies automatically.
pacman::p_load(dplyr, tidyr, data.table, ggplot2)
# Dependencies:
# - dplyr: For data manipulation (group_by, mutate, filter, etc.)
# - tidyr: Potentially for reshaping data (though not explicitly used in this snippet, often paired with dplyr)
# - data.table: For efficient data loading (fread) and potentially faster operations on large datasets.
# - ggplot2: For creating high-quality data visualizations.

# Loading data
# especially faster than `read.csv` or `read_csv`. The `.gz` extension indicates compressed files,
# which `fread` handles directly, a convenient feature.
Table1 <- fread("/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/01_MyD88-GFP+IRAK4-KO+IRAK4WT-mScarlet_Analysis.csv.gz")
Table2 <- fread("/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/02_MyD88-GFP+IRAK4-KO+IRAK4DD-mScarlet_Analysis.csv.gz")

# Setting directory paths
Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots"
Plot_Script_Directory <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/02_Figure_Scripts/01_Figure1/Figure_1C_1D"
Function_folder <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/02_Figure_Scripts/00_Functions"

# Merging tables
# the same columns and structure, which is typical for experimental replicates or conditions.
Table <- rbind(
  Table1,
  Table2
)

# Removing unnecessary objects
# This is a good practice when dealing with large datasets to manage memory, especially in R.
rm(
  Table1,
  Table2
)

# It's a common data cleaning/standardization step.
Table <- Table %>%
  mutate(
    SHORT_LABEL = case_when(
      SHORT_LABEL == "IRAK4 WT" ~ "WT", # Renames "IRAK4 WT" to "WT"
      SHORT_LABEL == "IRAK4 DD" ~ "IRAK4DD", # Renames "IRAK4 DD" to "IRAK4DD" (no actual change here, perhaps it was intended for something else or for consistency if other labels existed)
      TRUE ~ SHORT_LABEL # Keeps all other SHORT_LABEL values as they are.
    )
  ) %>%
  as.data.table()

# Creating Folder to save Images
# This ensures that plot saving operations don't fail due to missing directories.
Plot_Directory_Path <- file.path(Plot_Directory_Save_Path, "01_Figure 1")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}
# Creates a subfolder for "Figure 1".

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Path, "Figure_1C_1D")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}
# Creates another nested subfolder specific to "Figure 1C_1D".
# Clarity: Using `file.path` is good for cross-platform compatibility.

# Primary Protein Filtering
Table <- Table %>%
  filter(
    PROTEIN == "MyD88", # Only include rows where the protein is MyD88. This suggests the raw data contains multiple protein types.
    MAX_NORMALIZED_INTENSITY >= 1.5, # Filters for tracks where the maximum intensity reached a certain threshold.
    NORMALIZED_INTENSITY >= 0.75 # Filters for individual frames where the intensity meets a lower threshold.
  )

COLOCALISATION_THREHOLD = 1.5
# is used in the `Cell_Summary_by_Track_for_ColocalisationPercentage` function (as seen in the previous review).

# Dwell Time IRAK4 WT vs IRAK4 DD------------------------------------------------
Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "01_Dwell Time")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}
# Creates a dedicated subfolder for dwell time plots.

ScriptList <- c(
  "01_Dwell Time for all Replicates.R" # Dwell Time
)

# Command
for(x in 1:length(ScriptList)){
  tryCatch({
    print(paste("::::::::::::::::::::", x, "::::::::::::::::::::")) # Prints progress.
    setwd(Plot_Script_Directory) # Sets the working directory to where the scripts are located. This is important for `source()` to find the file.
    source(ScriptList[x], local = T) # Executes the R script. `local = T` ensures variables defined in the sourced script
  }, error = function(e) {print(paste("Error loading", ScriptList[x], ":", e$message))})
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
# Creates a dedicated subfolder for colocalization plots.

ScriptList <- c(
  "02_Colocalisation Percentage with Dwell Time greater than 30s.R" # Max Normalized Intensity MyD88 Table Creation.R
)

LOWER_LIMIT = 0 # The Plot x-axis Lower limit
UPPER_LIMIT = 46 # The Plot a-axis Upper limit
LOWER_LIMIT_AXIS_TICKS = 0 # Axis Lower Limit
X_AXIS_INTERVAL = 20 # Axis tick mark Interval
# - Global variables: Same as `COLOCALISATION_THREHOLD`, better to pass these as arguments to plotting functions.

# Command
for(x in 1:length(ScriptList)){
  tryCatch({
    print(paste("::::::::::::::::::::", x, "::::::::::::::::::::"))
    setwd(Plot_Script_Directory) # Again, sets working directory.
    source(ScriptList[x], local = T) # Executes the script. Same consideration about `local = T`.
  }, error = function(e) {print(paste("Error loading", ScriptList[x], ":", e$message))})
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
rm(list = ls()) # Removes all objects from the current R environment.
gc() # Forces garbage collection, freeing up memory.