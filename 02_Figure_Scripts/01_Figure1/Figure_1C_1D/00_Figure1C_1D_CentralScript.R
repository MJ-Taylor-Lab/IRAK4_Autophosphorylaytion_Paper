library(pacman)
# Comment: `pacman` is an excellent package for managing other R packages.
# `p_load` checks if a package is installed and loads it; if not, it installs and then loads it.
# This makes the script more portable as it handles dependencies automatically.
pacman::p_load(dplyr, tidyr, data.table, ggplot2)
# Dependencies:
# - dplyr: For data manipulation (group_by, mutate, filter, etc.)
# - tidyr: Potentially for reshaping data (though not explicitly used in this snippet, often paired with dplyr)
# - data.table: For efficient data loading (fread) and potentially faster operations on large datasets.
# - ggplot2: For creating high-quality data visualizations.

# Loading data
# Comment: `fread` from the `data.table` package is highly optimized for reading large delimited files,
# especially faster than `read.csv` or `read_csv`. The `.gz` extension indicates compressed files,
# which `fread` handles directly, a convenient feature.
Table1 <- fread("/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/01_MyD88-GFP+IRAK4-KO+IRAK4WT-mScarlet_Analysis.csv.gz")
Table2 <- fread("/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/02_MyD88-GFP+IRAK4-KO+IRAK4DD-mScarlet_Analysis.csv.gz")
# Critiques:
# - Hardcoded paths: The paths are absolute and specific to the user's local machine (`/Users/u_niranjan/Desktop/...`).
#   This makes the script non-portable. For sharing or running on different machines,
#   it's best practice to use relative paths, environment variables, or a configuration file.

# Setting directory paths
# Comment: Defines various directory paths where plots will be saved and where other R scripts/functions are located.
Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots"
Plot_Script_Directory <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/02_Figure_Scripts/01_Figure1/Figure_1C_1D"
Function_folder <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/02_Figure_Scripts/00_Functions"

# Merging tables
# Comment: Combines the two loaded data tables row-wise. This assumes both tables have
# the same columns and structure, which is typical for experimental replicates or conditions.
Table <- rbind(
  Table1,
  Table2
)

# Removing unnecessary objects
# Comment: Frees up memory by removing the individual table objects after they have been combined.
# This is a good practice when dealing with large datasets to manage memory, especially in R.
rm(
  Table1,
  Table2
)

# Comment: This block renames specific values within the 'SHORT_LABEL' column.
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
# Comment: This section checks if the specified output directories exist and creates them if not.
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
# Comment: Filters the `Table` to include only specific data points relevant for the analysis.
Table <- Table %>%
  filter(
    PROTEIN == "MyD88", # Only include rows where the protein is MyD88. This suggests the raw data contains multiple protein types.
    MAX_NORMALIZED_INTENSITY >= 1.5, # Filters for tracks where the maximum intensity reached a certain threshold.
    NORMALIZED_INTENSITY >= 0.75 # Filters for individual frames where the intensity meets a lower threshold.
  )

COLOCALISATION_THREHOLD = 1.5
# Comment: Defines a global variable for the colocalization threshold. This threshold
# is used in the `Cell_Summary_by_Track_for_ColocalisationPercentage` function (as seen in the previous review).

# Dwell Time IRAK4 WT vs IRAK4 DD------------------------------------------------
# Comment: This section prepares for and executes an R script to perform dwell time analysis.
Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "01_Dwell Time")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}
# Creates a dedicated subfolder for dwell time plots.

ScriptList <- c(
  "01_Dwell Time for all Replicates.R" # Dwell Time
)
# Comment: Defines a list of external R scripts to be executed.
# This structure is good for managing multiple related analysis scripts.

# Command
for(x in 1:length(ScriptList)){
  tryCatch({
    print(paste("::::::::::::::::::::", x, "::::::::::::::::::::")) # Prints progress.
    setwd(Plot_Script_Directory) # Sets the working directory to where the scripts are located. This is important for `source()` to find the file.
    source(ScriptList[x], local = T) # Executes the R script. `local = T` ensures variables defined in the sourced script
  }, error = function(e) {print(paste("Error loading", ScriptList[x], ":", e$message))})
  # Comment: Basic error handling using `tryCatch`. If a sourced script fails, it prints an error message
  # but allows the main script to continue, which can be useful for batch processing.
}

rm(
  x,
  ScriptList
)

# Colocalisation IRAK4 WT vs IRAK4 DD----------------------------------------------------------
# Comment: Similar block for colocalization percentage analysis.
Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "02_Colocalisation")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}
# Creates a dedicated subfolder for colocalization plots.

ScriptList <- c(
  "02_Colocalisation Percentage with Dwell Time greater than 30s.R" # Max Normalized Intensity MyD88 Table Creation.R
)
# Note: The comment seems to refer to a different script name ("Max Normalized Intensity MyD88 Table Creation.R")
# than the actual script being sourced ("02_Colocalisation Percentage with Dwell Time greater than 30s.R").
# This might be a copy-paste error or an outdated comment.

LOWER_LIMIT = 0 # The Plot x-axis Lower limit
UPPER_LIMIT = 46 # The Plot a-axis Upper limit
LOWER_LIMIT_AXIS_TICKS = 0 # Axis Lower Limit
X_AXIS_INTERVAL = 20 # Axis tick mark Interval
# Comment: Defines global plotting parameters. These are typically used by `ggplot2` in the sourced scripts.
# Critique:
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
# Comment: Cleans up global plotting parameters and loop variables.

# Cleanup -----------------------------------------------------------------
rm(list = ls()) # Removes all objects from the current R environment.
gc() # Forces garbage collection, freeing up memory.