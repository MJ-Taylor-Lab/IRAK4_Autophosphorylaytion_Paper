library(pacman)
pacman::p_load(dplyr, tidyr, data.table, ggplot2, ggprism, lemon, ggforce)

Table1 <- fread("/Users/u_niranjan/Desktop/Git Scripts/00_Myddosomal_internal_phosphorylation_cohort_table/08_IRAK1KI_DMSO_Compiled_Essential.csv.gz")
Table2 <- fread("/Users/u_niranjan/Desktop/Git Scripts/00_Myddosomal_internal_phosphorylation_cohort_table/09_IRAK1KI_PF06650833_20uM_Compiled_Essential.csv.gz")
Table3 <- fread("/Users/u_niranjan/Desktop/Git Scripts/00_Myddosomal_internal_phosphorylation_cohort_table/10_IRAK1KI_PF06650833_500nM_Compiled_Essential.csv.gz")

Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Possible paper figure/Git Figure_20230807" # Where to save Images
Plot_Script_Directory <- "/Users/u_niranjan/Desktop/Git Scripts/A_Myddosomal_internal_phosphorylation_cascade_regulates_assembly" #Where Plot Scripts are are

Table <- rbind(
  Table1,
  Table2,
  Table3
)

unique(Table$COHORT)
unique(Table$SHORT_LABEL)

rm(
  Table1,
  Table2,
  Table3
)

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Save_Path, "05_Figure 5")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}


# Colocalisation ----------------------------------------------------------
Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "01_Colocalisation")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

Plot_Script_Sub_Directory <- file.path(Plot_Script_Directory, "05_Figure5")

ScriptList <- c(
  "01_Colocalisation Percentage Sweep by Dwell Time Violin.R" # Max Normalized Intensity MyD88 Table Creation.R
)

# Command
for(x in 1:length(ScriptList)){
  tryCatch({
    print(paste("::::::::::::::::::::", x, "::::::::::::::::::::"))
    setwd(Plot_Script_Sub_Directory)
    source(ScriptList[x], local = T)
  }, error = function(e) {print(paste("Error loading", ScriptList[x]))})
}

rm(
  x,
  ScriptList
)


# Max Normalized Intensity ------------------------------------------------
Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "02_Max Normalised Intensity")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

Plot_Script_Sub_Directory <- file.path(Plot_Script_Directory, "05_Figure5")

ScriptList <- c(
  "02_Max Normalized Intensity single Replicate.R", # Max Normalized Intensity IRAK1 Individual Replicate.R
  "03_Max Normalized Intensity average.R" # Max Normalized Intensity IRAK1 Average.R
)

# Command
for(x in 1:length(ScriptList)){
  tryCatch({
    print(paste("::::::::::::::::::::", x, "::::::::::::::::::::"))
    setwd(Plot_Script_Sub_Directory)
    source(ScriptList[x], local = T)
  }, error = function(e) {print(paste("Error loading", ScriptList[x]))})
}

rm(
  x,
  ScriptList
)

# Dwell Time ------------------------------------------------
Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "03_Dwell Time")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

Plot_Script_Sub_Directory <- file.path(Plot_Script_Directory, "05_Figure5")

ScriptList <- c(
  "04_Dwell Time for single Replicate.R", # Dwell Time for Individual Replicate.R
  "05_Percentage of tracks with Dwell Time greater than 30s.R", 
  "06_Percentage of tracks that do recruit with Dwell Time greater than 30s.R"
)

# Command
for(x in 1:length(ScriptList)){
  tryCatch({
    print(paste("::::::::::::::::::::", x, "::::::::::::::::::::"))
    setwd(Plot_Script_Sub_Directory)
    source(ScriptList[x], local = T)
  }, error = function(e) {print(paste("Error loading", ScriptList[x]))})
}

rm(
  x,
  ScriptList
)


# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()

