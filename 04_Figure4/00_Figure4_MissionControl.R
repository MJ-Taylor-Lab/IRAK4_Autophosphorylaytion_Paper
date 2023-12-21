library(pacman)
pacman::p_load(dplyr, tidyr, data.table, ggplot2, ggprism, lemon, ggforce)

Table1 <- fread("/Users/u_niranjan/Desktop/Git Scripts/00_Myddosomal_internal_phosphorylation_cohort_table/05_IRAK4WT_Compiled_Essential.csv.gz")
Table2 <- fread("/Users/u_niranjan/Desktop/Git Scripts/00_Myddosomal_internal_phosphorylation_cohort_table/06_IRAK4DelKDo_Compiled_Essential.csv.gz")
Table3 <- fread("/Users/u_niranjan/Desktop/Git Scripts/00_Myddosomal_internal_phosphorylation_cohort_table/07_IRAK4KD_Compiled_Essential.csv.gz")

Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Possible paper figure/Git Figure_20231102" # Where to save Images
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

Table <- Table %>% 
  mutate(
    ORDER_NUMBER = case_when(
      SHORT_LABEL == "IRAK4 WT" ~ 1,
      SHORT_LABEL == "IRAK4 KD" ~ 2,
      SHORT_LABEL == "IRAK4 DelKDo" ~ 3
    )
  ) %>% 
  as.data.table()

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Save_Path, "04_Figure 4")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}


# Colocalisation ----------------------------------------------------------
Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "01_Colocalisation")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

Plot_Script_Sub_Directory <- file.path(Plot_Script_Directory, "04_Figure4")

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

Plot_Script_Sub_Directory <- file.path(Plot_Script_Directory, "04_Figure4")

ScriptList <- c(
  "02_Max Normalized Intensity MyD88.R", # Max Normalized Intensity MyD88.R
  "03_Max Normalized Intensity IRAK4.R" # Max Normalized Intensity IRAK4.R
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

Plot_Script_Sub_Directory <- file.path(Plot_Script_Directory, "04_Figure4")

ScriptList <- c(
  "04_Dwell TIme Histogram of Image Averages.R" # Dwell Time
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
