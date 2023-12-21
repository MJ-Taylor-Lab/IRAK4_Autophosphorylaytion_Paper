library(pacman)
pacman::p_load(dplyr, tidyr, data.table, ggplot2, ggprism, lemon, ggforce)

Table1 <- fread("/Users/u_niranjan/Desktop/Git Scripts/00_Myddosomal_internal_phosphorylation_cohort_table/17_IRAK3WT_DMSO_Compiled_Essential.csv.gz")
Table2 <- fread("/Users/u_niranjan/Desktop/Git Scripts/00_Myddosomal_internal_phosphorylation_cohort_table/18_IRAK3WT_PF06650833_20uM_Compiled_Essential.csv.gz")
Table3 <- fread("/Users/u_niranjan/Desktop/Git Scripts/00_Myddosomal_internal_phosphorylation_cohort_table/19_IRAK3WT_PF06650833_500nM_Compiled_Essential.csv.gz")

Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Possible paper figure/Git Figure_20230925" # Where to save Images
Plot_Script_Directory <- "/Users/u_niranjan/Desktop/Git Scripts/IRAK4_Autophosphorylaytion_Paper/03_Figure Plots/07_Figure7/02_IRAK3_DMSOvsIRAK4KinaseInhibitor" #Where Plot Scripts are are

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
  # filter(
  #   IMAGE != "20230827 plate01_well5B_5nM_cl356_IRAK3WT_DMSO_002",
  #   IMAGE != "20230831 plate01_well6G_5nM_cl356_IRAK3WT_DMSO_001"
  # ) %>%
  mutate(
    ORDER_NUMBER = case_when(
      SHORT_LABEL == "DMSO" ~ 1,
      SHORT_LABEL == "Kinase Inhibitor 500 nM" ~ 2,
      SHORT_LABEL == "Kinase Inhibitor 20 uM" ~ 3
    )
  ) %>% 
  as.data.table()

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Save_Path, "07_Figure 7")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}


# Colocalisation ----------------------------------------------------------
Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "05_IRAK3Colocalisation")
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
Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "06_IRAK3Max Normalised Intensity")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

Plot_Script_Sub_Directory <- file.path(Plot_Script_Directory, "05_Figure5")

ScriptList <- c(
  "02_Max Normalized Intensity average.R" # Max Normalized Intensity IRAK1 Average.R
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
Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "07_IRAK3Dwell Time")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

Plot_Script_Sub_Directory <- file.path(Plot_Script_Directory, "05_Figure5")

ScriptList <- c(
  "03_Dwell Time for single Replicate.R", # Dwell Time for Individual Replicate.R
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

