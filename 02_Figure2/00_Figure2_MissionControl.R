library(pacman)
pacman::p_load(dplyr, tidyr, data.table, ggplot2, ggprism, lemon, ggforce)

Table1 <- fread("/Users/u_niranjan/Desktop/Possible paper figure/New Figure 27_06_2023/Figure 2/Immunostaining Pipeline Pilot/01_Manual_Analysis/01_MyD88_pIRAK4_IRAK4/csv_file_path.csv")
Table2 <- fread("/Users/u_niranjan/Desktop/Possible paper figure/New Figure 27_06_2023/Figure 2/Immunostaining Pipeline Pilot/01_Manual_Analysis/02_MyD88_pIRAK4_IRAK1/csv_file_path.csv")

Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Possible paper figure/Git Figure_20230807" # Where to save Images
Plot_Script_Directory <- "/Users/u_niranjan/Desktop/Git Scripts/A_Myddosomal_internal_phosphorylation_cascade_regulates_assembly" #Where Plot Scripts are are

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Save_Path, "02_Figure 2")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}


# MyD88-IRAK4-pIRAK4 ----------------------------------------------------------
Table <- Table1

Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "01_SpotCount_MyD88_IRAK4_pIRAK4")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

Plot_Script_Sub_Directory <- file.path(Plot_Script_Directory, "02_Figure2")

ScriptList <- c(
  "01_SpotCount_MyD88_IRAK4_pIRAK4.R"# MNumber of MyD88-IRAK4-pIRAK4 colocalised tracks.R
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


# MyD88-IRAK4-pIRAK4 ----------------------------------------------------------
Table <- Table2

Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "02_SpotCount_MyD88_IRAK1_pIRAK4")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

Plot_Script_Sub_Directory <- file.path(Plot_Script_Directory, "02_Figure2")

ScriptList <- c(
  "02_SpotCount_MyD88_IRAK1_pIRAK4.R"# MNumber of MyD88-IRAK1-pIRAK4 colocalised tracks.R
)
# Commands
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


