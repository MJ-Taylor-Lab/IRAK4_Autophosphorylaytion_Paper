library(pacman)
pacman::p_load(dplyr, tidyr, data.table, ggplot2, ggprism, lemon, ggforce)

Table1 <- fread("/Users/u_niranjan/Desktop/Cell_Figures/00_New_Figures/Figure_4/Fixed_Cell_Staining_Manual_Analysis/01_MyD88_pIRAK4_IRAK4/csv_file_path.csv")
Table2 <- fread("/Users/u_niranjan/Desktop/Cell_Figures/00_New_Figures/Figure_4/Fixed_Cell_Staining_Manual_Analysis/02_MyD88_pIRAK4_IRAK1/csv_file_path.csv")

Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots"
Plot_Script_Directory <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/02_Figure_Scripts/04_Figure4/Figure 4A_4B" #Where Plot Scripts are are

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Save_Path, "04_Figure 4")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}


# MyD88-IRAK4-pIRAK4 ----------------------------------------------------------
Table <- Table1

Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "01_Figure_4A")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

ScriptList <- c(
  "01_SpotCount_MyD88_IRAK4_pIRAK4.R"# MNumber of MyD88-IRAK4-pIRAK4 colocalised tracks.R
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


# MyD88-IRAK4-pIRAK4 ----------------------------------------------------------
Table <- Table2

Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "02_Figure_4B")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

ScriptList <- c(
  "02_SpotCount_MyD88_IRAK1_pIRAK4.R"# MNumber of MyD88-IRAK1-pIRAK4 colocalised tracks.R
)
# Commands
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

# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()


