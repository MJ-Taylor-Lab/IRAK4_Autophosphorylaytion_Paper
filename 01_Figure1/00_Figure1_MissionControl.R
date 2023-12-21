library(pacman)
pacman::p_load(dplyr, tidyr, data.table, ggplot2, ggprism, lemon, ggforce)

Table1 <- fread("/Users/u_niranjan/Desktop/Git Scripts/00_Myddosomal_internal_phosphorylation_cohort_table/01_IRAK4KI_Compiled_Essential.csv.gz")
Table2 <- fread("/Users/u_niranjan/Desktop/Git Scripts/00_Myddosomal_internal_phosphorylation_cohort_table/02_IRAK1KI_Compiled_Essential.csv.gz")

Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Possible paper figure/Git Figure_20230925" # Where to save Images
Plot_Script_Directory <- "/Users/u_niranjan/Desktop/Git Scripts/A_Myddosomal_internal_phosphorylation_cascade_regulates_assembly" #Where Plot Scripts are are

Table <- rbind(
  Table1,
  Table2
)

unique(Table$COHORT)
unique(Table$SHORT_LABEL)

rm(
  Table1,
  Table2
)

Table <- Table %>% 
  mutate(
    ORDER_NUMBER = case_when(
      SHORT_LABEL == "IRAK4" ~ 1,
      SHORT_LABEL == "IRAK1" ~ 2
    )
  ) %>% 
  as.data.table()

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Save_Path, "01_Figure 1")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}


# Colocalisation ----------------------------------------------------------
Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "01_Colocalisation")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

Plot_Script_Sub_Directory <- file.path(Plot_Script_Directory, "01_Figure1")

ScriptList <- c(
  "01_MyD88 IRAK4 Colocalisation Percentage Over Time of Imaging.R",# MNumber of MyD88-IRAK4 colocalised tracks over time.R
  "02_MyD88 IRAK1 Colocalisation Percentage Over Time of Imaging.R" # MNumber of MyD88-IRAK1 colocalised tracks over time.R
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


# Western Blot Analysis ---------------------------------------------------
Table <- fread("/Users/u_niranjan/Desktop/Possible paper figure/New Figure 27_06_2023/Figure 1/p-IRAK4 blot_rep6.csv")
Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Possible paper figure/Git Figure_20230807" # Where to save Images

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Save_Path, "01_Figure 1")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}

Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "02_Western Blot Analysis")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()


