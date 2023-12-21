library(pacman)
pacman::p_load(dplyr, tidyr, data.table, ggplot2, ggprism, lemon, ggforce)

Table1 <- fread("/Users/u_niranjan/Desktop/Git Scripts/00_Myddosomal_internal_phosphorylation_cohort_table/20_IRAK4WT+IRAK1WT_Compiled_Essential.csv.gz")
Table2 <- fread("/Users/u_niranjan/Desktop/Git Scripts/00_Myddosomal_internal_phosphorylation_cohort_table/21_IRAK4KD+IRAK1WT_Compiled_Essential.csv.gz")
Table3 <- fread("/Users/u_niranjan/Desktop/Git Scripts/00_Myddosomal_internal_phosphorylation_cohort_table/22_IRAK4DD+IRAK1WT_Compiled_Essential.csv.gz")


Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Possible paper figure/Git Figure_20231213" # Where to save Images
Plot_Script_Directory <- "/Users/u_niranjan/Desktop/Git Scripts/IRAK4_Autophosphorylaytion_Paper/03_Figure Plots/06_Figure6" #Where Plot Scripts are are

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
      SHORT_LABEL == "IRAK4 WT + IRAK1 WT" ~ 1,
      SHORT_LABEL == "IRAK4 KD + IRAK1 WT" ~ 2,
      SHORT_LABEL == "IRAK4 DD + IRAK1 WT" ~ 3
    ),
    PROTEIN_1 = case_when(
      PROTEIN == "IRAK4WT" ~ "IRAK4",
      PROTEIN == "IRAK4KD" ~ "IRAK4",
      PROTEIN == "IRAK4DD" ~ "IRAK4"
    )
  ) %>% 
  as.data.table()

# Creating Folder to save Images
Plot_Directory_Path <- file.path(Plot_Directory_Save_Path, "06_Figure 6")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}


# Colocalisation ----------------------------------------------------------
Plot_Directory_Save_Path <- file.path(Plot_Directory_Path, "01_Colocalisation")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

Plot_Script_Sub_Directory <- file.path(Plot_Script_Directory, "03_Figure3")

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

# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()
