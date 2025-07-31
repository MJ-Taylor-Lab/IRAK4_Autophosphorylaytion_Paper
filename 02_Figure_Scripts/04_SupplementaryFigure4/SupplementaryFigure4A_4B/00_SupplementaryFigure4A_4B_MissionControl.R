library(pacman)
pacman::p_load(dplyr, tidyr, data.table, ggplot2, ggprism, lemon, ggforce)

Table1 <- fread("/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/38_MyD88-GFP_IRAK4-mScarlet_Analysis.csv.gz")
Table2 <- fread("/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/39_MyD88-GFP_IRAK1-mScarlet_Analysis.csv.gz")

Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots"
Plot_Script_Directory <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/02_Figure_Scripts/04_SupplementaryFigure4/SupplementaryFigure4A_4B"



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
Plot_Directory_Path <- file.path(Plot_Directory_Save_Path, "04_SupplementaryFigure4")
if(!file.exists(Plot_Directory_Path)){
  dir.create(Plot_Directory_Path)
}


# Colocalisation ----------------------------------------------------------
ScriptList <- c(
  "01_MyD88 IRAK4 Colocalisation Percentage Over Time of Imaging.R",# MNumber of MyD88-IRAK4 colocalised tracks over time.R
  "02_MyD88 IRAK1 Colocalisation Percentage Over Time of Imaging.R" # MNumber of MyD88-IRAK1 colocalised tracks over time.R
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

# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()


# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()


