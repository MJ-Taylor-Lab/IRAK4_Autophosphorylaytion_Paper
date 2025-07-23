library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/Analysis Output/20230727_IRAK4PhosphoPaper/20_IRAK4WT_IRAK1WT/20231212/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/Analysis Output/20230727_IRAK4PhosphoPaper/20_IRAK4WT_IRAK1WT/20231219/Output/Analysis.csv.gz")


Table <- rbind(
  Table1,
  Table2
)

# Checking the max mean and median values to test if calibration failed --------
max(Table1$MAX_NORMALIZED_INTENSITY)
median(Table1$MAX_NORMALIZED_INTENSITY)
mean(Table1$MAX_NORMALIZED_INTENSITY)

max(Table2$MAX_NORMALIZED_INTENSITY)
median(Table2$MAX_NORMALIZED_INTENSITY)
mean(Table2$MAX_NORMALIZED_INTENSITY)

unique(Table$COHORT)
unique(Table$PROTEIN)
unique(Table$IMAGE)

# Table Cleanup and save -----------------------------------------------------------
rm(
  Table1,
  Table2
)

Table <- Table %>%  as.data.table()
Table <- Table %>% 
  filter(
    COHORT == "TKO+MyD88-3xFLAG+IRAK4WT+IRAK1WT"
  ) %>% 
  mutate(
    SHORT_LABEL = "IRAK4 WT + IRAK1 WT"
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    LIFETIME = max(TIME_ADJUSTED),
    MAX_NORMALIZED_INTENSITY = max(NORMALIZED_INTENSITY),
    MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1 = max(COMPLEMENTARY_NORMALIZED_INTENSITY_1)
  ) %>% 
  arrange(
    UNIVERSAL_SPOT_ID
  ) %>%
  as.data.table()

# Exclusion List ----------------------------------------------------------
# Excluding data from cells that were wrongly identified
ManualExclusionList_1 <- data_frame(
  UNIVERSAL_CELL_ID = "20231212 plate01_well2F_10nM_cl_D11_IRAK4WT_IRAK1WT_001_",
  CELL = c(6,7,8,10,12,13,17,23,24)
) 
ManualExclusionList_1 <- ManualExclusionList_1 %>% 
  mutate(
    UNIVERSAL_CELL_ID = paste0(UNIVERSAL_CELL_ID, CELL)
  )


ManualExclusionList_2 <- data_frame(
  UNIVERSAL_CELL_ID = "20231212 plate02_well2F_10nM_cl_D11_IRAK4WT_IRAK1WT_001_",
  CELL = c(5,6,7,9,12,13,14,17,18,20,21)
) 
ManualExclusionList_2 <- ManualExclusionList_2 %>% 
  mutate(
    UNIVERSAL_CELL_ID = paste0(UNIVERSAL_CELL_ID, CELL)
  )

ManualExclusionList_3 <- data_frame(
  UNIVERSAL_CELL_ID = "20231219 plate01_well3B_10nM_cl_D11_IRAK4WT_IRAK1WT_002_",
  CELL = c(1,5,7,8,9,13,14,15,16)
) 
ManualExclusionList_3 <- ManualExclusionList_3 %>% 
  mutate(
    UNIVERSAL_CELL_ID = paste0(UNIVERSAL_CELL_ID, CELL)
  )


ManualExclusionList_4 <- data_frame(
  UNIVERSAL_CELL_ID = "20231219 plate02_well2G_10nM_cl_D11_IRAK4WT_IRAK1WT_001_",
  CELL = c(1,2,3,4,5,6,8,9,10,11,14,15,16,18,20,21,22,23,24,25,26,27,28,29)
) 
ManualExclusionList_4 <- ManualExclusionList_4 %>% 
  mutate(
    UNIVERSAL_CELL_ID = paste0(UNIVERSAL_CELL_ID, CELL)
  )

ManualExclusionList <- rbind(
  ManualExclusionList_1,
  ManualExclusionList_2,
  ManualExclusionList_3,
  ManualExclusionList_4
)

rm(
  ManualExclusionList_1,
  ManualExclusionList_2,
  ManualExclusionList_3,
  ManualExclusionList_4
)


# Removing incorrectly identified cells -----------------------------------
Table <- Table %>% 
  mutate(
    UNIVERSAL_CELL_ID = paste0(IMAGE, "_", CELL),
    UNIVERSAL_CELL_ID = UNIVERSAL_CELL_ID %in% ManualExclusionList$UNIVERSAL_CELL_ID
  ) %>% 
  filter(
    UNIVERSAL_CELL_ID == FALSE
  ) %>% 
  select(
    -UNIVERSAL_CELL_ID
  )

Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4_Autophosphorylaytion_Paper_Rewrite/00_Myddosomal_internal_phosphorylation_cohort_table/03_3xKO_MyD88-IRES-Puro_IRAK4WT_IRAK1WT_Compiled_Essential.csv.gz"
fwrite(Table, Table_path)

rm(list = ls())
