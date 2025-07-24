library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/62_IRAK3WT_PF06650833_500nM/20250715/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/62_IRAK3WT_PF06650833_500nM/20250722_1/Output/Analysis.csv.gz")
Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/62_IRAK3WT_PF06650833_500nM/20250722_2/Output/Analysis.csv.gz")

Table <- rbind(
  Table1,
  Table2,
  Table3
)

# Checking the max mean and median values to test if calibration failed --------
max(Table1$MAX_NORMALIZED_INTENSITY)
median(Table1$MAX_NORMALIZED_INTENSITY)
mean(Table1$MAX_NORMALIZED_INTENSITY)

max(Table2$MAX_NORMALIZED_INTENSITY)
median(Table2$MAX_NORMALIZED_INTENSITY)
mean(Table2$MAX_NORMALIZED_INTENSITY)

max(Table3$MAX_NORMALIZED_INTENSITY)
median(Table3$MAX_NORMALIZED_INTENSITY)
mean(Table3$MAX_NORMALIZED_INTENSITY)

unique(Table$COHORT)
unique(Table$PROTEIN)
unique(Table$IMAGE)

# Table Cleanup and save -----------------------------------------------------------
rm(
  Table1,
  Table2,
  Table3
)


Table <- Table %>%  as.data.table()
Table <- Table %>% 
  filter(
    COHORT == "IRAK3WT_PF06650833_500nM",
    IMAGE != "20250722 plate01_well2E_10nM_cl359_IRAK3WT_inhi500nm_001"
  ) %>% 
  mutate(
    SHORT_LABEL = "Kinase Inhibitor 500 nM"
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

unique(Table$IMAGE)

Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/35_MyD88-GFP_IRAK1-KO+IRAK3-mScarlet_Zimlovisertib500nM_Analysis.csv.gz"
fwrite(Table, Table_path)


# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()
