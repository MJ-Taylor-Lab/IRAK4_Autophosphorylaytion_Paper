library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

# Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/60_IRAK3WT_DMSO/20231010/Output/Analysis.csv.gz")
# Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/60_IRAK3WT_DMSO/20231013/Output/Analysis.csv.gz")
# Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/60_IRAK3WT_DMSO/20231019/Output/Analysis.csv.gz")
# Table4 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/60_IRAK3WT_DMSO/20250416/Output/Analysis.csv.gz")
# Table5 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/60_IRAK3WT_DMSO/20250422/Output/Analysis.csv.gz")
Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/60_IRAK3WT_DMSO/20250715/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/60_IRAK3WT_DMSO/20250722_1/Output/Analysis.csv.gz")
Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/60_IRAK3WT_DMSO/20250722_2/Output/Analysis.csv.gz")


Table <- rbind(
  Table1,
  Table2,
  Table3
  # Table4,
  # Table5
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

# max(Table4$MAX_NORMALIZED_INTENSITY)
# median(Table4$MAX_NORMALIZED_INTENSITY)
# mean(Table4$MAX_NORMALIZED_INTENSITY)
# 
# max(Table5$MAX_NORMALIZED_INTENSITY)
# median(Table5$MAX_NORMALIZED_INTENSITY)
# mean(Table5$MAX_NORMALIZED_INTENSITY)

unique(Table$COHORT)
unique(Table$PROTEIN)
unique(Table$IMAGE)

# Table Cleanup and save -----------------------------------------------------------
rm(
  Table1,
  Table2,
  Table3
  # Table4,
  # Table5
)

Table <- Table %>%  as.data.table()
Table <- Table %>% 
  filter(
    COHORT == "IRAK3WT+DMSO"
    # IMAGE != "20231013 plate02_well5F_5nM_cl359_IRAK3WT_DMSO_001"
  ) %>% 
  mutate(
    SHORT_LABEL = "DMSO"
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

Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4_Autophosphorylaytion_Paper_Rewrite/00_Myddosomal_internal_phosphorylation_cohort_table/62_IRAK3WT_DMSO_Compiled_Essential.csv.gz"
fwrite(Table, Table_path)


# Cleanup -----------------------------------------------------------------
rm(list = ls())
