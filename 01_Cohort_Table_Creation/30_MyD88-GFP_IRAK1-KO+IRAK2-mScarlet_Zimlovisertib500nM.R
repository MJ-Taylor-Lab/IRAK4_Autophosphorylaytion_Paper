library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/67_IRAK2WT_PF06650833_500nM/20250703/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/67_IRAK2WT_PF06650833_500nM/20250708_1/Output/Analysis.csv.gz")
Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/67_IRAK2WT_PF06650833_500nM/20250708_2/Output/Analysis.csv.gz")
Table4 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/cobra/Analysis Output/20230727_IRAK4PhosphoPaper/15_IRAK2WT_PF06650833_500nM/20230827/Output/Analysis.csv.gz")
Table5 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/cobra/Analysis Output/20230727_IRAK4PhosphoPaper/15_IRAK2WT_PF06650833_500nM/20230831/Output/Analysis.csv.gz")
Table6 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/cobra/Analysis Output/20230727_IRAK4PhosphoPaper/15_IRAK2WT_PF06650833_500nM/20230905/Output/Analysis.csv.gz")

Table <- rbind(
  Table1,
  Table2,
  Table3,
  Table4,
  Table5,
  Table6
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

max(Table4$MAX_NORMALIZED_INTENSITY)
median(Table4$MAX_NORMALIZED_INTENSITY)
mean(Table4$MAX_NORMALIZED_INTENSITY)

max(Table5$MAX_NORMALIZED_INTENSITY)
median(Table5$MAX_NORMALIZED_INTENSITY)
mean(Table5$MAX_NORMALIZED_INTENSITY)

max(Table6$MAX_NORMALIZED_INTENSITY)
median(Table6$MAX_NORMALIZED_INTENSITY)
mean(Table6$MAX_NORMALIZED_INTENSITY)


unique(Table$COHORT)
unique(Table$PROTEIN)
unique(Table$IMAGE)

# Table Cleanup and save -----------------------------------------------------------
rm(
  Table1,
  Table2,
  Table3,
  Table4,
  Table5,
  Table6
)


Table <- Table %>%  as.data.table()
Table <- Table %>% 
  mutate(
    COHORT = "IRAK2WT_PF06650833_500nM"
  ) %>% 
  # filter(
  #   COHORT == "IRAK2WT_PF06650833_500nM"
  # ) %>% 
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

Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/30_MyD88-GFP_IRAK1-KO+IRAK2-mScarlet_Zimlovisertib500nM_Analysis.csv.gz"
fwrite(Table, Table_path)

rm(list = ls())
gc()
