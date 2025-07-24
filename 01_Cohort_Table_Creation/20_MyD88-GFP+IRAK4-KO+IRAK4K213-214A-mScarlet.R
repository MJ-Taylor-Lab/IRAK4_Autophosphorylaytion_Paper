library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/cobra/Analysis Output/20230727_IRAK4PhosphoPaper/07_IRAK4_KinaseDead/20211129/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/cobra/Analysis Output/20230727_IRAK4PhosphoPaper/07_IRAK4_KinaseDead/20220824/Output/Analysis.csv.gz")
Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/cobra/Analysis Output/20230727_IRAK4PhosphoPaper/07_IRAK4_KinaseDead/20230103/Output/Analysis.csv.gz")
Table4 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/cobra/Analysis Output/20230727_IRAK4PhosphoPaper/07_IRAK4_KinaseDead/20231031/Output/Analysis.csv.gz")


Table <- rbind(
  Table1,
  Table2,
  Table3,
  Table4
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

unique(Table$COHORT)
unique(Table$PROTEIN)
unique(Table$IMAGE)

# Table Cleanup and save -----------------------------------------------------------
rm(
  Table1,
  Table2,
  Table3,
  Table4
)

Table <- Table %>%  as.data.table()
Table <- Table %>% 
  filter(
    COHORT == "IRAK4KinaseDead"
  ) %>% 
  mutate(
    SHORT_LABEL = "IRAK4 K213/214A"
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

Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/20_MyD88-GFP+IRAK4-KO+IRAK4K213-214A-mScarlet_Analysis.csv.gz"
fwrite(Table, Table_path)

rm(list = ls())
