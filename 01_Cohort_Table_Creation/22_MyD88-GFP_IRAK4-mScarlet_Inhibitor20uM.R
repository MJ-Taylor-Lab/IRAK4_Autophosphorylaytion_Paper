library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/cobra/Analysis Output/20230727_IRAK4PhosphoPaper/04_IRAK4KI_PF06650833_20uM/20211008/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/cobra/Analysis Output/20230727_IRAK4PhosphoPaper/04_IRAK4KI_PF06650833_20uM/20220110/Output/Analysis.csv.gz")

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

# Table -------------------------------------------------------------------
Table <- Table %>%  as.data.table()
Table <- Table %>% 
  filter(
    COHORT == "IRAK4KI+PF06650833_20uM",
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 != COMPLEMENTARY_TOTAL_INTENSITY_1
  ) %>% 
  mutate(
    SHORT_LABEL = "Kinase Inhibitor 20 uM"
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

Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/22_MyD88-GFP_IRAK4-mScarlet_Zimlovisertib20uM_Analysis.csv.gz"
fwrite(Table, Table_path)

rm(list = ls())
gc()
