library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/Analysis Output/20230727_IRAK4PhosphoPaper/09_IRAK1KI_PF06650833_20uM/20220211/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/Analysis Output/20230727_IRAK4PhosphoPaper/09_IRAK1KI_PF06650833_20uM/20220310/Output/Analysis.csv.gz")
Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/Analysis Output/20230727_IRAK4PhosphoPaper/09_IRAK1KI_PF06650833_20uM/20230414/Output/Analysis.csv.gz")

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

# Table -------------------------------------------------------------------
Table <- Table %>%  as.data.table()
Table <- Table %>% 
  filter(
    COHORT == "IRAK1KI+PF06650833_20uM"
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

Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/00_Myddosomal_internal_phosphorylation_cohort_table/09_IRAK1KI_PF06650833_20uM_Compiled_Essential.csv.gz"
fwrite(Table, Table_path)

rm(list = ls())
