library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/Analysis Output/20221022/12_IRAK2WildType/20230823/Output/Analysis.csv.gz")

Table <- rbind(
  Table1
)

# Checking the max mean and median values to test if calibration failed --------
max(Table1$MAX_NORMALIZED_INTENSITY)
median(Table1$MAX_NORMALIZED_INTENSITY)
mean(Table1$MAX_NORMALIZED_INTENSITY)

unique(Table$COHORT)
unique(Table$PROTEIN)
unique(Table$IMAGE)

# Table Cleanup and save -----------------------------------------------------------
rm(
  Table1
)

Table <- Table %>%  as.data.table()
Table <- Table %>% 
  filter(
    COHORT == "IRAK2WildType"
  ) %>% 
  mutate(
    SHORT_LABEL = "IRAK2 WT"
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

Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/00_Myddosomal_internal_phosphorylation_cohort_table/12_IRAK2WT_Compiled_Essential.csv.gz"
fwrite(Table, Table_path)


# Cleanup -----------------------------------------------------------------
rm(list = ls())
