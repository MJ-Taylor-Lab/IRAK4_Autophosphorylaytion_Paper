library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/36_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4DD_mIRAK2WT/20250114/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/36_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4DD_mIRAK2WT/20250128/Output/Analysis.csv.gz")
Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/36_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4DD_mIRAK2WT/20250131/Output/Analysis.csv.gz")
Table4 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/36_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4DD_mIRAK2WT/20250204/Output/Analysis.csv.gz")

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
# 
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

# Manually curating cells -------------------------------------------------
# Define the unique IMAGE values with their corresponding counts
image_values <- c(
  rep("20250114 plate1_well2B_20nM_cl572_IRAK4DD_IRAK2WT_001", 9),
  rep("20250114 plate2_well2C_20nM_cl572_IRAK4DD_IRAK2WT_001", 21),
  rep("20250128 plate2_well8G_20nM_cl572_IRAK4DD_IRAK2WT_001", 15),
  rep("20250128 plate2_well9F_20nM_cl572_IRAK4DD_IRAK2WT_001", 10),
  rep("20250131 plate1_well11B_20nM_cl572_IRAK4DD_IRAK2WT_001", 5),
  rep("20250204 plate2_well3F_20nM_cl572_IRAK4DD_IRAK2WT_001", 8),
  rep("20250204 plate2_well3G_20nM_cl572_IRAK4DD_IRAK2WT_001", 10)
)

# Define the corresponding CELL values
cell_values <- c(
  1, 2, 5, 6, 7, 12, 13, 16, 18, 
  5, 8, 10, 11, 13, 16, 20, 22, 24, 25, 27, 29, 31, 32, 34, 35, 36, 38, 39, 40, 41,
  3, 6, 8, 11, 13, 14, 15, 22, 23, 25, 26, 27, 28, 30, 31, 
  2, 3, 4, 7, 11, 14, 18, 19, 20, 22, 
  3, 5, 7, 9, 14, 
  8, 16, 20, 21, 22, 23, 26, 28,
  7, 8, 10, 15, 19, 20, 21, 23, 26, 28
)

# Create the Selection_Table data frame
Selection_Table <- data.frame(
  IMAGE = image_values,
  CELL = cell_values
)

rm(
  image_values,
  cell_values
)


# Making Final Table ------------------------------------------------------
Table <- Table %>%  as.data.table()
Table <- Table %>% 
  inner_join(
    Selection_Table, by = c("IMAGE", "CELL")
  ) %>% 
  mutate(
    SHORT_LABEL = "IRAK4 DD + IRAK2 WT"
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

Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/07_3xKO+MyD88-3xFLAG-IRES-Puro+mIRAK4DD-GFP+mIRAK2WT-mScarlet_Analysis.csv.gz"
fwrite(Table, Table_path)

rm(list = ls())
gc()
