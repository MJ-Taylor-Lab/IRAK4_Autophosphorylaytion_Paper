library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/38_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4DD_mIRAK3WT/20250114/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/38_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4DD_mIRAK3WT/20250128/Output/Analysis.csv.gz")
Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/38_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4DD_mIRAK3WT/20250131/Output/Analysis.csv.gz")
Table4 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/38_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4DD_mIRAK3WT/20250204/Output/Analysis.csv.gz")

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
  rep("20250114 plate1_well2C_20nM_cl574_IRAK4DD_IRAK3WT_001", 11),
  rep("20250128 plate1_well8C_20nM_cl574_IRAK4DD_IRAK3WT_001", 10),
  rep("20250128 plate1_well9B_20nM_cl574_IRAK4DD_IRAK3WT_001", 13),
  rep("20250131 plate2_well10D_20nM_cl574_IRAK4DD_IRAK3WT_001", 5),
  rep("20250204 plate1_well2C_20nM_cl574_IRAK4DD_IRAK3WT_001", 8),
  rep("20250204 plate1_well3B_20nM_cl574_IRAK4DD_IRAK3WT_001", 4)
)

# Define the corresponding CELL values
cell_values <- c(
  3, 5, 10, 12, 14, 15, 20, 24, 25, 27, 28,
  1, 6, 8, 9, 10, 15, 19, 20, 22, 23, 
  2, 4, 6, 7, 8, 9, 11, 12, 14, 17, 18, 21, 22, 
  11, 12, 19, 27, 29, 
  3, 7, 12, 18, 21, 26, 29, 30, 
  1, 3, 7, 12
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
  filter(
    COHORT == "TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4DD_mIRAK3WT"
  ) %>% 
  mutate(
    SHORT_LABEL = "IRAK4 DD + IRAK3 WT"
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

Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/09_3xKO+MyD88-3xFLAG-IRES-Puro+mIRAK4DD-GFP+mIRAK3WT-mScarlet_Analysis.csv.gz"
fwrite(Table, Table_path)

rm(list = ls())
gc()
