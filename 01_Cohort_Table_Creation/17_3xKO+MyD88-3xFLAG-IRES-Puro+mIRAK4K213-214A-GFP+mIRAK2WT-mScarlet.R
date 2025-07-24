library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/46_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4K213-214A_mIRAK2WT/20250225/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/46_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4K213-214A_mIRAK2WT/20250228/Output/Analysis.csv.gz")
Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/46_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4K213-214A_mIRAK2WT/20250304/Output/Analysis.csv.gz")
Table4 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/46_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4K213-214A_mIRAK2WT/20250314/Output/Analysis.csv.gz")
Table5 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/46_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4K213-214A_mIRAK2WT/20250325/Output/Analysis.csv.gz")
Table6 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/46_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4K213-214A_mIRAK2WT/20250328/Output/Analysis.csv.gz")

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


# # Manually curating cells -------------------------------------------------
# Define the unique IMAGE values with their corresponding counts
image_values <- c(
  rep("20250225 plate2_well3G_20mM_cl622_IRAK4K213-214A_IRAK2WT_", 9),
  rep("20250228 plate2_well5E_20mM_cl622_IRAK4K213-214A_IRAK2WT_00", 7),
  rep("20230304 plate2_well2F_20mM_cl622_IRAK4K213-214A_IRAK2WT_", 9),
  rep("20250314 plate2_well9E_40mM_cl622_IRAK4K213-214A_IRAK2WT_", 15),
  rep("20250325 plate2_well5G_40mM_cl622_IRAK4K213-214A_IRAK2WT_001", 2),
  rep("20250325 plate2_well5G_40mM_cl622_IRAK4K213-214A_IRAK2WT_002", 12),
  rep("20250328 plate1_well3D_20mM_cl622_IRAK4K213-214A_IRAK2WT_001", 4)
)
# Define the corresponding CELL values
cell_values <- c(
  4, 6, 10, 11, 15, 16, 18, 19, 20,
  4, 8, 9, 14, 27, 28, 29,
  4, 9, 10, 15, 17, 18, 25, 31, 34,
  3, 4, 6, 7, 8, 11, 12, 16, 17, 18, 19, 21, 22, 23, 24,
  24, 25,
  1, 2, 3, 5, 7, 10, 11, 12, 13, 14, 15, 16,
  1, 2, 7, 10
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
    SHORT_LABEL = "IRAK4 K213-214A + IRAK2 WT"
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

Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/17_3xKO+MyD88-3xFLAG-IRES-Puro+mIRAK4K213-214A-GFP+mIRAK2WT-mScarlet_Analysis.csv.gz"
fwrite(Table, Table_path)

rm(list = ls())
gc()