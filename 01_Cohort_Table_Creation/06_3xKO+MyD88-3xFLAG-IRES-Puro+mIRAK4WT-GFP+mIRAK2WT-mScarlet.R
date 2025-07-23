library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/35_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK2WT/20250114/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/35_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK2WT/20250128/Output/Analysis.csv.gz")
Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/35_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK2WT/20250131/Output/Analysis.csv.gz")
Table4 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/35_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK2WT/20250204/Output/Analysis.csv.gz")

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


# Manually curating cells -------------------------------------------------
# Define the unique IMAGE values with their corresponding counts
image_values <- c(
  rep("20250114 plate2_well2B_20nM_cl571_IRAK4WT_IRAK2WT_001", 12),
  rep("20250114 plate2_well3B_20nM_cl571_IRAK4WT_IRAK2WT_001", 11),
  rep("20250114 plate2_well3C_20nM_cl571_IRAK4WT_IRAK2WT_001", 14),
  rep("20250128 plate2_well8F_20nM_cl571_IRAK4WT_IRAK2WT_001", 13), # Corrected count
  rep("20250128 plate2_well9G_20nM_cl571_IRAK4WT_IRAK2WT_001", 19), # Corrected count
  # rep("20250131 plate1_well11D_20nM_cl571_IRAK4WT_IRAK2WT_001", 9),
  rep("20250204 plate1_well2B_20nM_cl571_IRAK4WT_IRAK2WT_001", 15),
  rep("20250204 plate2_well2F_20nM_cl571_IRAK4WT_IRAK2WT_001", 13), # Corrected count
  rep("20250204 plate2_well2G_20nM_cl571_IRAK4WT_IRAK2WT_001", 12)
)

# Define the corresponding CELL values
cell_values <- c(
  3, 4, 7, 8, 10, 13, 16, 19, 23, 25, 28, 31,
  1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13,
  1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 17, 19,
  2, 3, 5, 8, 10, 11, 13, 16, 17, 18, 19, 20, 21,
  1, 2, 3, 6, 7, 9, 10, 11, 13, 14, 15, 16, 17, 19, 26, 31, 32, 33, 35,
  # 2, 3, 7, 12, 15, 16, 17, 18, 21,
  3, 4, 5, 7, 8, 9, 10, 16, 17, 19, 20, 23, 24, 26, 29,
  6, 8, 9, 10, 14, 18, 19, 22, 23, 24, 26, 30, 31,
  1, 2, 3, 8, 9, 13, 15, 16, 20, 25, 26, 27
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
    IMAGE != "20250204 plate2_well2F_20nM_cl571_IRAK4WT_IRAK2WT_001"
  ) %>% 
  mutate(
    SHORT_LABEL = "IRAK4 WT + IRAK2 WT"
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


Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4_Autophosphorylaytion_Paper_Rewrite/00_Myddosomal_internal_phosphorylation_cohort_table/35_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK2WT.csv.gz"
fwrite(Table, Table_path)

rm(list = ls())
gc()
