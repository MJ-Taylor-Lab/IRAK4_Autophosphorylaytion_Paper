library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/44_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4K213-214A_mIRAK1WT/20250314/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/44_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4K213-214A_mIRAK1WT/20250318/Output/Analysis.csv.gz")
Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/44_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4K213-214A_mIRAK1WT/20250328/Output/Analysis.csv.gz")
Table4 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/44_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4K213-214A_mIRAK1WT/20250401_01/Output/Analysis.csv.gz")
Table5 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/44_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4K213-214A_mIRAK1WT/20250401_02/Output/Analysis.csv.gz")

Table <- rbind(
  Table1,
  Table2,
  Table3,
  Table4,
  Table5
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


unique(Table$COHORT)
unique(Table$PROTEIN)
unique(Table$IMAGE)

# Table Cleanup and save -----------------------------------------------------------
rm(
  Table1,
  Table2,
  Table3,
  Table4,
  Table5
)


# # Manually curating cells -------------------------------------------------
# Define the unique IMAGE values with their corresponding counts
image_values <- c(
  rep("20250314 plate1_well8F_40mM_cl620_IRAK4K213-214A_IRAK1WT_00", 6),
  rep("20250314 plate1_well8F_40mM_cl620_IRAK4K213-214A_IRAK1WT_001", 15),
  rep("20250318 plate2_well2E_40mM_cl620_IRAK4K213-214A_IRAK1WT_002", 6),
  rep("20250328 plate1_well2C_40mM_cl620_IRAK4K213-214A_IRAK1WT_001", 8),
  rep("20250328 plate2_well2G_40mM_cl620_IRAK4K213-214A_IRAK1WT_001", 7),
  rep("20250328 plate2_well3F_40mM_cl620_IRAK4K213-214A_IRAK1WT_001", 9),
  rep("20250401 plate1_well6D_40mM_cl620_IRAK4K213-214A_IRAK1WT_001", 13),
  rep("20250401 plate1_well5C_40mM_cl620_IRAK4K213-214A_IRAK1WT_", 20)
)
# # Define the corresponding CELL values
cell_values <- c(
  4, 6, 11, 17, 23, 26, 
  3, 4, 6, 7, 8, 11, 12, 14, 15, 19, 22, 23, 25, 29, 31, 
  2, 4, 5, 6, 8, 9, 
  2, 3, 6, 7, 8, 14, 17, 18, 
  6, 8, 10, 11, 12, 14, 19, 
  1, 11, 12, 14, 15, 17, 18, 19, 21, 
  4, 5, 6, 7, 10, 12, 13, 14, 16, 17, 19, 20, 22, 
  1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 19, 21, 24, 25, 26, 27, 28
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
    SHORT_LABEL = "IRAK4 K213-214A + IRAK1 WT"
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


Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4_Autophosphorylaytion_Paper_Rewrite/00_Myddosomal_internal_phosphorylation_cohort_table/51_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4K213-214A_mIRAK1WT.csv.gz"
fwrite(Table, Table_path)

rm(list = ls())
gc()
