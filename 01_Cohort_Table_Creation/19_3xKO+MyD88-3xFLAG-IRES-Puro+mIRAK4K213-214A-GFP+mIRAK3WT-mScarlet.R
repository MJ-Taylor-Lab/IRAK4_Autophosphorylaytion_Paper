library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/48_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4K213-214A_mIRAK3WT/20250314/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/48_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4K213-214A_mIRAK3WT/20250328/Output/Analysis.csv.gz")
Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/48_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4K213-214A_mIRAK3WT/20250401_01/Output/Analysis.csv.gz")
Table4 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/48_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4K213-214A_mIRAK3WT/20250401_02/Output/Analysis.csv.gz")
Table5 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/48_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4K213-214A_mIRAK3WT/20250416/Output/Analysis.csv.gz")

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
  rep("20250314 plate1_well9G_40mM_cl623_IRAK4K213-214A_IRAK3WT_", 1),
  rep("20250328 plate1_well3B_40mM_cl623_IRAK4K213-214A_IRAK3WT_001", 10),
  rep("20250328 plate2_well3G_40mM_cl623_IRAK4K213-214A_IRAK3WT_001", 12),
  rep("20250401 plate1_well6B_40mM_cl623_IRAK4K213-214A_IRAK3WT_", 9),
  rep("20250401 plate2_well3E_40mM_cl623_IRAK4K213-214A_IRAK3WT_", 10)
)

# # Define the corresponding CELL values
cell_values <- c(
  4,
  2, 5, 6, 8, 12, 13, 20, 22, 23, 29, 
  1, 6, 9, 11, 13, 14, 20, 21, 23, 24, 25, 26,
  1, 3, 4, 11, 13, 14, 18, 23, 24,
  1, 4, 5, 6, 7, 9, 11, 14, 16, 17
)

# # Create the Selection_Table data frame
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
    SHORT_LABEL = "IRAK4 K213-214A + IRAK3 WT"
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

Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/19_3xKO+MyD88-3xFLAG-IRES-Puro+mIRAK4K213-214A-GFP+mIRAK3WT-mScarlet_Analysis.csv.gz"
fwrite(Table, Table_path)

rm(list = ls())
gc()