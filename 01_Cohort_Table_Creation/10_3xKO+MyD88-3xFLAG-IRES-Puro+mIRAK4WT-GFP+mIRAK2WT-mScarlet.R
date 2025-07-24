library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/39_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK2WT/20250207/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/39_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK2WT/20250211/Output/Analysis.csv.gz")

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


# # Manually curating cells -------------------------------------------------
# Define the unique IMAGE values with their corresponding counts
image_values <- c(
  rep("20250207 plate1_well4D_20nM_cl571_IRAK4WT_IRAK2WT_001", 8),
  rep("20250207 plate1_well5B_20nM_cl571_IRAK4WT_IRAK2WT_001", 6),
  rep("20250207 plate2_well5E_20nM_cl571_IRAK4WT_IRAK2WT_001", 5),
  rep("20250207 plate2_well5F_20nM_cl571_IRAK4WT_IRAK2WT_001", 13),
  rep("20250207 plate2_well5G_20nM_cl571_IRAK4WT_IRAK2WT_001", 11),
  rep("20250211 plate1_well8B_20nM_cl571_IRAK4WT_IRAK2WT_001", 16),
  rep("20250211 plate2_well7F_20nM_cl571_IRAK4WT_IRAK2WT_001", 4),
  rep("20250211 plate2_well8F_20nM_cl571_IRAK4WT_IRAK2WT_001", 17)
)

# Define the corresponding CELL values
cell_values <- c(
  1, 2, 8, 9, 10, 11, 15, 18,
  7, 10, 11, 12, 18, 20,
  1, 8, 11, 15, 20,
  4, 5, 6, 7, 10, 12, 13, 16, 18, 25, 28, 30, 32,
  1, 4, 5, 6, 8, 9, 17, 21, 23, 25, 26, 
  2, 3, 4, 10, 11, 18, 19, 20, 21, 22, 23, 24, 27, 32, 33, 36,
  1, 6, 7, 13, 
  2, 4, 6, 9, 13, 14, 15, 17, 19, 20, 21, 22, 23, 24, 25, 29, 37
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

unique(Table$IMAGE)

Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/10_3xKO+MyD88-3xFLAG-IRES-Puro+mIRAK4WT-GFP+mIRAK2WT-mScarlet_Analysis.csv.gz"
fwrite(Table, Table_path)

rm(list = ls())
gc()
