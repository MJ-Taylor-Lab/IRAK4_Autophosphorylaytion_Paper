library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/40_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK2DD/20250207/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/40_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK2DD/20250211/Output/Analysis.csv.gz")

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
  rep("20250207 plate1_well4C_20nM_cl610_IRAK4WT_IRAK2DD_001", 12),
  rep("20250207 plate1_well5C_20nM_cl610_IRAK4WT_IRAK2DD_001", 16),
  rep("20250211 plate1_well7C_20nM_cl610_IRAK4WT_IRAK2DD_001", 19),
  rep("20250211 plate1_well8C_20nM_cl610_IRAK4WT_IRAK2DD_001", 12),
  rep("20250211 plate2_well7G_20nM_cl610_IRAK4WT_IRAK2DD_001", 9),
  rep("20250211 plate2_well8G_20nM_cl610_IRAK4WT_IRAK2DD_001", 19)
)

# Define the corresponding CELL values
cell_values <- c(
  2, 7, 13, 14, 15, 19, 22, 23, 24, 25, 28, 34,
  1, 2, 5, 6, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 20, 24,
  2, 3, 4, 5, 6, 9, 10, 11, 13, 14, 15, 17, 18, 21, 23, 24, 25, 26, 27,
  1, 2, 3, 6, 10, 11, 12, 13, 14, 16, 19, 20,
  1, 2, 3, 9, 11, 12, 13, 15, 17,
  1, 3, 4, 6, 8, 9, 10, 11, 15, 16, 19, 20, 24, 25, 26, 27, 28, 29, 31
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
    SHORT_LABEL = "IRAK4 WT + IRAK2 DD"
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


Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4_Autophosphorylaytion_Paper_Rewrite/00_Myddosomal_internal_phosphorylation_cohort_table/40_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK2DD.csv.gz"
fwrite(Table, Table_path)

rm(list = ls())
gc()
