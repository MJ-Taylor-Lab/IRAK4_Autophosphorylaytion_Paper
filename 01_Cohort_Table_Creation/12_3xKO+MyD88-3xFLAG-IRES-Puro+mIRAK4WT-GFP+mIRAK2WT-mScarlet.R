library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/41_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK3WT/20250225/Output/Analysis.csv.gz")
# Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/41_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK3WT/20250228/Output/Analysis.csv.gz")
Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/41_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK3WT/20250304/Output/Analysis.csv.gz")
Table4 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/41_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK3WT/20250311/Output/Analysis.csv.gz")
Table5 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/41_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK3WT/20250318/Output/Analysis.csv.gz")

Table <- rbind(
  Table1,
  # Table2,
  Table3,
  Table4,
  Table5
)
# Checking the max mean and median values to test if calibration failed --------
max(Table1$MAX_NORMALIZED_INTENSITY)
median(Table1$MAX_NORMALIZED_INTENSITY)
mean(Table1$MAX_NORMALIZED_INTENSITY)

# max(Table2$MAX_NORMALIZED_INTENSITY)
# median(Table2$MAX_NORMALIZED_INTENSITY)
# mean(Table2$MAX_NORMALIZED_INTENSITY)

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
  # Table2,
  Table3,
  Table4,
  Table5
)


# # Manually curating cells -------------------------------------------------
# Define the unique IMAGE values with their corresponding counts
image_values <- c(
  rep("20250225 plate1_well2B_20nM_cl611_IRAK4WT_IRAK3WT_001", 8),
  rep("20250225 plate1_well2D_20nM_cl611_IRAK4WT_IRAK3WT_001", 6),
  # rep("20250228 plate1_well5D_20nM_cl573_IRAK4WT_IRAK3WT_001", 3),
  rep("20230304 plate1_well2C_20nM_cl573_IRAK4WT_IRAK3WT_001", 6),
  rep("20230304 plate1_well2D_20nM_cl573_IRAK4WT_IRAK3WT_001", 6),
  rep("20250311 plate1_well5B_20nM_cl573_IRAK4WT_IRAK3WT_001", 4),
  rep("20250311 plate1_well5D_20nM_cl573_IRAK4WT_IRAK3WT_001", 6),
  rep("20250318 plate1_well3C_20nM_cl573_IRAK4WT_IRAK3WT_001", 1)
)
# Define the corresponding CELL values
cell_values <- c(
  3, 7, 9, 12, 15, 20, 21, 23,
  1, 3, 8, 9, 10, 22,
  # 1, 3, 8,
  2, 3, 4, 7, 11, 13,
  5, 7, 22, 23, 26, 27,
  1, 10, 12, 19, 
  5, 6, 7, 15, 17, 22,
  12
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
    SHORT_LABEL = "IRAK4 WT + IRAK3 WT"
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


Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4_Autophosphorylaytion_Paper_Rewrite/00_Myddosomal_internal_phosphorylation_cohort_table/48_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK3WT.csv.gz"
fwrite(Table, Table_path)

rm(list = ls())
gc()
