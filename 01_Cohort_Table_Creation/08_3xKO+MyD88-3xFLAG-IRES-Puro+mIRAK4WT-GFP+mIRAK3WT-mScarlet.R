library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/37_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK3WT/20250114/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/37_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK3WT/20250128/Output/Analysis.csv.gz")
Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/37_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK3WT/20250131/Output/Analysis.csv.gz")
Table4 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/37_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK3WT/20250204/Output/Analysis.csv.gz")

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
  rep("20250114 plate1_well2D_20nM_cl573_IRAK4WT_IRAK3WT_001", 12),
  rep("20250114 plate1_well3C_20nM_cl573_IRAK4WT_IRAK3WT_001", 13),
  rep("20250128 plate1_well8B_20nM_cl573_IRAK4WT_IRAK3WT_001", 25),
  rep("20250128 plate1_well9C_20nM_cl573_IRAK4WT_IRAK3WT_001", 15),
  rep("20250204 plate1_well2D_20nM_cl573_IRAK4WT_IRAK3WT_001", 6),
  rep("20250204 plate1_well3C_20nM_cl573_IRAK4WT_IRAK3WT_001", 17)
)

# Define the corresponding CELL values
cell_values <- c(
  1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 15, 18,
  2, 7, 8, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 28, 29,
  2, 3, 4, 9, 10, 11, 12, 14, 19, 20, 21, 22, 24, 25, 26,
  4, 5, 8, 12, 13, 14, 
  1, 2, 4, 5, 6, 7, 8, 10, 12, 13, 14, 15, 17, 18, 19, 20, 22
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

unique(Table$IMAGE)

Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/08_3xKO+MyD88-3xFLAG-IRES-Puro+mIRAK4WT-GFP+mIRAK3WT-mScarlet_Analysis.csv.gz"
fwrite(Table, Table_path)

rm(list = ls())
gc()
