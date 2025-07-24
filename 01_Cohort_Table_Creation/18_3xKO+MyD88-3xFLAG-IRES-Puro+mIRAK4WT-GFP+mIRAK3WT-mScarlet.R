library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/47_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK3WT/20250328/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/47_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK3WT/20250401_01/Output/Analysis.csv.gz")
Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/47_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK3WT/20250401_02/Output/Analysis.csv.gz")
Table4 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/47_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK3WT/20250416_01/Output/Analysis.csv.gz")
Table5 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/47_TKO+MyD88-3xFLAG-IRES-Puro_mIRAK4WT_mIRAK3WT/20250416_02/Output/Analysis.csv.gz")

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
  rep("20250328 plate1_well2B_40mM_cl573_IRAK4WT_IRAK3WT_001", 1),
  rep("20250401 plate1_well5B_40mM_cl573_IRAK4WT_IRAK3WT_001", 9),
  rep("20250401 plate2_well2E_40mM_cl573_IRAK4WT_IRAK3WT_001", 7),
  rep("20250416 plate1_well2D_40mM_cl573_IRAK4WT_IRAK3WT_001", 3),
  rep("20250416 plate1_well3D_40mM_cl573_IRAK4WT_IRAK3WT_001", 4)
)

# # Define the corresponding CELL values
cell_values <- c(
  3, 
  3, 8, 9, 13, 16, 17, 21, 22, 25,
  4, 7, 8, 16, 17, 21, 23,
  1, 3, 6, 
  2, 4, 6, 7 
)
# 
# # Create the Selection_Table data frame
Selection_Table <- data.frame(
  IMAGE = image_values,
  CELL = cell_values
)
# 
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

Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/18_3xKO+MyD88-3xFLAG-IRES-Puro+mIRAK4WT-GFP+mIRAK3WT-mScarlet_Analysis.csv.gz"
fwrite(Table, Table_path)

rm(list = ls())
gc()