library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/51_TKO+hMyD88-3xFLAG-IRES-TagBFP2_hIRAK4-D329A-Phosphomimetic_hIRAK1WT/20250325/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/51_TKO+hMyD88-3xFLAG-IRES-TagBFP2_hIRAK4-D329A-Phosphomimetic_hIRAK1WT/20250327/Output/Analysis.csv.gz")
Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/51_TKO+hMyD88-3xFLAG-IRES-TagBFP2_hIRAK4-D329A-Phosphomimetic_hIRAK1WT/20250401/Output/Analysis.csv.gz")

Table <- rbind(
  Table1,
  Table2,
  Table3
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


unique(Table$COHORT)
unique(Table$PROTEIN)
unique(Table$IMAGE)

# Table Cleanup and save -----------------------------------------------------------
rm(
  Table1,
  Table2,
  Table3
)


# # Manually curating cells -------------------------------------------------
# Define the unique IMAGE values with their corresponding counts
image_values <- c(
  rep("20250325 plate1_well4C_40nM_cl567_hIRAK4D329A-PM_hIRAK1WT_001", 7),
  rep("20250325 plate2_well4F_40nM_cl567_hIRAK4D329A-PM_hIRAK1WT_001", 16),
  rep("20250327 plate1_well7D_40nM_cl567_hIRAK4D329A-PM_hIRAK1WT_", 19),
  rep("20250401 plate2_well2F_40nM_cl567_hIRAK4D329A-PM_hIRAK1WT_", 29)
)
# # Define the corresponding CELL values
cell_values <- c(
  10, 11, 13, 17, 18, 24, 27,
  2, 6, 7, 8, 10, 11, 12, 13, 16, 17, 20, 21, 22, 23, 24, 25,
  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20,
  1, 2, 3, 6, 7, 8, 10, 11, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34
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
    SHORT_LABEL = "hIRAK4-D329A-Phosphomimetic"
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

Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/45_3xKO+hMyD88-3xFLAG-IRES-TagBFP2+hIRAK4D329A-Phosphomimetic-GFP+hIRAK1WT-mScarlet_Analysis.csv.gz"
fwrite(Table, Table_path)


# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()