library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/49_TKO+hMyD88-3xFLAG-IRES-TagBFP2_hIRAK4WT_hIRAK1WT/20250325/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/49_TKO+hMyD88-3xFLAG-IRES-TagBFP2_hIRAK4WT_hIRAK1WT/20250327/Output/Analysis.csv.gz")
Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/49_TKO+hMyD88-3xFLAG-IRES-TagBFP2_hIRAK4WT_hIRAK1WT/20250401/Output/Analysis.csv.gz")

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
  rep("20250325 plate1_well5B_40nM_cl411_hIRAK4WT_hIRAK1WT_001", 1),
  rep("20250327 plate1_well7C_40nM_cl411_hIRAK4WT_hIRAK1WT_001", 7),
  rep("20250327 plate1_well8D_40nM_cl411_hIRAK4WT_hIRAK1WT_001", 8),
  rep("20250401 plate2_well3F_40nM_cl411_hIRAK4WT_hIRAK1WT_001", 9)
)

# # Define the corresponding CELL values
cell_values <- c(
  7,
  7, 9, 14, 15, 17, 22, 25, 
  2, 5, 9, 12, 14, 15, 17, 19,
  6, 9, 13, 16, 25, 26, 28, 29, 30
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
    SHORT_LABEL = "hIRAK4-WT"
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

Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/43_3xKO+hMyD88-3xFLAG-IRES-TagBFP2+hIRAK4WT-GFP+hIRAK1WT-mScarlet_Analysis.csv.gz"
fwrite(Table, Table_path)


# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()