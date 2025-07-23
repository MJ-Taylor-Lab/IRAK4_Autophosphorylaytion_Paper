library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/50_TKO+hMyD88-3xFLAG-IRES-TagBFP2_hIRAK4-D329A_hIRAK1WT/20250325_01/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/50_TKO+hMyD88-3xFLAG-IRES-TagBFP2_hIRAK4-D329A_hIRAK1WT/20250325_02/Output/Analysis.csv.gz")
Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/50_TKO+hMyD88-3xFLAG-IRES-TagBFP2_hIRAK4-D329A_hIRAK1WT/20250327/Output/Analysis.csv.gz")
Table4 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/raven/Analysis_Output/IRAK4_Kinase_Paper_InProgress/50_TKO+hMyD88-3xFLAG-IRES-TagBFP2_hIRAK4-D329A_hIRAK1WT/20250401/Output/Analysis.csv.gz")

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


# # Manually curating cells -------------------------------------------------
# Define the unique IMAGE values with their corresponding counts
image_values <- c(
  rep("20250325 plate1_well5G_40nM_cl621_hIRAK4D329A_hIRAK1WT_001", 11),
  rep("20250325 plate1_well6B_40nM_cl621_hIRAK4D329A_hIRAK1WT_001", 7),
  rep("20250327 plate1_well7B_40nM_cl621_hIRAK4D329A_hIRAK1WT_", 12),
  rep("20250327 plate1_well8B_40nM_cl621_hIRAK4D329A_hIRAK1WT_", 11),
  rep("20250401 plate2_well3G_40nM_cl621_hIRAK4D329A_hIRAK1WT_", 10)
)
# # Define the corresponding CELL values
cell_values <- c(
  1, 2, 4, 5, 7, 12, 21, 22, 23, 24, 25,
  3, 7, 9, 25, 29, 33, 34,
  3, 5, 12, 13, 14, 18, 21, 22, 25, 26, 30, 32, 
  1, 6, 8, 10, 14, 15, 18, 19, 21, 22, 23,
  1, 2, 5, 6, 7, 11, 13, 14, 17, 18
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
    SHORT_LABEL = "hIRAK4-D329A"
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


Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4_Autophosphorylaytion_Paper_Rewrite/00_Myddosomal_internal_phosphorylation_cohort_table/57_TKO+hMyD88-3xFLAG-IRES-Puro_hIRAK4-D329A_hIRAK1WT.csv.gz"
fwrite(Table, Table_path)

rm(list = ls())
gc()
