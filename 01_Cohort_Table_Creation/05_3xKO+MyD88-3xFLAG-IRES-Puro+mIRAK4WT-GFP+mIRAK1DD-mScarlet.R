# Load required packages
library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

# # Load data from three CSV files into separate tables
# Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/Analysis Output/20221022/23_IRAK4WT_IRAK1DD/20231130/Output/Analysis.csv.gz")
# Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/Analysis Output/20221022/23_IRAK4WT_IRAK1DD/20231212/Output/Analysis.csv.gz")
# Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/Analysis Output/20221022/23_IRAK4WT_IRAK1DD/20231212_2/Output/Analysis.csv.gz")
Table4 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/Analysis Output/20221022/23_IRAK4WT_IRAK1DD/20240618_1/Output/Analysis.csv.gz")
Table5 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/Analysis Output/20221022/23_IRAK4WT_IRAK1DD/20240618_2/Output/Analysis.csv.gz")

# Combine the three tables into a single table
Table <- rbind(
  # Table1,
  # Table2,
  # Table3,
  Table4,
  Table5
)

# Checking the max, mean, and median values to test if calibration failed --------
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

# Print unique values of COHORT, PROTEIN, and IMAGE columns
unique(Table$COHORT)
unique(Table$PROTEIN)
unique(Table$IMAGE)

# Remove unnecessary objects
rm(
  Table1,
  Table2,
  Table3,
  Table4,
  Table5
)

# Convert the combined table to data.table class
Table <- Table %>%  as.data.table()

# Filter data for a specific cohort and mutate column values
Table <- Table %>% 
  filter(
    COHORT == "TKO+MyD88-3xFLAG+IRAK4WT+IRAK1DD",
    IMAGE != "20231130 plate01_well7C_10nM_cl_D11_IRAK4WT_IRAK1DD_001"
  ) %>% 
  mutate(
    SHORT_LABEL = "IRAK4 WT + IRAK1 DD"
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


# Exclusion List ----------------------------------------------------------
# Excluding data from cells that were wrongly identified
ManualExclusionList_1 <- data_frame(
  UNIVERSAL_CELL_ID = "20231130 plate02_well6C_10nM_cl_D11_IRAK4WT_IRAK1DD_002",
  CELL = c(1,2,3,4,5,6,7,8,9,10,11,13,14,15,17,18,19,20)
) 
ManualExclusionList_1 <- ManualExclusionList_1 %>% 
  mutate(
    UNIVERSAL_CELL_ID = paste0(UNIVERSAL_CELL_ID, CELL)
  )


ManualExclusionList_2 <- data_frame(
  UNIVERSAL_CELL_ID = "20231130 plate02_well6C_10nM_cl_D11_IRAK4WT_IRAK1DD_003",
  CELL = c(1,3,4,5,7,10,11,12,13,14,16,18,19,21,24,26,28)
)
ManualExclusionList_2 <- ManualExclusionList_2 %>% 
  mutate(
    UNIVERSAL_CELL_ID = paste0(UNIVERSAL_CELL_ID, CELL)
  )

ManualExclusionList_3 <- data_frame(
  UNIVERSAL_CELL_ID = "20231130 plate02_well6C_10nM_cl_D11_IRAK4WT_IRAK1DD_004",
  CELL = c(4,7)
)
ManualExclusionList_3 <- ManualExclusionList_3 %>% 
  mutate(
    UNIVERSAL_CELL_ID = paste0(UNIVERSAL_CELL_ID, CELL)
  )

ManualExclusionList_4 <- data_frame(
  UNIVERSAL_CELL_ID = "20231212 plate01_well2G_10nM_cl_D11_IRAK4WT_IRAK1DD_001",
  CELL = c(1,3,4,5,8,10,11,12,14,15,16,18,19,23,24)
)
ManualExclusionList_4 <- ManualExclusionList_4 %>% 
  mutate(
    UNIVERSAL_CELL_ID = paste0(UNIVERSAL_CELL_ID, CELL)
  )

ManualExclusionList_5 <- data_frame(
  UNIVERSAL_CELL_ID = "20231212 plate01_well2G_10nM_cl_D11_IRAK4WT_IRAK1DD_002",
  CELL = c(4,7,11,14,15,16,17,19,21,23,25,26)
) 
ManualExclusionList_5 <- ManualExclusionList_5 %>% 
  mutate(
    UNIVERSAL_CELL_ID = paste0(UNIVERSAL_CELL_ID, CELL)
  )


ManualExclusionList <- rbind(
  ManualExclusionList_1,
  ManualExclusionList_2,
  ManualExclusionList_3,
  ManualExclusionList_4,
  ManualExclusionList_5
)

rm(
  ManualExclusionList_1,
  ManualExclusionList_2,
  ManualExclusionList_3,
  ManualExclusionList_4,
  ManualExclusionList_5
)


# Removing incorrectly identified cells -----------------------------------
Table <- Table %>% 
  mutate(
    UNIVERSAL_CELL_ID = paste0(IMAGE, "_", CELL),
    UNIVERSAL_CELL_ID = UNIVERSAL_CELL_ID %in% ManualExclusionList$UNIVERSAL_CELL_ID
  ) %>% 
  filter(
    UNIVERSAL_CELL_ID == FALSE
  ) %>% 
  select(
    -UNIVERSAL_CELL_ID
  )

# Save the processed table to a CSV file
Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4_Autophosphorylaytion_Paper_Rewrite/00_Myddosomal_internal_phosphorylation_cohort_table/05_3xKO_MyD88-IRES-Puro_IRAK4WT_IRAK1DD_Compiled_Essential.csv.gz"
fwrite(Table, Table_path)

# Clean up
rm(list = ls())
