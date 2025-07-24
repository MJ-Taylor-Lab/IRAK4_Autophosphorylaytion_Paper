# Load required packages
library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

# Load data from two CSV files into separate tables
Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/cobra/Analysis Output/20230727_IRAK4PhosphoPaper/22_IRAK4DD_IRAK1WT/20231212/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/cobra/Analysis Output/20230727_IRAK4PhosphoPaper/22_IRAK4DD_IRAK1WT/20231219/Output/Analysis.csv.gz")

# Combine the two tables into a single table
Table <- rbind(
  Table1,
  Table2
)

# Checking the max, mean, and median values to test if calibration failed --------
max(Table1$MAX_NORMALIZED_INTENSITY)
median(Table1$MAX_NORMALIZED_INTENSITY)
mean(Table1$MAX_NORMALIZED_INTENSITY)

max(Table2$MAX_NORMALIZED_INTENSITY)
median(Table2$MAX_NORMALIZED_INTENSITY)
mean(Table2$MAX_NORMALIZED_INTENSITY)

# Print unique values of COHORT, PROTEIN, and IMAGE columns
unique(Table$COHORT)
unique(Table$PROTEIN)
unique(Table$IMAGE)

# Table Cleanup and save -----------------------------------------------------------
rm(
  Table1,
  Table2
)

# Selecting Cells Manually ----------------------------------------------------------
# Define the unique IMAGE values with their corresponding counts
image_values <- c(
  rep("20231212 plate01_well3G_10nM_cl_D11_IRAK4DD_IRAK1WT_001", 18),
  rep("20231219 plate01_well2B_10nM_cl_D11_IRAK4DD_IRAK1WT_001", 7),
  rep("20231219 plate01_well3C_10nM_cl_D11_IRAK4DD_IRAK1WT_001", 4)
)

# Define the corresponding CELL values
cell_values <- c(
  1, 2, 3, 4, 6, 7, 9, 11, 13, 14, 15, 17, 18, 19, 23, 24, 28, 30,
  1, 2, 5, 7, 10, 11, 13,
  3, 4, 9, 10
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
    SHORT_LABEL = "IRAK4 DD + IRAK1 WT"
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


# Define an exclusion list for incorrectly identified cells ---------------
ManualExclusionList_1 <- data_frame(
  UNIVERSAL_CELL_ID = "20231212 plate01_well3G_10nM_cl_D11_IRAK4DD_IRAK1WT_001_",
  CELL = c(5,8,10,12,16,20,21,22,25,26,27,29)
) 
ManualExclusionList_1 <- ManualExclusionList_1 %>% 
  mutate(
    UNIVERSAL_CELL_ID = paste0(UNIVERSAL_CELL_ID, CELL)
  )

ManualExclusionList_2 <- data_frame(
  UNIVERSAL_CELL_ID = "20231219 plate01_well2B_10nM_cl_D11_IRAK4DD_IRAK1WT_001_",
  CELL = c(3,4,6,8,9,12,14)
) 
ManualExclusionList_2 <- ManualExclusionList_2 %>% 
  mutate(
    UNIVERSAL_CELL_ID = paste0(UNIVERSAL_CELL_ID, CELL)
  )

ManualExclusionList_3 <- data_frame(
  UNIVERSAL_CELL_ID = "20231219 plate01_well3C_10nM_cl_D11_IRAK4DD_IRAK1WT_001_",
  CELL = c(1,2,5,6,7,8,11)
) 
ManualExclusionList_3 <- ManualExclusionList_3 %>% 
  mutate(
    UNIVERSAL_CELL_ID = paste0(UNIVERSAL_CELL_ID, CELL)
  )

# Combine the exclusion lists
ManualExclusionList <- rbind(
  ManualExclusionList_1,
  ManualExclusionList_2,
  ManualExclusionList_3
)

# Remove unnecessary objects
rm(
  ManualExclusionList_1,
  ManualExclusionList_2,
  ManualExclusionList_3
)


# Remove incorrectly identified cells based on the exclusion list ---------
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
Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/04_3xKO+MyD88-3xFLAG-IRES-Puro+mIRAK4DD-GFP+mIRAK1WT-mScarlet_Analysis.csv.gz"
fwrite(Table, Table_path)

# Clean up
rm(list = ls())
