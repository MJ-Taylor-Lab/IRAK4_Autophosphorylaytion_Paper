# Load required packages
library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

# # Load data from three CSV files into separate tables
# Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/cobra/Analysis Output/20221022/23_IRAK4WT_IRAK1DD/20231130/Output/Analysis.csv.gz")
# Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/cobra/Analysis Output/20221022/23_IRAK4WT_IRAK1DD/20231212/Output/Analysis.csv.gz")
# Table3 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/cobra/Analysis Output/20221022/23_IRAK4WT_IRAK1DD/20231212_2/Output/Analysis.csv.gz")
Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/cobra/Analysis Output/20221022/23_IRAK4WT_IRAK1DD/20240618_1/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/cobra/Analysis Output/20221022/23_IRAK4WT_IRAK1DD/20240618_2/Output/Analysis.csv.gz")

# Combine the three tables into a single table
Table <- rbind(
  Table1,
  Table2
  # # Table3,
  # Table4,
  # Table5
)

# Checking the max, mean, and median values to test if calibration failed --------
max(Table1$MAX_NORMALIZED_INTENSITY)
median(Table1$MAX_NORMALIZED_INTENSITY)
mean(Table1$MAX_NORMALIZED_INTENSITY)

max(Table2$MAX_NORMALIZED_INTENSITY)
median(Table2$MAX_NORMALIZED_INTENSITY)
mean(Table2$MAX_NORMALIZED_INTENSITY)

# max(Table3$MAX_NORMALIZED_INTENSITY)
# median(Table3$MAX_NORMALIZED_INTENSITY)
# mean(Table3$MAX_NORMALIZED_INTENSITY)
# 
# max(Table4$MAX_NORMALIZED_INTENSITY)
# median(Table4$MAX_NORMALIZED_INTENSITY)
# mean(Table4$MAX_NORMALIZED_INTENSITY)
# 
# max(Table5$MAX_NORMALIZED_INTENSITY)
# median(Table5$MAX_NORMALIZED_INTENSITY)
# mean(Table5$MAX_NORMALIZED_INTENSITY)

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
  rep("20240618 plate01_well2D_10nM_cl_D11_IRAK4WT_IRAK1DD_001", 11), #analyse all images again
  rep("20240618 plate01_well2D_10nM_cl_D11_IRAK4WT_IRAK1DD_002", 10),  
  rep("20240618 plate01_well3C_10nM_cl_D11_IRAK4WT_IRAK1DD_001", 23),  
  rep("20240618 plate01_well3C_10nM_cl_D11_IRAK4WT_IRAK1DD_002", 19)
)

# Define the corresponding CELL values
cell_values <- c(
  2, 5, 6, 14, 15, 16, 18, 20, 21, 22, 24,
  1, 2, 3, 5, 6, 7, 8, 9, 10, 11,
  1, 2, 3, 4, 6, 7, 8, 9, 10, 13, 15, 17, 19, 20, 21, 22, 23, 25, 26, 27, 28, 30, 31,
  1, 2, 3, 5, 6, 7, 8, 10, 11, 12, 14, 17, 18, 20, 22, 24, 25, 26, 28
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

unique(Table$IMAGE)

Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/05_3xKO+MyD88-3xFLAG-IRES-Puro+mIRAK4WT-GFP+mIRAK1DD-mScarlet_Analysis.csv.gz"
fwrite(Table, Table_path)

rm(list = ls())
gc()