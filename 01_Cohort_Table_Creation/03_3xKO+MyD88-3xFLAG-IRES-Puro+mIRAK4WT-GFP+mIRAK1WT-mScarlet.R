library(pacman)
pacman::p_load(dplyr, tidyr, data.table)

Table1 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/cobra/Analysis Output/20230727_IRAK4PhosphoPaper/20_IRAK4WT_IRAK1WT/20231212/Output/Analysis.csv.gz")
Table2 <- fread("/Volumes/TAYLOR-LAB/Niranjan/04 Image Analysis/cobra/Analysis Output/20230727_IRAK4PhosphoPaper/20_IRAK4WT_IRAK1WT/20231219/Output/Analysis.csv.gz")


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


# Selecting Cells Manually ----------------------------------------------------------
# Define the unique IMAGE values with their corresponding counts
image_values <- c(
  rep("20231212 plate01_well2F_10nM_cl_D11_IRAK4WT_IRAK1WT_001", 15),
  rep("20231212 plate02_well2F_10nM_cl_D11_IRAK4WT_IRAK1WT_001", 12),
  rep("20231219 plate01_well3B_10nM_cl_D11_IRAK4WT_IRAK1WT_002", 11),
  rep("20231219 plate02_well2G_10nM_cl_D11_IRAK4WT_IRAK1WT_001", 5)
)

# Define the corresponding CELL values
cell_values <- c(
  1, 2, 3, 4, 5, 9, 11, 14, 15, 16, 18, 19, 20, 21, 22,
  1, 2, 3, 4, 8, 10, 11, 15, 16, 19, 22, 23,
  1, 2, 3, 4, 6, 10, 11, 12, 17, 18, 19,
  7, 12, 13, 17, 19
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
    SHORT_LABEL = "IRAK4 WT + IRAK1 WT"
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

Table_path <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/00_Cohort_table/03_3xKO+MyD88-3xFLAG-IRES-Puro+mIRAK4WT-GFP+mIRAK1WT-mScarlet_Analysis.csv.gz"
fwrite(Table, Table_path)

rm(list = ls())
gc()
