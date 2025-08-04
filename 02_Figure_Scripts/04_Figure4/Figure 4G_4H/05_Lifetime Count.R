# Load necessary libraries using pacman for easy package management
library(pacman)
pacman::p_load(dplyr, data.table, ggplot2, lemon, scales)

# Read data from a CSV file into a data table
Table <- fread("/Users/u_niranjan/Desktop/Image_Analysis_Grid/IRAK2KinaseDomainOnly/TrackMate_Information.csv")
# Define the folder path where plots will be saved
Plot_Save_Folder <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots"

# Create directories to save the plot if they don't already exist
Plot_Directory_Save_Path <- file.path(Plot_Save_Folder, "04_Figure 4")
if (!file.exists(Plot_Directory_Save_Path)) {
  dir.create(Plot_Directory_Save_Path)
}

# Create directories to save the plot if they don't already exist
Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "06_Figure4H")
if (!file.exists(Plot_Directory_Save_Path)) {
  dir.create(Plot_Directory_Save_Path)
}

Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "01_IRAK1KinaseDomainOnly_Lifetime_Plot")
if (!file.exists(Plot_Directory_Save_Path)) {
  dir.create(Plot_Directory_Save_Path)
}

# Process the data: group by IMAGE, CELL_NUMBER, and COHORT, then summarise and mutate
Table <- Table %>% 
  filter(
    LIFETIME > 0,
    IMAGE != "20241120 grid01_20nM_cl500_IRAK2KinaseDomain_DMSO_007",
    IMAGE != "20241127 grid03_20nM_cl500_IRAK2KinaseDomain_DMSO_004",
    IMAGE != "20241127 grid03_20nM_cl500_IRAK2KinaseDomain_DMSO_005",
    IMAGE != "20241127 grid03_20nM_cl500_IRAK2KinaseDomain_DMSO_007",
    IMAGE != "20241210 grid01_20nM_cl500_IRAK2KinaseDomain_DMSO_001",
    IMAGE != "20241210 grid01_20nM_cl500_IRAK2KinaseDomain_DMSO_002",
    IMAGE != "20241210 grid01_20nM_cl500_IRAK2KinaseDomain_DMSO_003",
    IMAGE != "20241210 grid01_20nM_cl500_IRAK2KinaseDomain_DMSO_002",
    IMAGE != "20241120 grid02_20nM_cl500_IRAK2KinaseDomain_inhib20uM_004",
    IMAGE != "20241120 grid02_20nM_cl500_IRAK2KinaseDomain_inhib20uM_006",
    IMAGE != "20241127 grid04_20nM_cl500_IRAK2KinaseDomain_inhib20uM_006",
    IMAGE != "20241210 grid02_20nM_cl500_IRAK2KinaseDomain_inhib20uM_003",
    IMAGE != "20241210 grid02_20nM_cl500_IRAK2KinaseDomain_inhib20uM_001",
    IMAGE != "20241120 grid01_20nM_cl500_IRAK2KinaseDomain_DMSO_004",
    IMAGE != "20241127 grid03_20nM_cl500_IRAK2KinaseDomain_DMSO_002",
    IMAGE != "20241127 grid04_20nM_cl500_IRAK2KinaseDomain_inhib20uM_005",
    IMAGE != "20241210 grid02_20nM_cl500_IRAK2KinaseDomain_inhib20uM_004" & CELL_NUMBER != 12,
    IMAGE != "20241127 grid04_20nM_cl500_IRAK2KinaseDomain_inhib20uM_007" & CELL_NUMBER != 13,
    IMAGE != "20241210 grid02_20nM_cl500_IRAK2KinaseDomain_inhib20uM_004" & CELL_NUMBER != 12,
    IMAGE != "20241210 grid02_20nM_cl500_IRAK2KinaseDomain_inhib20uM_007" & CELL_NUMBER != 11,
    IMAGE != "20250530 grid02_20nM_cl500_IRAK2KinaseDomain_inhib20uM_006",
    IMAGE != "20250530 grid02_20nM_cl500_IRAK2KinaseDomain_inhib20uM_007",
    IMAGE != "20250530 grid01_20nM_cl500_IRAK2KinaseDomain_DMSO_006",
    IMAGE != "20250530 grid02_20nM_cl500_IRAK2KinaseDomain_inhib20uM_002" & CELL_NUMBER != 9,
    IMAGE != "20241127 grid04_20nM_cl500_IRAK2KinaseDomain_inhib20uM_003" & CELL_NUMBER != 1,
    IMAGE != "20250530 grid02_20nM_cl500_IRAK2KinaseDomain_inhib20uM_005" & CELL_NUMBER != 3,
    IMAGE != "20250530 grid02_20nM_cl500_IRAK2KinaseDomain_inhib20uM_003" & CELL_NUMBER != 28,
    IMAGE != "20250530 grid02_20nM_cl500_IRAK2KinaseDomain_inhib20uM_004" & CELL_NUMBER != 17,
    IMAGE != "20250530 grid03_20nM_cl500_IRAK2KinaseDomain_inhib20uM_007"
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID,
    COHORT
  ) %>% 
  summarise(
    LIFETIME = mean(LIFETIME)  # Calculate mean number of events
  ) %>% 
  mutate(
    # Assign labels based on COHORT values
    SHORT_LABEL = case_when(
      COHORT == "DMSO" ~ "DMSO",
      COHORT == "PF06650833_20uM" ~ "IRAK4 Kinase\nInhibitor 20 µM"
    )
  )


Table_BIN <- Table %>% 
  group_by(
    COHORT,
    SHORT_LABEL
  ) %>% 
  mutate(
      BIN = (trunc(LIFETIME / 0.2)) * 0.2,  # Bin dwell time into intervals of 3
  ) %>% 
  group_by(
    SHORT_LABEL, 
    COHORT, 
    BIN
  ) %>%
  mutate(
    BIN_COUNT = n()  # Count the number of occurrences in each BIN
  ) %>% 
  distinct(
    BIN_COUNT, 
    BIN
  ) %>% 
  arrange(
    SHORT_LABEL,
    BIN
  ) %>% 
  filter(
    BIN < 100
  ) %>% 
  as.data.table()
  

Table_BIN$SHORT_LABEL <- factor(
  Table_BIN$SHORT_LABEL,
  levels = c("DMSO", "IRAK4 Kinase\nInhibitor 20 µM")
)

# Write csv ---------------------------------------------------------------
Source_Data_Path <- paste0(basename(Plot_Directory_Save_Path), ".csv")
Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Source_Data_Path)
write.csv(Table, Plot_Save_Path)


# GGPLOT
LOWER_LIMIT = 0.1           #Axis Lower Limit
UPPER_LIMIT = 1.6       #Axis Upper Limit
AXIS_BREAK_SEQ = 0.6        #X axis tick marks
FACET_ROW_NUMBERS = 2     #Number of facet rows
X_LABEL = "Lifetime (s)"
Y_LABEL = "# of Events"

Plot <- ggplot(
  data = Table
) +
  geom_bar(
    data = Table_BIN,
    aes(
      y = BIN_COUNT, 
      x = BIN,
      fill = SHORT_LABEL
    ),
    width = 0.2,
    stat = "identity",
    color = "black",
    linewidth = 0.1
  ) +
  scale_fill_manual(
    values = c("#FED09E", "#CC7B16")
  ) +
  facet_rep_wrap(
    ~SHORT_LABEL,
    scales = "fixed",
    repeat.tick.labels = TRUE,
    nrow = FACET_ROW_NUMBERS
  ) +
  scale_x_continuous(
    limits =  c(LOWER_LIMIT,UPPER_LIMIT), ## Sets x axis limts
    breaks = seq(0.3,UPPER_LIMIT, by = AXIS_BREAK_SEQ),
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = pretty_breaks(n = 3)
  ) +
  labs(
    x = X_LABEL,
    y = Y_LABEL
  ) +
  theme_classic()
  
  
Plot +
  theme(
    axis.title.y = element_text(size = 40),  # Style Y-axis title
    axis.text = element_text(size = 40),  # Style axis text
    legend.position = "none",  # Remove legend
    strip.text.x = element_blank()
  )


# Save the plot to a file
Plot_Save_Path <- file.path(Plot_Directory_Save_Path, "02_Lifetime Plot of IRAK1 Kinase Domain only for presentation.pdf")
ggsave(
  Plot_Save_Path,
  plot = Plot,
  height = 3*5,
  width = 5*5
)


Plot +
  theme(
    axis.title = element_text(size = 7),  # Style Y-axis title
    axis.text = element_text(size = 7),  # Style axis text
    legend.position = "none",  # Remove legend
    strip.text.x = element_blank()
  )


# Save the plot to a file
Plot_Save_Path <- file.path(Plot_Directory_Save_Path, "02_Lifetime Plot of IRAK1 Kinase Domain only.pdf")
ggsave(
  Plot_Save_Path,
  plot = last_plot(),
  height = 45,
  width = 38,
  units = "mm"
)


# # Clean up the environment by removing all objects and running garbage collection
rm(list = ls())
gc()