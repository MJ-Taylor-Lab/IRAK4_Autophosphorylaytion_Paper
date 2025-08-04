# Load necessary libraries using pacman for easy package management
library(pacman)
pacman::p_load(dplyr, data.table, ggbeeswarm, ggplot2)

# Read data from a CSV file into a data table
Table <- fread("/Users/u_niranjan/Desktop/Image_Analysis_Grid/IRAK3KinaseDomainOnly/TrackMate_Information.csv")

# Define the folder path where plots will be saved
Plot_Save_Folder <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots"

# Create directories to save the plot if they don't already exist
Plot_Directory_Save_Path <- file.path(Plot_Save_Folder, "04_Figure 4")
if (!file.exists(Plot_Directory_Save_Path)) {
  dir.create(Plot_Directory_Save_Path)
}

# Create directories to save the plot if they don't already exist
Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "07_Figure4J")
if (!file.exists(Plot_Directory_Save_Path)) {
  dir.create(Plot_Directory_Save_Path)
}

Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "01_IRAK1KinaseDomainOnly_Violin_Plot")
if (!file.exists(Plot_Directory_Save_Path)) {
  dir.create(Plot_Directory_Save_Path)
}

# Process the data: filter, group by IMAGE, CELL_NUMBER, and COHORT, then summarise and mutate
Table <- Table %>% 
  # filter(
  #   LIFETIME > 0
  # ) %>%
  group_by(
    IMAGE,
    REPLICATE,
    CELL_NUMBER,
    COHORT
  ) %>% 
  summarise(
    NUMBER_OF_EVENTS = mean(NUMBER_OF_EVENTS)  # Calculate mean number of events
  ) %>% 
  arrange(
    IMAGE
  )  %>% 
  mutate(
    # Assign labels based on COHORT values
    SHORT_LABEL = case_when(
      COHORT == "DMSO" ~ "DMSO",
      COHORT == "PF06650833_20uM" ~ "IRAK4 Kinase\nInhibitor 20 µM"
    )
  ) %>% 
  group_by(
    IMAGE
  ) %>% 
  mutate(
    MAX_CELL_NUMBER_BY_FOV = max(CELL_NUMBER) # Get the number of cells per field of view
  ) %>% 
  filter(
    IMAGE != "20241210 grid04_20nM_cl501_IRAK3KinaseDomain_inhib20uM_005",
    IMAGE != "20241210 grid04_20nM_cl501_IRAK3KinaseDomain_inhib20uM_006",
    IMAGE != "20241204 grid01_20nM_cl501_IRAK3KinaseDomain_DMSO_001",
    IMAGE != "20241204 grid01_20nM_cl501_IRAK3KinaseDomain_DMSO_003",
    IMAGE != "20241127 grid02_20nM_cl501_IRAK3KinaseDomain_inhib20uM_003",
    IMAGE != "20241127 grid02_20nM_cl501_IRAK3KinaseDomain_inhib20uM_006"
  )

# Set the order of factors for SHORT_LABEL
Table$SHORT_LABEL <- factor(
  Table$SHORT_LABEL,
  levels = c("IRAK4 Kinase\nInhibitor 20 µM", "DMSO")
)


# Further processing to group and summarise data by IMAGE and REPLICATE
Table_byREPLICATE <- Table %>% 
  group_by(
    IMAGE,
    SHORT_LABEL,
    REPLICATE,
    COHORT
  ) %>% 
  summarise(
    MAX_CELL_NUMBER_BY_FOV = max(MAX_CELL_NUMBER_BY_FOV),
    MEAN_NUMBER_OF_EVENTS = mean(NUMBER_OF_EVENTS)  # Calculate mean number of events by Image 
  ) %>% 
  group_by(
    SHORT_LABEL,
    COHORT,
    REPLICATE
  ) %>% 
  summarise(
    MAX_CELL_NUMBER_BY_REPLICATE = sum(MAX_CELL_NUMBER_BY_FOV),
    MEAN_NUMBER_OF_EVENTS = mean(MEAN_NUMBER_OF_EVENTS)  # Calculate mean number of events by Replicate 
  )  %>% 
  mutate(
    # Assign colors based on SHORT_LABEL
    COLOR = case_when(
      SHORT_LABEL == "DMSO" ~ "#FED09E",
      SHORT_LABEL == "IRAK4 Kinase\nInhibitor 20 µM" ~ "#CC7B16"
    )
  ) %>% 
  group_by(
    SHORT_LABEL
  ) %>% 
  mutate(
    REPLICATE_NUMBER = n_distinct(REPLICATE)
  )

# Calculate p-value for the difference in MEAN_NUMBER_OF_EVENTS between cohorts
p_value_Result <- Table_byREPLICATE
p_value_Result_1 <- wilcox.test(data = p_value_Result, MEAN_NUMBER_OF_EVENTS ~ SHORT_LABEL)$p.value
p_value_Result_1 <- signif(p_value_Result_1, digits = 3)
rm(p_value_Result)

# P value calculation
p_value_Result <- Table_byREPLICATE

p_value_Result_1 <- wilcox.test(
  data = p_value_Result,
  MEAN_NUMBER_OF_EVENTS ~ SHORT_LABEL
)$p.value

p_value_Result_1 <- signif(p_value_Result_1, digits = 3)

rm(p_value_Result) 


# Further processing to arrange and group data by COHORT, then summarise
Table_byCOHORT <- Table_byREPLICATE %>% 
  group_by(
    COHORT,
    SHORT_LABEL
  ) %>%
  summarise(
    SEM = sd(MEAN_NUMBER_OF_EVENTS), # Calculate Stanard Error of mean
    REPLICATE_NUMBER = mean(REPLICATE_NUMBER),
    SEM = SEM/ REPLICATE_NUMBER,
    MEDIAN_NUMBER_OF_EVENTS = median(MEAN_NUMBER_OF_EVENTS),  # Calculate median number of events
    MEAN_NUMBER_OF_EVENTS = mean(MEAN_NUMBER_OF_EVENTS),  # Calculate mean number of events
    MEAN_PLUS_SEM = MEAN_NUMBER_OF_EVENTS + SEM,
    MEAN_MINUS_SEM = case_when(
      MEAN_NUMBER_OF_EVENTS - SEM < 0 ~ 0,
      TRUE ~ MEAN_NUMBER_OF_EVENTS - SEM
    )
  ) %>% 
  mutate(
    # Assign colors based on SHORT_LABEL
    COLOR = case_when(
      SHORT_LABEL == "DMSO" ~ "#FED09E",
      SHORT_LABEL == "IRAK4 Kinase\nInhibitor 20 µM" ~ "#CC7B16"
    )
  )


# Write csv ---------------------------------------------------------------
Source_Data_Path <- paste0(basename(Plot_Directory_Save_Path), ".csv")
Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Source_Data_Path)
write.csv(Table, Plot_Save_Path)


# GGPLOT ------------------------------------------------------------------
# Create a beeswarm plot
Plot <- ggplot(data = Table) +
  geom_violin(
    scale = "width",
    linewidth = 0.2,
    aes(
      y = SHORT_LABEL, 
      x = NUMBER_OF_EVENTS, 
      color = SHORT_LABEL
    ),
    width = 0.5,
    color = "grey",
    fill = "grey"
  ) +
  geom_errorbar(
    data = Table_byCOHORT,
    aes(
      y = SHORT_LABEL, 
      x = MEAN_NUMBER_OF_EVENTS, 
      xmin = MEAN_MINUS_SEM, 
      xmax = MEAN_PLUS_SEM
    ),
    width = 0.1,
    size = 0.2
  ) + 
  geom_crossbar(
    data = Table_byCOHORT,
    aes(
      y = SHORT_LABEL, 
      x = MEAN_NUMBER_OF_EVENTS, 
      xmin = MEAN_NUMBER_OF_EVENTS, 
      xmax = MEAN_NUMBER_OF_EVENTS
    ),
    width = 0.1,
    size = 0.1,
    color = "black"
  ) +
  geom_jitter(
    data = Table_byREPLICATE,
    aes(
      y = SHORT_LABEL, 
      x = MEAN_NUMBER_OF_EVENTS
    ),
    color = "black",
    fill = Table_byREPLICATE$COLOR,
    shape = 21,
    size = 2,
    stroke = 0.2,
    height = 0.1
  ) +
  labs(
    x = "# Number of Events"
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(-0.1, 150),
    position = "top"
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 7),
    axis.text = element_text(size = 6),
    legend.position = "none",
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.background = element_blank(),
    axis.line.x.top = element_line(color = "black"),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )

# Create labels for p-value and cell information
Cell_Information = "Y"
if (Cell_Information == "Y") {
  Cell_Number_Label <- c("Total Number of Cells treated with DMSO = ")
  Temp <- Table_byREPLICATE %>% filter(SHORT_LABEL == "DMSO")
  Cell_Number_Label <- paste0(Cell_Number_Label, " ", Temp$MAX_CELL_NUMBER_BY_REPLICATE[1])
  for (index in 2:nrow(Temp)) {
    Cell_Number_Label <- paste0(Cell_Number_Label, " /", Temp$MAX_CELL_NUMBER_BY_REPLICATE[index])
  }
  Cell_Number_Label <- paste0(Cell_Number_Label, " ==> Total =", sum(Temp$MAX_CELL_NUMBER_BY_REPLICATE), "\nTotal Number of Cells treated with inhibitor = ")
  Temp <- Table_byREPLICATE %>% filter(SHORT_LABEL != "DMSO")
  Cell_Number_Label <- paste0(Cell_Number_Label, " ", Temp$MAX_CELL_NUMBER_BY_REPLICATE[1])
  for (index in 2:nrow(Temp)) {
    Cell_Number_Label <- paste0(Cell_Number_Label, " /", Temp$MAX_CELL_NUMBER_BY_REPLICATE[index])
  }
  Cell_Number_Label <- paste0(Cell_Number_Label, " ==> Total =", sum(Temp$MAX_CELL_NUMBER_BY_REPLICATE))
}
rm(Temp, Cell_Information)

# Add p-value and cell information to the plot
Plot_pvalue <- Plot +
  geom_segment(
    aes(
      y = 1, 
      yend = 2, 
      x = 125,
      xend = 125
    ),
    size = 0.05,
    color = "black"
  ) +
  geom_segment(
    aes(
      y = 1, 
      yend = 1, 
      x = 120, 
      xend = 125
    ),
    size = 0.05,
    color = "black"
  ) +
  geom_segment(
    aes(
      y = 2, 
      yend = 2, 
      x = 120, 
      xend = 125
    ),
    size = 0.05,
    color = "black"
  ) +
  annotate(
    "text",
    y = 1.5,
    x = 128,
    label = p_value_Result_1,
    size = 3
  ) +
  annotate(
    "text",
    y = 1.5,
    x = 60,
    label = Cell_Number_Label,
    size = 3
  )

# Save the plot with p-value to a file
Plot_Save_Path <- file.path(Plot_Directory_Save_Path, "01_Violin Plot of Number of Events with p value.pdf")
ggsave(Plot_Save_Path, plot = Plot_pvalue)

# Save the plot for publication
Plot_Save_Path <- file.path(Plot_Directory_Save_Path, "02_Violin Plot of Number of Events for publication.pdf")
ggsave(
  Plot_Save_Path, 
  plot = Plot, 
  height = 50, 
  width = 50, 
  units = "mm"
)


# Clean up the environment by removing all objects and running garbage collection
rm(list = ls())
gc()
