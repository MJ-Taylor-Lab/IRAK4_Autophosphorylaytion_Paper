library(pacman)
pacman::p_load(readxl, dplyr, data.table, ggplot2)

# Setting directory paths
Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots_short"

Primary_Data <- "/Users/u_niranjan/Desktop/Cell_Figures/00_New_Figures/Primary_Data.xlsx"

Function_folder <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/03_PlotFromPrimaryData/Functions"

Y_AXIS_LOWER_LIMIT = -0.5
Y_AXIS_UPPER_LIMIT = 130

# Opening Data for relevant plot ----------------------------------------
Cell_Summary_by_Day <- read_excel(Primary_Data, sheet = "Supplementary Figure 3G")

# Creating folder to save data --------------------------------------------
# Creating Folder to save Images
Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "03_Supplementary Figure 3")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

# Creating Folder to save Images
Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "Supplementary Figure_3G")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}


Cell_Summary_by_Cohort <- Cell_Summary_by_Day %>% 
  group_by(
    Treatment_Condition_short
  ) %>%
  summarise(
    Mean_Fold_Change = mean(Fold_Change, na.rm = TRUE),
    SD_Fold_Change = sd(Fold_Change, na.rm = TRUE),
    .groups = 'drop'
  ) %>% 
  mutate(
    # Set YMIN to 0 if the calculation results in a negative value
    YMIN = pmax(0, Mean_Fold_Change - SD_Fold_Change),
    YMAX = Mean_Fold_Change + SD_Fold_Change,
    # Define colors for the bars.
    COLOR = case_when(
      Treatment_Condition_short == "DMSO"  ~ "#FED09E",
      Treatment_Condition_short == "500nM" ~ "#FAA73F",
      TRUE ~ "#CC7B16"
    )
  )

# Convert SHORT_LABEL to a factor with specified levels
Cell_Summary_by_Cohort$Treatment_Condition_short <- factor(
  Cell_Summary_by_Cohort$Treatment_Condition_short,
  levels = c("DMSO", "500nM", "20uM")
)

# -----------------------------------------------------------------------------
# STATISTICAL ANALYSIS: PAIRWISE ONE WAY ANOVA TEST
# -----------------------------------------------------------------------------

# Perform pairwise One way anova test to compare fold changes between treatment groups.
p_value_test_confirmation <- "Y"
if (p_value_test_confirmation == "Y") {
  
  # Get the unique treatment conditions to compare.
  treatments_to_compare <- as.character(unique(Cell_Summary_by_Day$Treatment_Condition_short))
  
  # Generate all unique pairs of treatments.
  treatment_combinations <- combn(sort(treatments_to_compare), 2, simplify = FALSE)
  
  # Define a function to perform the test for a single pair.
  run_one_way_anova_test <- function(pair) {
    # Filter the data for the two groups in the current pair.
    test_data <- Cell_Summary_by_Day %>%
      filter(
        Treatment_Condition_short %in% pair
      )
    
    # Perform the One-Way ANOVA
    # The formula Fold_Change ~ Treatment_Condition_short is suitable for ANOVA
    anova_result <- aov(Fold_Change ~ Treatment_Condition_short, data = test_data)
    
    # Get the summary of the ANOVA, which includes the p-value
    anova_summary <- summary(anova_result)
    
    # Extract the p-value from the ANOVA summary
    # The p-value for the 'Treatment_Condition_short' factor is typically in the 5th column
    # and the 1st row of the summary table's data part.
    p_value_anova <- anova_summary[[1]]$`Pr(>F)`[1]
    
    # Return a tidy data.table with the results.
    data.table(
      COHORT1 = pair[1],
      COHORT2 = pair[2],
      p_value = signif(p_value_anova, 3)
    )
  }
  
  # Apply the function to all pairs and combine the results.
  p_value_Table <- rbindlist(lapply(treatment_combinations, run_one_way_anova_test))
  
  # Print the results table.
  print(p_value_Table)
}

# EXPORT P_VALUE DATA

# Write the summarized data tables to CSV files for record-keeping or other analyses.
write.csv(p_value_Table, file.path(Plot_Directory_Save_Path, "00_PairwiseOnewayAnova_p_values.csv"))

rm(p_value_Table)

# -----------------------------------------------------------------------------
# GGPLOT
# -----------------------------------------------------------------------------

Plot <- ggplot(
  data = Cell_Summary_by_Cohort, 
  aes(
    x = Treatment_Condition_short, 
    y = Mean_Fold_Change
  )
) +
  geom_bar(
    data = Cell_Summary_by_Cohort,
    stat = "identity", 
    aes(
      x = Treatment_Condition_short,
      y = Mean_Fold_Change,
      fill = COLOR
    ), 
    color = "black", 
    width = 0.4, 
    size = 0.2
  ) +
  scale_fill_identity() + # Use the colors defined in the COLOR column.
  geom_errorbar(
    aes(
      ymin = YMIN, 
      ymax = YMAX
    ), 
    width = 0.2, 
    linewidth = 0.4
  ) +
  geom_jitter(
    data = Cell_Summary_by_Day, 
    aes(
      y = Fold_Change
    ), 
    width = 0.05,
    size = 1, 
    stroke = 0.2, 
    shape = 21, 
    fill = "grey", 
    alpha = 0.8
  ) +
  labs(
    y = "IL-2 Fold Change\n(Stimulated/Unstimulated)"
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 3), 
    limits = c(Y_AXIS_LOWER_LIMIT, Y_AXIS_UPPER_LIMIT), 
    expand = c(0, 0)
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 7),
    legend.position = "none"
  )

Plot_Save_Path_1 <- "01_ELISA Plot.pdf"
Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Plot_Save_Path_1)
ggsave(
  Plot_Save_Path,
  plot = Plot,
  height = 35,
  width = 40,
  units = "mm"
)

# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()
