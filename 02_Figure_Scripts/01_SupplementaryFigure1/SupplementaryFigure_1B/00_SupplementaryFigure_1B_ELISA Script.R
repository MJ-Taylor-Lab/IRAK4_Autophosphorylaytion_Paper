# Load required packages using pacman for automatic installation and loading.
library(pacman)
pacman::p_load(data.table, tidyr, dplyr, ggplot2)

# --- User-defined Parameters ---

# Path to the main input table which lists all the data files for each plate.
Input_Table <- fread("/Users/u_niranjan/Desktop/Cell_Figures/00_New_Figures/Supplementary Figure 1/ELISA - MyD88-GFP_IRAK4-WT_or_K213A214A_DD_KO-mScarlet/TablestoAnalyse.csv")

# Directory where the final plots will be saved.
Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots"

# Dilution factor to apply to the final concentration calculations.
Dilution_Factor <- 5

# --- Directory Creation ---

# Construct and create the specific output directories for the plots if they don't already exist.
Plot_Directory_Path <- file.path(Plot_Directory_Save_Path, "01_SupplementaryFigure1", "SupplementaryFigure_1B")
if (!dir.exists(Plot_Directory_Path)) {
  dir.create(Plot_Directory_Path, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# ELISA DATA PROCESSING FUNCTION
# -----------------------------------------------------------------------------

#' @title Process a single ELISA plate.
#' @description This function reads raw data for one plate, generates a standard
#'              curve, calculates concentrations, and returns a processed data table.
#' @param Plate_Count An integer representing the row number from `Input_Table` to process.
#' @return A data.table with calculated concentrations for the samples.

ELISA_Fx <- function(Plate_Count) {
  
  # --- 1. Read Plate Data ---
  # Get file paths for the current plate from the main input table.
  Values_Measured_Path <- Input_Table$Values_Measured[Plate_Count]
  Treatment_Condition_Path <- Input_Table$Treatment_Conditions[Plate_Count]
  Stimulation_Condition_Path <- Input_Table$Stimulation_Condition[Plate_Count]
  Sample_Day_Path <- Input_Table$Sample_Day[Plate_Count]
  
  # Read the data from the CSV files and combine them into a single data.table.
  Plate <- data.table(
    Values_Measured = as.vector(as.matrix(fread(Values_Measured_Path, header = FALSE))),
    Treatment_Condition = as.vector(as.matrix(fread(Treatment_Condition_Path, header = FALSE))),
    Stimulation_Condition = as.vector(as.matrix(fread(Stimulation_Condition_Path, header = FALSE))),
    Sample_Day = as.vector(as.matrix(fread(Sample_Day_Path, header = FALSE)))
  )
  
  # Remove wells that are marked as "Blank".
  Plate <- Plate[Treatment_Condition != "Blank"]
  
  # --- 2. Standard Curve Generation and Fitting ---
  # Isolate the calibration standards from the plate data.
  Plate_Standards <- Plate[Stimulation_Condition == "Calibration"] %>%
    group_by(
      Treatment_Condition
    ) %>%
    summarise(
      Values_Measured_mean = mean(as.numeric(Values_Measured), na.rm = TRUE)
    ) %>%
    mutate(
      Treatment_Condition = as.numeric(Treatment_Condition)
    ) %>%
    arrange(Treatment_Condition)
  
  # Fit a linear model: Concentration ~ Measured_Value
  Fit_1 <- lm(Treatment_Condition ~ Values_Measured_mean, data = Plate_Standards)
  R_1 <- summary(Fit_1)$r.squared
  
  # Check if removing the highest standard point results in a better fit.
  Plate_Standards_Temp <- head(Plate_Standards, -1)
  Fit_2 <- lm(Treatment_Condition ~ Values_Measured_mean, data = Plate_Standards_Temp)
  R_2 <- summary(Fit_2)$r.squared
  
  # Select the model with the higher R-squared value.
  if (R_2 - 0.01 > R_1) {
    Fit <- Fit_2
    Rsquare <- R_2
    Plate_Standards <- Plate_Standards_Temp
  } else {
    Fit <- Fit_1
    Rsquare <- R_1
  }
  
  # --- 3. Plot and Save the Standard Curve ---
  # Add the fitted values to the standards data for plotting the regression line.
  Plate_Standards <- Plate_Standards %>%
    mutate(Fit_Test = predict(Fit, newdata = .))
  
  # Create a plot of the standard curve.
  Calibration_Plot <- ggplot(
    data = Plate_Standards, 
    aes(
      x = Values_Measured_mean, 
      y = Treatment_Condition
    )
  ) +
    geom_point(
      size = 5
    ) +
    geom_line(
      aes(y = Fit_Test), 
      linetype = "dashed"
    ) +
    annotate(
      'text', 
      x = 0.15, 
      y = 700, 
      label = paste0("R^2 = ", signif(Rsquare, 4)), 
      size = 10
    ) +
    annotate(
      'text',
      x = max(Plate_Standards$Values_Measured_mean) * 0.75,
      y = 150,
      label = paste0("IL-Amount = \n", signif(coefficients(Fit)[2], 4), "*Intensity \n+ (", signif(coefficients(Fit)[1], 4), ")")
    ) +
    labs(
      x = "Measured Values", 
      y = "IL-Concentration (pg/mL)", 
      title = paste("Plate", Plate_Count, "Standard Curve")
    ) +
    theme_classic(
      base_size = 20
    )
  
  # Save the plot to a PDF file.
  ggsave(
    filename = file.path(Plot_Directory_Path, paste0("Plate_Input_Number_", Plate_Count, "_Calibration.pdf")),
    plot = Calibration_Plot,
    height = 9,
    width = 20
  )
  
  # --- 4. Calculate Concentrations for Samples ---
  # Filter out the calibration standards to process only the experimental samples.
  Plate_Samples <- Plate[Stimulation_Condition != "Calibration"] %>%
    # Use the coefficients from the fitted linear model (y = mx + c) to calculate
    # the concentration for each sample from its measured value.
    # Concentration = slope * Measured_Value + intercept
    mutate(
      Values_Measured = as.numeric(Values_Measured),
      IL_conc_Estimation = (coefficients(Fit)[2] * Values_Measured) + coefficients(Fit)[1]
    )
  
  return(Plate_Samples)
}


# -----------------------------------------------------------------------------
# DATA AGGREGATION AND PREPARATION
# -----------------------------------------------------------------------------

# Apply the ELISA_Fx function to each plate defined in Input_Table and combine results.
# `lapply` processes each row, and `rbindlist` efficiently stacks the resulting data.tables.
SplitPlateTable <- rbindlist(lapply(1:nrow(Input_Table), ELISA_Fx))

# Define the desired order for treatments. This replaces the complex Cohort_Order_Table logic.
cohort_order_base <- c("KO", "WT", "K213A/214A", "DD")

# Process the combined data.
Processed_Table <- SplitPlateTable %>%
  # Filter out any irrelevant conditions.
  filter(Treatment_Condition != "MyD88-GFP") %>%
  mutate(
    # Set any negative estimated concentrations to a small positive number (e.g., 1) to avoid math errors (e.g., in log transforms or ratios).
    IL_conc_Estimation = ifelse(IL_conc_Estimation < 0, 1, IL_conc_Estimation),
    # Apply the dilution factor.
    IL_conc_Estimation = IL_conc_Estimation * Dilution_Factor,
    
    # Create short, clean labels for plotting.
    Treatment_Condition_short = case_when(
      Treatment_Condition == "IRAK4 KO" ~ "KO",
      Treatment_Condition == "IRAK4 WT" ~ "WT",
      Treatment_Condition == "IRAK4 KD" ~ "K213A/214A",
      Treatment_Condition == "IRAK4 DD" ~ "DD"
    ),
    Stimulation_Condition_short = case_when(
      Stimulation_Condition == "Unstimulated" ~ "_US",
      TRUE ~ "_S"
    ),
    # Combine treatment and stimulation labels for unique group identification.
    Treatment_Condition_short_Extended = paste0(Treatment_Condition_short, Stimulation_Condition_short)
  )

# Set the factor levels for plotting to ensure the bars appear in the desired order.
plot_order_levels <- unlist(lapply(cohort_order_base, function(x) paste0(x, c("_US", "_S"))))

Processed_Table <- Processed_Table %>%
  mutate(
    Treatment_Condition_short_Extended = factor(Treatment_Condition_short_Extended, levels = plot_order_levels),
    Treatment_Condition_short = factor(Treatment_Condition_short, levels = cohort_order_base)
  ) %>%
  # Arrange the table based on the new factor levels for cleaner viewing.
  arrange(
    Treatment_Condition_short_Extended
  )


# -----------------------------------------------------------------------------
# PLOT 1: ABSOLUTE IL-2 CONCENTRATION
# -----------------------------------------------------------------------------

# Summarize data by day/replicate for individual data points on the plot.
PlateSummarybyDay <- Processed_Table %>%
  group_by(Treatment_Condition_short, Treatment_Condition_short_Extended, Sample_Day) %>%
  summarise(IL_conc_Estimation_mean = mean(IL_conc_Estimation), .groups = "drop")

# Summarize data across all replicates for the main bar plot values and error bars.
PlateSummarybyTreatment <- Processed_Table %>%
  group_by(Treatment_Condition_short, Treatment_Condition_short_Extended) %>%
  summarise(
    IL_conc_Estimation_mean = mean(IL_conc_Estimation),
    IL_conc_Estimation_sd = sd(IL_conc_Estimation),
    .groups = "drop"
  ) %>%
  mutate(
    # Set Y_MIN to 0 if the calculation results in a negative value
    Y_MIN = pmax(0, IL_conc_Estimation_mean - IL_conc_Estimation_sd),
    Y_MAX = IL_conc_Estimation_mean + IL_conc_Estimation_sd,
    # Define colors for the bars.
    COLOR = case_when(
      Treatment_Condition_short == "WT"  ~ "#F99108",
      TRUE ~ "#33AFCC"
    )
  )

# Create the bar plot with error bars and individual points.
ggplot(
  data = PlateSummarybyTreatment, 
  aes(
    x = Treatment_Condition_short_Extended, 
    y = IL_conc_Estimation_mean
  )
) +
  geom_bar(
    stat = "identity", 
    aes(
      fill = COLOR
    ), 
    color = "black"
  ) +
  scale_fill_identity() + # Use the colors defined in the COLOR column.
  geom_errorbar(
    aes(
      ymin = Y_MIN, 
      ymax = Y_MAX
    ), 
    width = 0.5, 
    linewidth = 0.5
  ) +
  geom_point(
    data = PlateSummarybyDay, 
    aes(
      y = IL_conc_Estimation_mean
    ), 
    size = 2, 
    shape = 16, 
    fill = "grey"
  ) +
  labs(
    y = "IL-2 conc. (pg/mL)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 20, angle = 45, vjust = 0.55),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 25),
    legend.position = "none"
  )

# Save the plot.
ggsave(
  file.path(Plot_Directory_Path, "01_Real_IL2_pg_per_ml_conc.pdf"),
  plot = last_plot(),
  height = 9,
  width = 20
)

# -----------------------------------------------------------------------------
# PLOT 2: FOLD CHANGE (STIMULATED / UNSTIMULATED)
# -----------------------------------------------------------------------------

# Calculate the mean IL-2 concentration for each condition and replicate.
FoldChangeSummary_Replicates <- Processed_Table %>%
  group_by(Treatment_Condition_short, Sample_Day, Stimulation_Condition) %>%
  summarise(IL_conc_mean = mean(IL_conc_Estimation), .groups = 'drop') %>%
  # Reshape data so 'Stimulated' and 'Unstimulated' are columns.
  pivot_wider(
    names_from = Stimulation_Condition,
    values_from = IL_conc_mean
  ) %>%
  # Calculate the fold change.
  mutate(Fold_Change = Stimulated / Unstimulated)

# Summarize the fold change data across all replicates for the bar plot.
FoldChangeSummary_Overall <- FoldChangeSummary_Replicates %>%
  group_by(Treatment_Condition_short) %>%
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
      Treatment_Condition_short == "WT"  ~ "#F99108",
      TRUE ~ "#33AFCC"
    )
  )

# Create the fold change plot.
ggplot(
  data = FoldChangeSummary_Overall, 
  aes(
    x = Treatment_Condition_short, 
    y = Mean_Fold_Change
  )
) +
  geom_bar(
    stat = "identity", 
    aes(
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
    data = FoldChangeSummary_Replicates, 
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
    limits = c(-0.5, 23), 
    expand = c(0, 0)
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 7),
    legend.position = "none"
  )

# Save the plot.
ggsave(
  file.path(Plot_Directory_Path, "02_Fold_Change_IL2_Unstimulated_to_Stimulated_IRAK4KI.pdf"),
  plot = last_plot(),
  height = 40,
  width = 65,
  units = "mm"
)

# -----------------------------------------------------------------------------
# STATISTICAL ANALYSIS: PAIRWISE ONE WAY ANOVA TEST
# -----------------------------------------------------------------------------

# Perform pairwise One way anova test to compare fold changes between treatment groups.
p_value_test_confirmation <- "Y"
if (p_value_test_confirmation == "Y") {
  
  # Get the unique treatment conditions to compare.
  treatments_to_compare <- as.character(unique(FoldChangeSummary_Replicates$Treatment_Condition_short))
  
  # Generate all unique pairs of treatments.
  treatment_combinations <- combn(sort(treatments_to_compare), 2, simplify = FALSE)
  
  # Define a function to perform the test for a single pair.
  run_one_way_anova_test <- function(pair) {
    # Filter the data for the two groups in the current pair.
    test_data <- FoldChangeSummary_Replicates %>%
      filter(Treatment_Condition_short %in% pair)
    
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

# -----------------------------------------------------------------------------
# EXPORT SUMMARY DATA
# -----------------------------------------------------------------------------

# Write the summarized data tables to CSV files for record-keeping or other analyses.
write.csv(FoldChangeSummary_Replicates, file.path(Plot_Directory_Path, "00_ELISASummary_by_Replicate.csv"))
write.csv(FoldChangeSummary_Overall, file.path(Plot_Directory_Path, "00_ELISASummary_by_Cohort.csv"))
write.csv(p_value_Table, file.path(Plot_Directory_Path, "00_Pairwise_p_values.csv"))


# Clean up ----------------------------------------------------------------
rm(list = ls())
gc()
