library(pacman)
pacman::p_load(data.table, tidyr, dplyr, ggplot2,interleave)

Input_Table <- fread("/Users/u_niranjan/Desktop/Possible paper figure/Git Figure_20231213/12_Supplementary Figure 2/20230607/20230607_1/TablestoAnalyse.csv")
Output_Directory <- "/Users/u_niranjan/Desktop/Possible paper figure/Git Figure_20231213/12_Supplementary Figure 2/20230607/20230607_1"

Dilution_Factor <- 5

Cohort_Order_Table_Path <- "/Users/u_niranjan/Desktop/Possible paper figure/Git Figure_20231213/12_Supplementary Figure 2/20230607/20230607_1/CohortOrder.csv" #Put nonsense string if you dont want order 
if(file.exists(Cohort_Order_Table_Path)){
  Cohort_Order_Table <- read.csv(Cohort_Order_Table_Path, header = FALSE)
  Cohort_Order_Table <- Cohort_Order_Table %>% as.vector()
  
  Cohort_Order_Table_Stim <- Cohort_Order_Table %>% as.data.table()
  colnames(Cohort_Order_Table_Stim) <- c("Cohort")
  Cohort_Order_Table_Stim <- Cohort_Order_Table_Stim %>% 
    mutate(
      Cohort = paste0(Cohort, "_Stimulated")
    )
  
  Cohort_Order_Table_Unstim <- Cohort_Order_Table %>% as.data.table()
  colnames(Cohort_Order_Table_Unstim) <- c("Cohort")
  Cohort_Order_Table_Unstim <- Cohort_Order_Table_Unstim %>% 
    mutate(
      Cohort = paste0(Cohort, "_Unstimulated")
    ) 
  Cohort_Order_Table_2 <- rbind(Cohort_Order_Table_Unstim, Cohort_Order_Table_Stim)
  
  MAX_ROWS <- max(nrow(Cohort_Order_Table_Unstim), nrow(Cohort_Order_Table_Stim))
  Cohort_Order_Table_3 <- data.frame(ID = numeric(0), Value = character(0))
  
  # Loop through the rows of both tables and alternately add them to the combined table
  for (index in 1:MAX_ROWS) {
    if (index <= nrow(Cohort_Order_Table_Unstim)) {
      Cohort_Order_Table_3 <- rbind(Cohort_Order_Table_3, Cohort_Order_Table_Unstim[index, ])
    }
    if (index <= nrow(table2)) {
      Cohort_Order_Table_3 <- rbind(Cohort_Order_Table_3, Cohort_Order_Table_Stim[index, ])
    }
  }
  
  Cohort_Order_Table <- Cohort_Order_Table %>% as.data.table()
  
  rm(
    index,
    Cohort_Order_Table_2,
    Cohort_Order_Table_Stim,
    Cohort_Order_Table_Unstim
  )
  print("Cohort Order Table Ready")
} else {
  print("No Cohort Table available")
}


#Creating Output Folder
Output_Directory <- file.path(Output_Directory, "Output")
if(!file.exists(Output_Directory)){
  dir.create(Output_Directory)
}

ELISA_Fx <- function(Plate_Count){

# Reading Data ------------------------------------------------------------
  #Getting Path of treatment conditions and values
  Values_Measured_Path <- Input_Table$Values_Measured[Plate_Count]
  Treatment_Condition_Path <- Input_Table$Treatment_Conditions[Plate_Count]
  Stimulation_Condition_Path <- Input_Table$Stimulation_Condition[Plate_Count]
  Sample_Day_Path <- Input_Table$Sample_Day[Plate_Count]
  
  #Reading Plate Treatement 
  Values_Measured <- fread(Values_Measured_Path, header = F)
  Treatment_Condition <- fread(Treatment_Condition_Path, header = F)
  Stimulation_Condition <- fread(Stimulation_Condition_Path, header = F)
  Sample_Day <- fread(Sample_Day_Path, header = F)
  
  rm(
    Values_Measured_Path,
    Treatment_Condition_Path,
    Stimulation_Condition_Path,
    Sample_Day_Path
  )
  
  #Converting tables into vector for to make a single table
  Values_Measured <- as.vector(as.matrix(Values_Measured))
  Treatment_Condition <- as.vector(as.matrix(Treatment_Condition))
  Stimulation_Condition <- as.vector(as.matrix(Stimulation_Condition))
  Sample_Day <- as.vector(as.matrix(Sample_Day))
  
  #Creating Table containing all plate Information
  Plate <- NULL
  Plate$Values_Measured <- Values_Measured
  Plate$Treatment_Condition <- Treatment_Condition
  Plate$Stimulation_Condition <- Stimulation_Condition
  Plate$Sample_Day <- Sample_Day
  
  rm(
    Values_Measured,
    Treatment_Condition,
    Stimulation_Condition,
    Sample_Day
  )
  
  Plate <- Plate %>% as.data.table()
  
  #Removing Empty Wells
  Plate <- Plate %>% 
    filter(
      Treatment_Condition != "Blank"
    ) %>% as.data.table()
  

# Standard Curve ----------------------------------------------------------
  #Creating a Standard Curve
  Plate_Standards <- Plate %>% 
    filter(
      Stimulation_Condition == "Calibration"
    ) %>% 
    group_by(
      Treatment_Condition
    ) %>% 
    summarise(
      Values_Measured_mean = mean(Values_Measured)
      # Values_Measured_median = median(Values_Measured)
    ) %>%  
    mutate(
      Treatment_Condition = as.numeric(Treatment_Condition)
    ) %>% 
    arrange(
      Treatment_Condition
    )
    
  Fit_1 <- lm(Treatment_Condition ~ Values_Measured_mean, data = Plate_Standards)
  R_1 <- summary(Fit_1)$r.squared
  
  Plate_Standards_Temp <- head(Plate_Standards, -1)
  Fit_2 <- lm(Treatment_Condition ~ Values_Measured_mean, data = Plate_Standards_Temp)
  R_2 <- summary(Fit_2)$r.squared
  
  if (R_2-0.01 > R_1) {
    Plate_Standards <- Plate_Standards_Temp
    Fit <- Fit_2
    Rsquare <- R_2
  } else {
    Plate_Standards <- Plate_Standards
    Fit <- Fit_1
    Rsquare <- R_1
  }
  
  Rsquare <- signif(Rsquare, digits = 4)
  
  rm(
    Plate_Standards_Temp,
    Fit_2,
    Fit_1,
    R_2,
    R_1
  )
  
  
  print(paste0("IL-Amount = slope*Intensity + X-Intercept"))
  print(paste0("IL-Amount = ", Fit$coefficients[2],"*Intensity + ", Fit$coefficients[1]))
  
  Plate_Standards <- Plate_Standards %>% 
    mutate(
      Fit_Test = (Fit$coefficients[2]*Values_Measured_mean) + Fit$coefficients[1]
    )
  
  ggplot(
    data = Plate_Standards,
  ) +
    geom_point(
      aes(
        x = Values_Measured_mean,
        y = Treatment_Condition
      ),
      size = 5
    ) +
    geom_line(
      aes(
        x = Values_Measured_mean,
        y = Fit_Test
      ),
      linetype = "dashed"
    ) +
    annotate(
      'text',
      x = 0.15,
      y = 700,
      label = paste0("R^2 = ",Rsquare),
      size = 10
    ) +
    annotate(
      'text',
      x = max(Plate_Standards$Values_Measured_mean) - (0.25*max(Plate_Standards$Values_Measured_mean)),
      y = 150,
      label = paste0("IL-Amount = \n", signif(Fit$coefficients[2], digits = 4),"*Intensity \n+ (", signif(Fit$coefficients[1], digits = 4), ")")
    ) +
    labs(
      x = "Measured Values",
      y = "IL-Concentration (pg/mL)"
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 20)
    )
  
  Save_Name <- paste0("Plate_Input_Number_" , Plate_Count, "_Calibration.pdf")
  Save_Name <- file.path(Output_Directory, Save_Name)

  ggsave(
    Save_Name,
    plot = last_plot(),
    height = 3*3,
    width = 5*4
  )
  
  rm(
    Save_Name,
    Plate_Standards,
    Rsquare
  )

  
# Fitting Data To Standarad Curve -----------------------------------------
  Plate <- Plate %>% 
    filter(
      Stimulation_Condition != "Calibration"
    ) %>% 
    mutate(
      Values_Measured = as.numeric(Values_Measured),
      IL_conc_Estimation = (Fit$coefficients[2]*Values_Measured) + Fit$coefficients[1],
    )
  
  rm(
    Fit
  )
  
  return(Plate)
}

SplitPlateTable <- lapply(1:NROW(Input_Table), ELISA_Fx)
SplitPlateTable <- rbindlist(SplitPlateTable)
  
Cohort_Order_Table_3 <- Cohort_Order_Table_3 %>% 
  mutate(
    ORDER_NUMBER = row_number()
  )

SplitPlateTable <- SplitPlateTable %>% 
  mutate(
    Stimulation_Condition_Extended = paste0(Treatment_Condition, "_", Stimulation_Condition),
    ORDER_NUMBER = Cohort_Order_Table_3$ORDER_NUMBER[match(Stimulation_Condition_Extended, Cohort_Order_Table_3$Cohort)],
    Values_Measured = case_when(
      Values_Measured < 0 ~ 0,
      TRUE ~ Values_Measured
    ),
    IL_conc_Estimation = case_when(
      IL_conc_Estimation < 0 ~ 1,
      TRUE ~ IL_conc_Estimation
    ),
    Treatment_Condition_short = case_when(
      Treatment_Condition == "WildType_NoTreatment" ~ "WT",
      Treatment_Condition == "WildType+DMSO"  ~ "WT+DMSO",
      Treatment_Condition == "WildType+500nM" ~ "WT+500nM",
      Treatment_Condition == "WildType+20uM" ~ "WT+20uM",
      Treatment_Condition == "IRAK4KI_NoTreatment" ~ "IRAK4KI",
      Treatment_Condition == "IRAK4KI+DMSO"  ~ "IRAK4KI+DMSO",
      Treatment_Condition == "IRAK4KI+500nM" ~ "IRAK4KI+500nM",
      Treatment_Condition == "IRAK4KI+20uM" ~ "IRAK4KI+20uM",
      Treatment_Condition == "IRAK1KI_NoTreatment" ~ "IRAK1KI",
      Treatment_Condition == "IRAK1KI+DMSO"  ~ "IRAK1KI+DMSO",
      Treatment_Condition == "IRAK1KI+500nM" ~ "IRAK1KI+500nM",
      Treatment_Condition == "IRAK1KI+20uM" ~ "IRAK1KI+20uM"
    ),
    Stimulation_Condition_short = case_when(
      Stimulation_Condition == "Unstimulated" ~ "_US",
      TRUE ~ "_S"
    ),
    Treatment_Condition_short_Extended = paste0(Treatment_Condition_short, Stimulation_Condition_short),
    IL_conc_Estimation = IL_conc_Estimation*Dilution_Factor
  ) %>% 
  arrange(
    ORDER_NUMBER,
    Treatment_Condition,
    Stimulation_Condition
  )


# Graph with real pg/ml ---------------------------------------------------
# Getting mean of each column on a day to day basis
PlateSummarybyDay <- SplitPlateTable %>% 
  group_by(
    Treatment_Condition_short,
    Treatment_Condition_short_Extended,
    ORDER_NUMBER,
    Sample_Day,
    Stimulation_Condition
  ) %>% 
  summarise(
    IL_conc_Estimation_mean = mean(IL_conc_Estimation),
    IL_conc_Estimation_median = median(IL_conc_Estimation)
  ) %>% 
  as.data.table()
  
#Getting overall value of a single Cohort
PlateSummarybyTreatment_Condition <- SplitPlateTable %>% 
  group_by(
    Treatment_Condition_short,
    Treatment_Condition_short_Extended,
    ORDER_NUMBER,
    Stimulation_Condition
  ) %>% 
  summarise(
    IL_conc_Estimation_mean = mean(IL_conc_Estimation),
    IL_conc_Estimation_median = median(IL_conc_Estimation),
    IL_conc_Estimation_sd = sd(IL_conc_Estimation),
    Y_MAX = IL_conc_Estimation_median + IL_conc_Estimation_sd,
    Y_MIN = case_when(
      IL_conc_Estimation_median - IL_conc_Estimation_sd < 0 ~ 0,
      TRUE ~ IL_conc_Estimation_median - IL_conc_Estimation_sd,
    )
  ) %>% 
  mutate(
    COLOR = case_when(
      Stimulation_Condition == "Stimulated" ~ "darkred",
      TRUE ~ "darkblue"
    )
  ) %>%
  as.data.table()

PlateSummarybyTreatment_Condition$Treatment_Condition_short_Extended <- factor(
  PlateSummarybyTreatment_Condition$Treatment_Condition_short_Extended, 
  levels = PlateSummarybyTreatment_Condition$Treatment_Condition_short_Extended[order(PlateSummarybyTreatment_Condition$ORDER_NUMBER)]
)


ggplot(
  data = PlateSummarybyTreatment_Condition,
  aes(
    x = Treatment_Condition_short_Extended,
    y = IL_conc_Estimation_mean
  )
) +
  geom_bar(
    stat = "identity",
    fill = PlateSummarybyTreatment_Condition$COLOR
  ) +
  geom_errorbar(
    data = PlateSummarybyTreatment_Condition,
    aes(
      x = Treatment_Condition_short_Extended,
      ymin = Y_MAX,
      ymax = Y_MIN
    ),
    width = 0.5,
    linewidth = 0.5
  ) +
  geom_point(
    data = PlateSummarybyDay,
    aes(
      x = Treatment_Condition_short_Extended,
      y = IL_conc_Estimation_mean
    ),
    size = 2,
    color = "black",
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
    axis.title.y = element_text(size = 25)
  )

Plot_Save <- file.path(Output_Directory, "01_Real_IL2_pg_per_ml_conc.pdf")
ggsave(
  Plot_Save,
  plot = last_plot(),
  height = 3*3,
  width = 5*4
)

# Graph with ratio to unstimulated---------------------------------------------------
# Getting mean of each column on a day to day basis
PlateSummarybyDay <- SplitPlateTable %>% 
  group_by(
    Treatment_Condition_short,
    Stimulation_Condition,
    ORDER_NUMBER,
    Sample_Day
  ) %>% 
  summarise(
    IL_conc_Estimation_mean = mean(IL_conc_Estimation),
    # IL_conc_Estimation_median = median(IL_conc_Estimation)
  ) %>%
  as.data.table()


PlateSummarybyDay_Unstimulated <- PlateSummarybyDay %>% 
  filter(
    Stimulation_Condition == "Unstimulated"
  ) %>% 
  select(
    -Stimulation_Condition,
    -ORDER_NUMBER
  )

names(PlateSummarybyDay_Unstimulated)[names(PlateSummarybyDay_Unstimulated) == "IL_conc_Estimation_mean"] <- "IL_conc_Estimation_mean_Unstimulated"

PlateSummarybyDay_Stimulated <- PlateSummarybyDay %>% 
  filter(
    Stimulation_Condition == "Stimulated"
  ) %>% 
  select(
    -Stimulation_Condition
  )

PlateSummarybyDay <- PlateSummarybyDay_Stimulated %>% 
  inner_join(
    PlateSummarybyDay_Unstimulated,
    by = c("Treatment_Condition_short" = "Treatment_Condition_short", "Sample_Day" = "Sample_Day")
  ) %>% 
  mutate(
    Unstimulated_by_Stimulated = IL_conc_Estimation_mean/ IL_conc_Estimation_mean_Unstimulated
  ) %>% 
  arrange(
    ORDER_NUMBER
  )
  

rm(
  PlateSummarybyDay_Stimulated,
  PlateSummarybyDay_Unstimulated
)

#Getting overall value of a single Cohort
PlateSummarybyTreatment_Condition <- PlateSummarybyDay %>% 
  group_by(
    Treatment_Condition_short,
    ORDER_NUMBER
  ) %>% 
  summarise(
    Ratio_sd = sd(Unstimulated_by_Stimulated),
    Unstimulated_by_Stimulated = mean(Unstimulated_by_Stimulated)
  ) %>% 
  as.data.table()

PlateSummarybyTreatment_Condition$Treatment_Condition_short <- factor(
  PlateSummarybyTreatment_Condition$Treatment_Condition_short, 
  levels = PlateSummarybyTreatment_Condition$Treatment_Condition_short[order(PlateSummarybyTreatment_Condition$ORDER_NUMBER)]
)



# Seperating Table by CellLine WT--------------------------------------------
PlateSummarybyDay_CellLine <- PlateSummarybyDay %>% 
  filter(
    # Treatment_Condition_short == "WT" | 
   Treatment_Condition_short == "WT+DMSO" | Treatment_Condition_short == "WT+500nM" | Treatment_Condition_short == "WT+20uM"
  )

PlateSummarybyTreatment_Condition_CellLine <- PlateSummarybyTreatment_Condition %>% 
  filter(
    # Treatment_Condition_short == "WT" | 
    Treatment_Condition_short == "WT+DMSO" | Treatment_Condition_short == "WT+500nM" | Treatment_Condition_short == "WT+20uM"
  ) %>% 
  mutate(
    COLOR = case_when(
      # Treatment_Condition_short == "WT" ~ "#FCCC33",  #yellowish orange,
      Treatment_Condition_short == "WT+DMSO" ~ "orange",
      Treatment_Condition_short == "WT+500nM" ~ "#33FFCC",# greenish lightblue
      Treatment_Condition_short == "WT+20uM" ~ "lightblue" 
    )
  )


ggplot(
  data = PlateSummarybyTreatment_Condition_CellLine,
  aes(
    x = Treatment_Condition_short,
    y = Unstimulated_by_Stimulated
  )
) +
  geom_bar(
    stat = "identity",
    fill = PlateSummarybyTreatment_Condition_CellLine$COLOR,
    color = "black",
    width = 0.6,
  ) +
  geom_errorbar(
    data = PlateSummarybyTreatment_Condition_CellLine,
    aes(
      x = Treatment_Condition_short,
      ymin = Unstimulated_by_Stimulated - Ratio_sd,
      ymax = Unstimulated_by_Stimulated + Ratio_sd
    ),
    width = 0.4,
    linewidth = 0.4
  ) +
  geom_point(
    data = PlateSummarybyDay_CellLine,
    aes(
      x = Treatment_Condition_short,
      y = Unstimulated_by_Stimulated
    ),
    size = 1,
    shape = 21,
    fill = "grey",
    alpha = 0.8,
    position = position_jitter(0.05)
  ) +
  labs(
    y = "IL-2 Fold Change\n(Stimulated/Unstimulated)"
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 6)
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 6)
  )

Plot_Save <- file.path(Output_Directory, "02_Fold_Change_IL2_Unstimulated_to_Stimulated_WT.pdf")
ggsave(
  Plot_Save,
  plot = last_plot(),
  height = 35,
  width = 60,
  units = "mm"
)

# Seperating Table by CellLine IRAK4Ki--------------------------------------------
PlateSummarybyDay_CellLine <- PlateSummarybyDay %>% 
  filter(
    # Treatment_Condition_short == "IRAK4KI" | 
    Treatment_Condition_short == "IRAK4KI+DMSO" | Treatment_Condition_short == "IRAK4KI+500nM" | Treatment_Condition_short == "IRAK4KI+20uM"
  )

PlateSummarybyTreatment_Condition_CellLine <- PlateSummarybyTreatment_Condition %>% 
  filter(
    # Treatment_Condition_short == "IRAK4KI" | 
    Treatment_Condition_short == "IRAK4KI+DMSO" | Treatment_Condition_short == "IRAK4KI+500nM" | Treatment_Condition_short == "IRAK4KI+20uM"
  ) %>% 
  mutate(
    COLOR = case_when(
      # Treatment_Condition_short == "IRAK4KI" ~ "#FCCC33",  #yellowish orange,
      Treatment_Condition_short == "IRAK4KI+DMSO" ~ "orange",
      Treatment_Condition_short == "IRAK4KI+500nM" ~ "#33FFCC",# greenish lightblue
      Treatment_Condition_short == "IRAK4KI+20uM" ~ "lightblue" 
    )
  )


ggplot(
  data = PlateSummarybyTreatment_Condition_CellLine,
  aes(
    x = Treatment_Condition_short,
    y = Unstimulated_by_Stimulated
  )
) +
  geom_bar(
    stat = "identity",
    fill = PlateSummarybyTreatment_Condition_CellLine$COLOR,
    color = "black",
    width = 0.6,
  ) +
  geom_errorbar(
    data = PlateSummarybyTreatment_Condition_CellLine,
    aes(
      x = Treatment_Condition_short,
      ymin = Unstimulated_by_Stimulated - Ratio_sd,
      ymax = Unstimulated_by_Stimulated + Ratio_sd
    ),
    width = 0.4,
    linewidth = 0.4
  ) +
  geom_point(
    data = PlateSummarybyDay_CellLine,
    aes(
      x = Treatment_Condition_short,
      y = Unstimulated_by_Stimulated
    ),
    size = 1,
    shape = 21,
    fill = "grey",
    alpha = 0.8,
    position = position_jitter(0.05)
  ) +
  labs(
    y = "IL-2 Fold Change\n(Stimulated/Unstimulated)"
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 6)
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 6)
  )

Plot_Save <- file.path(Output_Directory, "03_Fold_Change_IL2_Unstimulated_to_Stimulated_IRAK4KI.pdf")
ggsave(
  Plot_Save,
  plot = last_plot(),
  height = 35,
  width = 60,
  units = "mm"
)

# Seperating Table by CellLine IRAK1Ki--------------------------------------------
PlateSummarybyDay_CellLine <- PlateSummarybyDay %>% 
  filter(
    # Treatment_Condition_short == "IRAK1KI" | 
    Treatment_Condition_short == "IRAK1KI+DMSO" | Treatment_Condition_short == "IRAK1KI+500nM" | Treatment_Condition_short == "IRAK1KI+20uM"
  )

PlateSummarybyTreatment_Condition_CellLine <- PlateSummarybyTreatment_Condition %>% 
  filter(
    # Treatment_Condition_short == "IRAK1KI" | 
    Treatment_Condition_short == "IRAK1KI+DMSO" | Treatment_Condition_short == "IRAK1KI+500nM" | Treatment_Condition_short == "IRAK1KI+20uM"
  ) %>% 
  mutate(
    COLOR = case_when(
      # Treatment_Condition_short == "IRAK1KI" ~ "#FCCC33",  #yellowish orange,
      Treatment_Condition_short == "IRAK1KI+DMSO" ~ "orange",
      Treatment_Condition_short == "IRAK1KI+500nM" ~ "#33FFCC",# greenish lightblue
      Treatment_Condition_short == "IRAK1KI+20uM" ~ "lightblue" 
    )
  )


ggplot(
  data = PlateSummarybyTreatment_Condition_CellLine,
  aes(
    x = Treatment_Condition_short,
    y = Unstimulated_by_Stimulated
  )
) +
  geom_bar(
    stat = "identity",
    fill = PlateSummarybyTreatment_Condition_CellLine$COLOR,
    color = "black",
    width = 0.6,
  ) +
  geom_errorbar(
    data = PlateSummarybyTreatment_Condition_CellLine,
    aes(
      x = Treatment_Condition_short,
      ymin = Unstimulated_by_Stimulated - Ratio_sd,
      ymax = Unstimulated_by_Stimulated + Ratio_sd
    ),
    width = 0.4,
    linewidth = 0.4
  ) +
  geom_point(
    data = PlateSummarybyDay_CellLine,
    aes(
      x = Treatment_Condition_short,
      y = Unstimulated_by_Stimulated
    ),
    size = 1,
    shape = 21,
    fill = "grey",
    alpha = 0.8,
    position = position_jitter(0.05)
  ) +
  labs(
    y = "IL-2 Fold Change\n(Stimulated/Unstimulated)"
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 6)
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 6)
  )

Plot_Save <- file.path(Output_Directory, "04_Fold_Change_IL2_Unstimulated_to_Stimulated_IRAK1KI.pdf")
ggsave(
  Plot_Save,
  plot = last_plot(),
  height = 35,
  width = 60,
  units = "mm"
)


# Seperating Table by CellLine No Treatment--------------------------------------------
PlateSummarybyDay_CellLine <- PlateSummarybyDay %>% 
  filter(
    Treatment_Condition_short == "WT" | Treatment_Condition_short == "IRAK4KI" | Treatment_Condition_short == "IRAK1KI"
  )

PlateSummarybyTreatment_Condition_CellLine <- PlateSummarybyTreatment_Condition %>% 
  filter(
    # Treatment_Condition_short == "IRAK4KI" | 
    Treatment_Condition_short == "WT" | Treatment_Condition_short == "IRAK4KI" | Treatment_Condition_short == "IRAK1KI"
  ) %>% 
  mutate(
    COLOR = case_when(
      # Treatment_Condition_short == "IRAK4KI" ~ "#FCCC33",  #yellowish orange,
      Treatment_Condition_short == "WT" ~ "orange",
      Treatment_Condition_short == "IRAK4KI" ~ "gold2",# greenish lightblue
      Treatment_Condition_short == "IRAK1KI" ~ "gold4" 
    )
  )


ggplot(
  data = PlateSummarybyTreatment_Condition_CellLine,
  aes(
    x = Treatment_Condition_short,
    y = Unstimulated_by_Stimulated
  )
) +
  geom_bar(
    stat = "identity",
    fill = PlateSummarybyTreatment_Condition_CellLine$COLOR,
    color = "black",
    width = 0.6,
  ) +
  geom_errorbar(
    data = PlateSummarybyTreatment_Condition_CellLine,
    aes(
      x = Treatment_Condition_short,
      ymin = Unstimulated_by_Stimulated - Ratio_sd,
      ymax = Unstimulated_by_Stimulated + Ratio_sd
    ),
    width = 0.4,
    linewidth = 0.4
  ) +
  geom_point(
    data = PlateSummarybyDay_CellLine,
    aes(
      x = Treatment_Condition_short,
      y = Unstimulated_by_Stimulated
    ),
    size = 1,
    shape = 21,
    fill = "grey",
    alpha = 0.8,
    position = position_jitter(0.05)
  ) +
  labs(
    y = "IL-2 Fold Change\n(Stimulated/Unstimulated)"
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 6)
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 6)
  )

Plot_Save <- file.path(Output_Directory, "05_Fold_Change_IL2_Unstimulated_to_Stimulated_NoTreatment.pdf")
ggsave(
  Plot_Save,
  plot = last_plot(),
  height = 35,
  width = 60,
  units = "mm"
)



# Cleanup -----------------------------------------------------------------
rm(list = ls())

