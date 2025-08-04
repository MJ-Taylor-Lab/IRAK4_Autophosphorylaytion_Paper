Number_of_Association_Events_Data_Processing <- function(Cell_Summary_by_Cell){
  # Further processing to group and summarise data by IMAGE and REPLICATE
  Cell_Summary_by_Image <- Cell_Summary_by_Cell %>% 
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
    group_by(
      SHORT_LABEL
    ) %>% 
    mutate(
      REPLICATE_NUMBER = n_distinct(REPLICATE)
    )
  
  # Further processing to arrange and group data by COHORT, then summarise
  Cell_Summary_by_Cohort <- Cell_Summary_by_Image %>% 
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
    )
  
  # Create labels for p-value and cell information
  Cell_Information = "Y"
  if (Cell_Information == "Y") {
    Cell_Number_Label <- c("Total Number of Cells treated with DMSO =\n")
    Temp <- Cell_Summary_by_Image %>% filter(SHORT_LABEL == "DMSO")
    Cell_Number_Label <- paste0(Cell_Number_Label, " ", Temp$MAX_CELL_NUMBER_BY_REPLICATE[1])
    for (index in 2:nrow(Temp)) {
      Cell_Number_Label <- paste0(Cell_Number_Label, " /", Temp$MAX_CELL_NUMBER_BY_REPLICATE[index])
    }
    Cell_Number_Label <- paste0(Cell_Number_Label, " ==> Total =", sum(Temp$MAX_CELL_NUMBER_BY_REPLICATE), "\nTotal Number of Cells treated with inhibitor =\n")
    Temp <- Cell_Summary_by_Image %>% filter(SHORT_LABEL != "DMSO")
    Cell_Number_Label <- paste0(Cell_Number_Label, " ", Temp$MAX_CELL_NUMBER_BY_REPLICATE[1])
    for (index in 2:nrow(Temp)) {
      Cell_Number_Label <- paste0(Cell_Number_Label, " /", Temp$MAX_CELL_NUMBER_BY_REPLICATE[index])
    }
    Cell_Number_Label <- paste0(Cell_Number_Label, " ==> Total =", sum(Temp$MAX_CELL_NUMBER_BY_REPLICATE))
  }
  rm(Temp, Cell_Information)
  
  # Calculate p-value for the difference in MEAN_NUMBER_OF_EVENTS between cohorts
  p_value_Result <- Cell_Summary_by_Image
  p_value_Result_1 <- wilcox.test(data = p_value_Result, MEAN_NUMBER_OF_EVENTS ~ SHORT_LABEL)$p.value
  p_value_Result_1 <- signif(p_value_Result_1, digits = 3)
  rm(p_value_Result)
  
  # P value calculation
  p_value_Result <- Cell_Summary_by_Image
  
  p_value_Result_1 <- wilcox.test(
    data = p_value_Result,
    MEAN_NUMBER_OF_EVENTS ~ SHORT_LABEL
  )$p.value
  
  p_value_Result_1 <- signif(p_value_Result_1, digits = 3)
  
  rm(p_value_Result) 
  
  return(
    list(
      Cell_Summary_by_Cell = Cell_Summary_by_Cell,
      Cell_Summary_by_Image = Cell_Summary_by_Image,
      Cell_Summary_by_Cohort = Cell_Summary_by_Cohort,
      Cell_Number_Label = Cell_Number_Label,
      p_value_Result_1 = p_value_Result_1
    )
  )
}


Lifetime_Events_Data_Processing <- function(Cell_Summary_by_Track){
  Cell_Summary_by_Track <- Cell_Summary_by_Track %>% 
    group_by(
      UNIVERSAL_TRACK_ID,
      COHORT,
      SHORT_LABEL
    ) %>% 
    summarise(
      LIFETIME = mean(LIFETIME)  # Calculate mean number of events
    ) %>% 
    as.data.table()
  
  
  Cell_Summary_by_Bin <- Cell_Summary_by_Track %>% 
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
  
  
  Cell_Summary_by_Image <- Cell_Summary_by_Track %>% 
    mutate(
      IMAGE = sub("(_\\d+\\.+.*)$", "", UNIVERSAL_TRACK_ID)
    ) %>% 
    group_by(
      IMAGE,
      SHORT_LABEL
    ) %>% 
    summarise(
      MEAN_LIFETIME = mean(LIFETIME)
    ) %>% 
    transform(
      ID_COHORT = as.numeric(factor(SHORT_LABEL))
    ) %>% 
    group_by(
      ID_COHORT
    ) %>% 
    mutate(
      ID_IMAGE = as.numeric(factor(IMAGE)),
      TOTAL_NUMBER_OF_IMAGES = max(ID_IMAGE)
    ) %>% 
    distinct(
      SHORT_LABEL,
      TOTAL_NUMBER_OF_IMAGES
    ) %>% 
    select(
      SHORT_LABEL,
      TOTAL_NUMBER_OF_IMAGES
    )
  
  
  Cell_Summary_by_Cohort <- Cell_Summary_by_Track %>% 
    group_by(
      SHORT_LABEL
    ) %>% 
    summarise(
      SD = sd(LIFETIME),
      MEAN_LIFETIME = mean(LIFETIME)
    )
  
  Cell_Summary_by_Cohort <- left_join(Cell_Summary_by_Cohort, Cell_Summary_by_Image, by = "SHORT_LABEL")
  
  Cell_Summary_by_Cohort <- Cell_Summary_by_Cohort %>% 
    mutate(
      SEM = SD/ TOTAL_NUMBER_OF_IMAGES
    )
  
  Label <- paste0(
    Cell_Summary_by_Cohort$SHORT_LABEL[1], "\nMean Lifetime +/- SEM = ", signif(Cell_Summary_by_Cohort$MEAN_LIFETIME[1], digits = 4), " +/- ",signif(Cell_Summary_by_Cohort$SEM[1], digits = 4),"\n",
    Cell_Summary_by_Cohort$SHORT_LABEL[2], "\nMean Lifetime = ", signif(Cell_Summary_by_Cohort$MEAN_LIFETIME[2], digits = 4), " +/- ",signif(Cell_Summary_by_Cohort$SEM[1], digits = 4)
  )
  
  return(
    list(
      Cell_Summary_by_Track = Cell_Summary_by_Track,
      Cell_Summary_by_Bin = Cell_Summary_by_Bin,
      Label = Label
    )
  )
}