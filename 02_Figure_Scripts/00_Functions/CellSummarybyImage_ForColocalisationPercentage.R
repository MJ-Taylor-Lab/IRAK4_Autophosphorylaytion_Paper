# Define a function to summarize cell tracking data for colocalization percentage analysis
Cell_Summary_by_Track_for_ColocalisationPercentage <- function(Table) {
  # Filter the input table for specific conditions
  Cell_Summary_by_Track <- Table %>%
    group_by(
      UNIVERSAL_TRACK_ID                    # Group by universal track ID
    ) %>%
    arrange(
      FRAME,                                # Arrange data by frame within each group
      .by_group = TRUE
    ) %>%
    group_by(
      UNIVERSAL_TRACK_ID                    # Group by universal track ID again (to ensure group consistency)
    ) %>%
    mutate(
      LIFETIME = max(TIME_ADJUSTED) - min(TIME_ADJUSTED)  # Calculate lifetime of each track
    ) %>%
    mutate(
      COLOCALIZATION = COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= COLOCALISATION_THREHOLD # Determine colocalization based on a threshold
    ) %>%
    group_by(
      UNIVERSAL_TRACK_ID                    # Group by universal track ID again
    ) %>%
    mutate(
      STREAK = cumsum(!COLOCALIZATION),     # Create a grouping variable for continuous frames above threshold
      DWELL_TIME = lead(TIME_ADJUSTED, default = last(TIME_ADJUSTED)) -  TIME_ADJUSTED # Calculate dwell time
    ) %>%
    filter(
      COLOCALIZATION == TRUE                # Only keep rows where colocalization is TRUE
    ) %>%
    group_by(
      COHORT,                               # Group by cohort
      SHORT_LABEL,                          # Group by short label
      IMAGE,                                # Group by image
      UNIVERSAL_TRACK_ID,                   # Group by universal track ID
      STREAK                                # Group continuous frames above threshold together
    ) %>%
    summarise(
      DWELL_FRAMES = sum(COLOCALIZATION),   # Sum of frames where colocalization is TRUE
      DWELL_TIME = sum(DWELL_TIME),         # Sum of dwell times
      MAX_NORMALIZED_INTENSITY = max(NORMALIZED_INTENSITY), # Maximum normalized intensity
      LIFETIME = max(LIFETIME)              # Maximum lifetime
    ) %>%
    filter(
      DWELL_FRAMES >= 2                     # Filter for dwell frames >= 2 (transient recruitment events)
    ) %>%
    mutate(
      DWELL_TIME_BIN = round(DWELL_TIME)    # Round dwell time to nearest integer
    ) %>%
    ungroup() %>% 
    arrange(
      SHORT_LABEL,
      UNIVERSAL_TRACK_ID                    # Arrange by universal track ID
    ) %>%
    as.data.table()                         # Convert to data.table
  
  return(Cell_Summary_by_Track)              # Return the resulting data table
}


#Finding the Colocalisation_Percentage_byCell, Image and Cohort based on Dwell Time Filter
Colocalisation_Percentage_byCell <- function(Table, Cell_Summary_by_Track, Count){
  # Filter the Cell_Summary_by_Track data frame and create a distinct list of universal track IDs based on dwell frame
  Colocalisation_Percentage_List <- Cell_Summary_by_Track %>% 
    filter(
      DWELL_FRAMES >= Dwell_Test_List$DWELL_FRAMES[Count] # Filter rows where dwell frames are greater than or equal to a specified threshold
    ) %>%
    distinct(
      UNIVERSAL_TRACK_ID, # Select distinct universal track IDs
      SHORT_LABEL
    )
  
  
  # Calculate the percentage of colocalization per cell based on specific conditions
  Colocalisation_Percentage_byCell <- Table %>%
    # Select distinct universal track IDs, keeping all columns
    distinct(
      UNIVERSAL_TRACK_ID,
      .keep_all = TRUE                         # Retain all columns while selecting distinct UNIVERSAL_TRACK_IDs
    ) %>%
    # Add a new column indicating if the track ID is in the colocalisation list
    mutate(
      COLOCLIZED_SPOT_TEST = UNIVERSAL_TRACK_ID %in% Colocalisation_Percentage_List$UNIVERSAL_TRACK_ID,
      # Assign 100 if the spot is colocalized, otherwise assign 0
      COLOCLIZED_SPOT_TEST = case_when(
        COLOCLIZED_SPOT_TEST == TRUE ~ 100,
        COLOCLIZED_SPOT_TEST != TRUE ~ 0
      )
    ) %>%
    # Arrange rows by UNIVERSAL_TRACK_ID
    arrange(
      UNIVERSAL_TRACK_ID,
      SHORT_LABEL
    ) %>%
    # Group by multiple columns to summarize at the cell level
    group_by(
      IMAGE,
      COHORT,
      SHORT_LABEL,
      CELL
    ) %>%
    # Calculate the mean colocalization percentage per cell
    summarise(
      MEAN_COLOCLIZED_SPOT_TEST = mean(COLOCLIZED_SPOT_TEST)   # Find the mean percentage of colocalization per cell
    ) %>%
    # Further group by IMAGE, COHORT, and SHORT_LABEL
    group_by(
      IMAGE,
      COHORT,
      SHORT_LABEL
    ) %>%
    # Extract and convert the date from the IMAGE column
    mutate(
      DATE = strsplit(IMAGE, " ", fixed = TRUE)[[1]][1],
      DATE = as.integer(DATE)
    ) %>%
    # Arrange rows by DATE, COHORT, and IMAGE
    arrange(
      DATE,
      COHORT,
      IMAGE
    ) %>%
    # Transform the date into a numeric factor for ID
    transform(
      ID_DATE = as.numeric(factor(DATE))
    ) %>%
    # Group by DATE to assign colors based on date ID
    group_by(
      DATE
    ) %>%
    # Assign colors based on the ID_DATE value
    mutate(
      COLOR_DATE = case_when(
        ID_DATE == 1 ~ "#7CAE00",
        ID_DATE == 2 ~ "#CD9600",
        ID_DATE == 3 ~ "#F8766D",
        ID_DATE == 4 ~ "#00BE67",
        ID_DATE == 5 ~ "#00BFC4",
        ID_DATE == 6 ~ "#00A9FF",
        ID_DATE == 7 ~ "#C77CFF",
        ID_DATE == 8 ~ "#FF61CC"
      )
    ) %>%
    # Ungroup the data for further operations
    ungroup() %>%
    # Arrange rows by SHORT_LABEL
    arrange(
      SHORT_LABEL
    ) %>%
    # Convert the resulting data frame to a data table
    as.data.table()
  
  
  
  # Colocalisation by Image
  Colocalisation_Percentage_byImage <- Colocalisation_Percentage_byCell %>% 
    group_by(
      IMAGE,
      COHORT,
      SHORT_LABEL
    ) %>% 
    summarise(
      MEAN_COLOCLIZED_SPOT_TEST = mean(MEAN_COLOCLIZED_SPOT_TEST)
    ) %>% 
    arrange(
      SHORT_LABEL
    ) %>% 
    as.data.table()
  
  
  #Calculating Max and Min Number of Cells in replicates to get range
  NUMBER_OF_CELLS_PER_REPLICATE_COUNT <- unique(Colocalisation_Percentage_byCell[, .(IMAGE, SHORT_LABEL, CELL)])
  NUMBER_OF_CELLS_PER_REPLICATE_COUNT <- NUMBER_OF_CELLS_PER_REPLICATE_COUNT %>% 
    group_by(
      IMAGE,
      SHORT_LABEL
    ) %>% 
    summarise(
      NUMBER_OF_CELLS_PER_REPLICATE = n()
    ) %>% 
    group_by(
      SHORT_LABEL
    ) %>% 
    summarise(
      MIN_NUMBER_OF_CELLS_PER_REPLICATE = min(NUMBER_OF_CELLS_PER_REPLICATE),
      MAX_NUMBER_OF_CELLS_PER_REPLICATE = max(NUMBER_OF_CELLS_PER_REPLICATE),
      TOTAL_NUMBER_OF_CELLS_PER_COHORT = sum(NUMBER_OF_CELLS_PER_REPLICATE)
    )
  
  # COlocalisation Percentage by Cohort
  Colocalisation_Percentage_byCOHORT <- Colocalisation_Percentage_byCell %>% 
    arrange(
      SHORT_LABEL,
      IMAGE
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
    ungroup() %>% 
    group_by(
      IMAGE,
      COHORT,
      SHORT_LABEL,
      TOTAL_NUMBER_OF_IMAGES
    ) %>% 
    summarise(
      MEAN_COLOCLIZED_SPOT_TEST = mean(MEAN_COLOCLIZED_SPOT_TEST)
    ) %>% 
    ungroup() %>% 
    group_by(
      SHORT_LABEL,
      COHORT
    ) %>%
    summarise(
      TOTAL_NUMBER_OF_IMAGES = mean(TOTAL_NUMBER_OF_IMAGES),
      STANDARD_DEVIATION_OF_MEAN_COLOCLIZED_SPOT_TEST = sd(MEAN_COLOCLIZED_SPOT_TEST),
      STANDARD_ERROR_Of_MEAN_SPOT_TEST = STANDARD_DEVIATION_OF_MEAN_COLOCLIZED_SPOT_TEST/TOTAL_NUMBER_OF_IMAGES,
      MEAN_COLOCLIZED_SPOT_TEST = mean(MEAN_COLOCLIZED_SPOT_TEST),
      YMAX_COLOCLIZED_SPOT_TEST = MEAN_COLOCLIZED_SPOT_TEST + STANDARD_ERROR_Of_MEAN_SPOT_TEST,
      YMIN_COLOCLIZED_SPOT_TEST = case_when(
        MEAN_COLOCLIZED_SPOT_TEST - STANDARD_ERROR_Of_MEAN_SPOT_TEST < 0 ~ 0,
        TRUE ~ MEAN_COLOCLIZED_SPOT_TEST - STANDARD_ERROR_Of_MEAN_SPOT_TEST
      )
    )
  
  Colocalisation_Percentage_byCOHORT <- inner_join(Colocalisation_Percentage_byCOHORT, NUMBER_OF_CELLS_PER_REPLICATE_COUNT, by = "SHORT_LABEL") 
  
  
  
  return(
    list(
      Colocalisation_Percentage_byCell = Colocalisation_Percentage_byCell, 
      Colocalisation_Percentage_byImage = Colocalisation_Percentage_byImage,
      Colocalisation_Percentage_byCOHORT = Colocalisation_Percentage_byCOHORT
    )
  )
  
  rm(Colocalisation_Percentage_List)
}