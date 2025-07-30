Size_Summary_ByImage <- function(
    Table,
    SIZE_THRESHOLD
  ){
  
  # Finding number of spots per cell that are greater or less than size of Threshold of 4
  Size_Summary_by_Track <- Table %>% 
    #Remove rows with missing complementary protien data
    filter(
      !is.na(COMPLEMENTARY_NORMALIZED_INTENSITY_1)
    ) %>% 
    group_by(
      UNIVERSAL_TRACK_ID
    ) %>% 
    # Selecting for IRAK4 tracks that recruit even a little bit of complementary protein
    filter(
      COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 0.75
    ) %>% 
    group_by(
      UNIVERSAL_TRACK_ID,
      CELL,
      IMAGE,
      SHORT_LABEL
    ) %>% 
    # Max Size of Complementary protein recruited
    summarise(
      MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1 = max(COMPLEMENTARY_NORMALIZED_INTENSITY_1)
    ) %>% 
    mutate(
      SIZE_CATEGORY = case_when(
        MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1 < SIZE_THRESHOLD ~ 0,
        TRUE ~ 100
      )
    ) %>% 
    as.data.table()
  
  #Calculating Max and Min Number of Cells in replicates to get range
  NUMBER_OF_CELLS_PER_REPLICATE_COUNT <- unique(Size_Summary_by_Track[, .(IMAGE, SHORT_LABEL, CELL, SIZE_CATEGORY)])
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
  
  #Caluclating number of association events and tracks.
  UNIQUE_PRIMARY_PROTIEN_COUNT <- Size_Summary_by_Track %>% 
    group_by(
      SHORT_LABEL
    ) %>% 
    mutate(
      PRIMARY_PROTIEN_ASSOCIATION_EVENT_COUNT = n()
    ) %>% 
    distinct(
      UNIVERSAL_TRACK_ID,
      PRIMARY_PROTIEN_ASSOCIATION_EVENT_COUNT
    ) %>% 
    mutate(
      PRIMARY_PROTIEN_TRACK_COUNT = n()
    ) %>% 
    summarise(
      PRIMARY_PROTIEN_ASSOCIATION_EVENT_COUNT_byCOHORT = max(PRIMARY_PROTIEN_ASSOCIATION_EVENT_COUNT),
      PRIMARY_PROTIEN_TRACK_COUNT_byCOHORT = max(PRIMARY_PROTIEN_TRACK_COUNT)
    )
  
  
  UNIQUE_PRIMARY_PROTIEN_COUNT <- inner_join(NUMBER_OF_CELLS_PER_REPLICATE_COUNT, UNIQUE_PRIMARY_PROTIEN_COUNT, by = "SHORT_LABEL") # For Each Cohort we get the total number of PRIMARY_PROTIEN tracks and Association Events
  
  Size_Summary_by_Cell <- Size_Summary_by_Track %>% 
    group_by(
      CELL,
      IMAGE,
      SHORT_LABEL
    ) %>% 
    summarise(
      PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD = mean(SIZE_CATEGORY)
    ) %>% 
    arrange(
      SHORT_LABEL,
      IMAGE,
      CELL
    ) %>%
    as.data.table()
  
  rm(Size_Summary_by_Track)
  
  # Summary statistics by image
  Size_Summary_ByImage <- Size_Summary_by_Cell %>% 
    group_by(
      IMAGE,
      SHORT_LABEL
    ) %>% 
    summarise(
      MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD = mean(PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD)
    ) %>% 
    arrange(
      SHORT_LABEL,
      IMAGE
    )
  
  # Summary statistics by cohort
  Size_Summary_ByCohort <- Size_Summary_ByImage %>% 
    arrange(
      SHORT_LABEL,
      IMAGE
    ) %>% 
    # Counting number of replicates per cohort
    transform(
      ID_COHORT = as.numeric(factor(SHORT_LABEL))
    ) %>% 
    group_by(
      ID_COHORT
    ) %>% 
    # Assigning highest number as number of replicates per cohort
    mutate(
      ID_IMAGE = as.numeric(factor(IMAGE)),
      TOTAL_NUMBER_OF_IMAGES = max(ID_IMAGE)
    ) %>% 
    ungroup() %>% 
    group_by(
      SHORT_LABEL
    ) %>% 
    summarise(
      TOTAL_NUMBER_OF_IMAGES = mean(TOTAL_NUMBER_OF_IMAGES),
      #Finding SD, SEM, MEAN, YMAX and YMIN
      SD = sd(MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD),
      SEM = SD/TOTAL_NUMBER_OF_IMAGES,
      MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD = mean(MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD),
      Y_MAX = MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD + SEM,
      Y_MIN = case_when(
        MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD - SEM < 0 ~ 0,
        TRUE ~ MEAN_OF_MEAN_PERCENTAGE_OF_TRACKS_LARGER_THAN_THRESHOLD - SEM
      )
    )
  
  Size_Summary_ByCohort <- inner_join(Size_Summary_ByCohort, UNIQUE_PRIMARY_PROTIEN_COUNT, by = "SHORT_LABEL") # For Each Cohort we get the total number of PRIMARY_PROTIEN tracks and Association Events
  
  return(
    list(
      Size_Summary_ByImage = Size_Summary_ByImage,
      Size_Summary_ByCohort = Size_Summary_ByCohort
    )
  )
}

