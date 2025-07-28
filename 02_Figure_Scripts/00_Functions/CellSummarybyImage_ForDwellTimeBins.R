Cell_Summary_by_Image <- function(Table){
  Cell_Summary_by_Image <- Table %>%
    group_by(
      UNIVERSAL_TRACK_ID
    ) %>%  # Group data by universal track ID
    arrange(
      FRAME, 
      .by_group = TRUE
    ) %>%  # Arrange rows by frame within each group
    mutate(
      LIFETIME = max(TIME_ADJUSTED) - min(TIME_ADJUSTED)  # Calculate lifetime of each track
    ) %>%
    mutate(
      COLOCALIZATION = COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= COLOCALISATION_THREHOLD  # Determine if colocalization occurs
    ) %>%
    mutate(
      STREAK = cumsum(!COLOCALIZATION),  # Calculate streaks of colocalization
      DWELL_TIME = lead(TIME_ADJUSTED, default = last(TIME_ADJUSTED)) - TIME_ADJUSTED  # Calculate dwell time
    ) %>%
    filter(
      COLOCALIZATION == TRUE
    ) %>%  # Only keep rows where colocalization is TRUE
    group_by(
      COHORT, 
      SHORT_LABEL, 
      IMAGE, 
      UNIVERSAL_TRACK_ID, 
      STREAK,
      CELL
    ) %>%
    summarise(
      DWELL_FRAMES = sum(COLOCALIZATION),  # Sum of colocalization frames
      DWELL_TIME = sum(DWELL_TIME),  # Sum of dwell time
      MAX_NORMALIZED_INTENSITY = max(NORMALIZED_INTENSITY),  # Maximum normalized intensity
      LIFETIME = max(LIFETIME)  # Maximum lifetime
    ) %>%
    filter(
      DWELL_FRAMES >= 2
    ) %>%  # Filter rows with at least 2 dwell frames
    mutate(
      DWELL_TIME_BIN = round(DWELL_TIME)  # Round dwell time to the nearest integer
    ) %>%
    arrange(
      UNIVERSAL_TRACK_ID
    ) %>%
    as.data.table()
  
  #Calculating Max and Min Number of Cells in replicates to get range
  NUMBER_OF_CELLS_PER_REPLICATE_COUNT <- unique(Cell_Summary_by_Image[, .(IMAGE, SHORT_LABEL, CELL)])
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
  UNIQUE_PRIMARY_PROTIEN_COUNT <- Cell_Summary_by_Image %>% 
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
  
  
  Cell_Summary_by_Image <- Cell_Summary_by_Image %>% 
    mutate(
      BIN = (trunc(DWELL_TIME_BIN / 3)) * 3,  # Bin dwell time into intervals of 3
      BIN = case_when(
        BIN < 30 ~ BIN,
        BIN >= 30 ~ 33  # Cap BIN values at 33
      )
    ) %>%
    group_by(
      IMAGE, 
      SHORT_LABEL, 
      COHORT, 
      BIN
    ) %>%
    mutate(
      BIN_COUNT = n()  # Count the number of occurrences in each BIN
    ) %>%
    filter(
      BIN > 0
    ) %>%  # Filter rows with BIN values greater than 0
    distinct(
      BIN_COUNT, 
      BIN
    ) %>%
    arrange(
      IMAGE, 
      BIN, 
      SHORT_LABEL, 
      COHORT
    ) %>%
    ungroup() %>%
    group_by(
      IMAGE, 
      SHORT_LABEL, 
      COHORT
    ) %>%
    mutate(
      FRACTION = sum(BIN_COUNT)  # Calculate the fraction of BIN counts
    ) %>%
    group_by(
      IMAGE, 
      SHORT_LABEL, 
      COHORT, 
      BIN
    ) %>%
    mutate(
      FRACTION = BIN_COUNT / FRACTION  # Calculate the fraction within each BIN
    ) %>%
    filter(
      BIN > 30
    ) %>%  # Filter rows with BIN values greater than 30
    as.data.table()
  
  #Images that dont possess BINS above 0 or 30 are removed in the filter above and need to be restored with values set at 0
  Cell_Summary_by_Image_ImageList <- unique(Cell_Summary_by_Image$IMAGE) %>% as.data.table()
  Missing_Images <- unique(Table$IMAGE) %>% as.data.table()
  if(nrow(Cell_Summary_by_Image_ImageList) != nrow(Missing_Images)){
    colnames(Missing_Images) <- "IMAGE"
    Missing_Images <- Table %>% 
      dplyr::distinct(
        IMAGE,
        COHORT,
        .keep_all = TRUE
      ) %>% 
      dplyr::select(
        IMAGE, SHORT_LABEL, COHORT
      ) %>% 
      mutate(
        MISSING_IMAGE_TEST = IMAGE %in% unique(Cell_Summary_by_Image$IMAGE)
      ) %>% 
      dplyr::filter(
        MISSING_IMAGE_TEST == FALSE
      ) %>% 
      dplyr::select(
        -c(MISSING_IMAGE_TEST)
      ) %>% 
      mutate(
        BIN_COUNT = 0,
        BIN = 33,
        FRACTION = 0
      )
    
    Cell_Summary_by_Image <- rbind(Cell_Summary_by_Image, Missing_Images)
  }
  
  rm(
    Cell_Summary_by_Image_ImageList,
    Missing_Images
  )
  
 
  Cell_Summary_by_Image <- inner_join(Cell_Summary_by_Image, UNIQUE_PRIMARY_PROTIEN_COUNT, by = "SHORT_LABEL") # For Each Cohort we get the total number of PRIMARY_PROTIEN tracks and Association Events
  Cell_Summary_by_Image <- Cell_Summary_by_Image %>% 
    rowwise() %>% 
    mutate(
      FRACTION = FRACTION*100
    ) %>% 
    arrange(
      SHORT_LABEL
    )
  
  return(Cell_Summary_by_Image)
}