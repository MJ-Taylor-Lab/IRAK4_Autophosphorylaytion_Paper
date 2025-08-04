Colocalisation_Percentage_byCellFx <- function(Colocalisation_Percentage_byCell){
  
  Colocalisation_Percentage_byCell <- Colocalisation_Percentage_byCell %>% as.data.table()
  
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
}