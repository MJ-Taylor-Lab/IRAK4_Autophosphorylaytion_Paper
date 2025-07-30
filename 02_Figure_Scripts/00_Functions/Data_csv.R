Data_Dwell_Time <- function(
    Cell_Summary_by_Track, 
    Cell_Summary_by_Cohort, 
    Plot_Script_Directory, 
    Plot_Directory_Save_Path,
    Save_Name
  ){
  
  # Writing Source Data
  Source_Data_Path <- paste0(basename(Plot_Directory_Save_Path), " ",  basename(Plot_Script_Directory), ".csv")
  Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Source_Data_Path)
  write.csv(Cell_Summary_by_Track, Plot_Save_Path)
  
  #Writing Legend Data
  Legend_Data <- Cell_Summary_by_Cohort %>% 
    select(
      SHORT_LABEL,
      COHORT,
      IMAGE_COUNT,
      FRACTION,
      SEM,
      MIN_NUMBER_OF_CELLS_PER_REPLICATE,
      MAX_NUMBER_OF_CELLS_PER_REPLICATE,
      TOTAL_NUMBER_OF_CELLS_PER_COHORT,
      PRIMARY_PROTIEN_ASSOCIATION_EVENT_COUNT_byCOHORT,
      PRIMARY_PROTIEN_TRACK_COUNT_byCOHORT
    )
  
  # Writing Source Data
  Source_Data_Path <- paste0(basename(Plot_Directory_Save_Path), " ",  basename(Plot_Script_Directory), " Legend Summary.csv")
  Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Source_Data_Path)
  write.csv(Legend_Data, Plot_Save_Path)
  
  rm(
    Legend_Data,
    Source_Data_Path
  )
  
}


Data_Colocalisation_Percentage <- function(
    Colocalisation_Percentage_byCell, 
    Colocalisation_Percentage_byCOHORT,
    Plot_Script_Directory, 
    Plot_Directory_Save_Path,
    Save_Name
  ){
  
  # Writing Source Data
  Source_Data_Path <- paste0(Save_Name, " ",  basename(Plot_Script_Directory), ".csv")
  Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Source_Data_Path)
  write.csv(Colocalisation_Percentage_byCell, Plot_Save_Path)
  
  
  
  #Writing Legend Data
  Legend_Data <- Colocalisation_Percentage_byCOHORT %>% 
    dplyr::select(
      SHORT_LABEL,
      COHORT,
      MEAN_COLOCLIZED_SPOT_TEST,
      STANDARD_ERROR_Of_MEAN_SPOT_TEST,
      TOTAL_NUMBER_OF_IMAGES,
      TOTAL_NUMBER_OF_CELLS_PER_COHORT,
      MIN_NUMBER_OF_CELLS_PER_REPLICATE,
      MAX_NUMBER_OF_CELLS_PER_REPLICATE
    )
  
  # Writing Source Data
  Source_Data_Path <- paste0(Save_Name, " ",  basename(Plot_Script_Directory), " Legend Summary.csv")
  Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Source_Data_Path)
  write.csv(Legend_Data, Plot_Save_Path)
  
  rm(
    Legend_Data,
    Source_Data_Path
  )
}