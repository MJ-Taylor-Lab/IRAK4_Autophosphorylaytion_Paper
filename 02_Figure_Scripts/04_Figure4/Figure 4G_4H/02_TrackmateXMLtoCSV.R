# Load necessary libraries using pacman for easy installation and loading
pacman::p_load(dplyr, data.table, purrr, xml2, XML, parallel)

# Define the path to the CSV file and load it into a data.table
CSV_PATH <- "/Users/u_niranjan/Desktop/Image_Analysis_Grid/IRAK2KinaseDomainOnly/Cell_List.csv"
Table <- fread(CSV_PATH, header = TRUE)  # fread is used for fast data loading

# Define a function to generate paths for Trackmate files based on the table index
Trackmate_Path <- function(Table_index){
  Temp_Tables_List <- list()  # Initialize an empty list to store modified tables
  
  # Loop through each cell in the table to generate TrackMate paths
  for(Cell_index in 1:Table$NUMBER_OF_CELLS[Table_index]){
    Temp_Table <- Table[Table_index]  # Copy the row for manipulation
    # Add TrackMate_Path and CELL_NUMBER columns to Temp_Table
    Temp_Table <- Temp_Table %>% 
      mutate(
        TrackMate_Path = file.path(PATH, paste0("Cell", Cell_index), "Cell"),
        CELL_NUMBER = Cell_index,
        TrackMate_Path = paste0(TrackMate_Path, CELL_NUMBER, ".xml")
      )
    
    Temp_Tables_List[[Cell_index]] <- Temp_Table  # Append modified table to the list
  }
  
  Temp_Table_Final <- rbindlist(Temp_Tables_List)  # Combine all modified tables into one
  return(Temp_Table_Final)
}

# Apply the Trackmate_Path function to each row in Table and combine the results
Table_Final <- lapply(1:nrow(Table), Trackmate_Path)
Table_Final <- rbindlist(Table_Final)
Table_Final <- Table_Final %>% 
  arrange(COHORT)  # Arrange the final table by COHORT

# Clean up the environment by removing intermediate variables
rm(Table, Trackmate_Path)

XMLtoCSV <- function(XML_index){
  print(paste0(":::::::Starting Image ", XML_index, ":::::::"))
  Input_XML <- Table_Final$TrackMate_Path[XML_index]  # Define the input XML path
  
  # Create the output file path
  Output_Path <- strsplit(Input_XML, split = "/")[[1]]
  Output_Directory <- Output_Path[1]
  for(index in 2:(length(Output_Path)-1)){
    Output_Directory <- file.path(Output_Directory, Output_Path[index])
  }
  Output_Directory <- file.path(Output_Directory, "TrackmateData.csv")
  rm(index)
  
  if(file.exists(Input_XML) == TRUE){
    xml_data <- read_xml(Input_XML)  # Read the XML file
    
    # Extract spots and tracks from the XML data
    spots <- xml_find_all(xml_data, ".//Spot")
    tracks <- xml_find_all(xml_data, ".//Track")
    
    if(length(tracks) != 0){
      # Process each track to extract information
      Track_Function <- function(index){
        track_xml <- tracks[[index]]
        track_name <- xml_attr(track_xml, "name")
        track_id <- xml_attr(track_xml, "TRACK_ID")
        spot_source_id <- xml_attr(xml_find_first(track_xml, ".//Edge"), "SPOT_SOURCE_ID")
        spot_target_id <- xml_attr(xml_find_first(track_xml, ".//Edge"), "SPOT_TARGET_ID")
        track_duration <- xml_attr(track_xml, "TRACK_DURATION")
        
        # Create a data frame with track information
        track_info_df <- data.frame(
          TRACK_NAME = track_name,
          TRACK_ID = track_id,
          TRACK_DURATION = track_duration,
          SPOT_TARGET_ID = as.numeric(spot_target_id),
          SPOT_SOURCE_ID = as.numeric(spot_source_id),
          stringsAsFactors = FALSE
        )
        
        return(track_info_df)
      }
      
      track_info <- lapply(1:length(tracks), Track_Function)
      track_info <- rbindlist(track_info)
      
      track_info <- track_info %>% 
        mutate(
          SPOT_SOURCE_ID_MIN = pmin(SPOT_TARGET_ID, SPOT_SOURCE_ID),
          SPOT_TARGET_ID_MAX = pmax(SPOT_TARGET_ID, SPOT_SOURCE_ID)
        ) %>% 
        select(
          -c("SPOT_TARGET_ID", "SPOT_SOURCE_ID")
        )
      
      # Expand the table based on the DIFFERENCE column
      track_info <- track_info[, .(SPOT = seq(SPOT_SOURCE_ID_MIN, SPOT_TARGET_ID_MAX - 1)),
                               by = .(TRACK_NAME, TRACK_ID, TRACK_DURATION, SPOT_SOURCE_ID_MIN, SPOT_TARGET_ID_MAX)]
      
      # Add the SPOT_TARGET_ID_MAX to the SPOT column to complete the sequence
      track_info <- rbind(track_info, track_info[, .(SPOT = SPOT_TARGET_ID_MAX), by = .(TRACK_NAME, TRACK_ID, TRACK_DURATION, SPOT_SOURCE_ID_MIN, SPOT_TARGET_ID_MAX)])
      track_info <- track_info %>% 
        arrange(
          SPOT
        ) %>% 
        distinct(
          SPOT,
          .keep_all = TRUE
        )
      
      # Initialize an empty data frame with the specified column names
      spot_data <- data.frame(
        LABEL = character(),
        ID = integer(),
        QUALITY = numeric(),
        POSITION_X = numeric(),
        POSITION_Y = numeric(),
        POSITION_Z = numeric(),
        POSITION_T = numeric(),
        FRAME = integer(),
        RADIUS = numeric(),
        VISIBILITY = integer(),
        MEAN_INTENSITY_CH1 = numeric(),
        MEDIAN_INTENSITY_CH1 = numeric(),
        MIN_INTENSITY_CH1 = numeric(),
        MAX_INTENSITY_CH1 = numeric(),
        TOTAL_INTENSITY_CH1 = numeric(),
        STD_INTENSITY_CH1 = numeric(),
        CONTRAST_CH1 = numeric(),
        SNR_CH1 = numeric(),
        stringsAsFactors = FALSE
      )
      
      for (spot in spots) {
        spot_data <- rbind(
          spot_data,
          data.frame(
            LABEL = xml_attr(spot, "name"),
            ID = as.integer(xml_attr(spot, "ID")),
            QUALITY = as.numeric(xml_attr(spot, "QUALITY")),
            POSITION_X = as.numeric(xml_attr(spot, "POSITION_X")),
            POSITION_Y = as.numeric(xml_attr(spot, "POSITION_Y")),
            POSITION_Z = as.numeric(xml_attr(spot, "POSITION_Z")),
            POSITION_T = as.numeric(xml_attr(spot, "POSITION_T")),
            FRAME = as.integer(xml_attr(spot, "FRAME")),
            RADIUS = as.numeric(xml_attr(spot, "RADIUS")),
            VISIBILITY = as.integer(xml_attr(spot, "VISIBILITY")),
            MEAN_INTENSITY_CH1 = as.numeric(xml_attr(spot, "MEAN_INTENSITY_CH1")),
            MEDIAN_INTENSITY_CH1 = as.numeric(xml_attr(spot, "MEDIAN_INTENSITY_CH1")),
            MIN_INTENSITY_CH1 = as.numeric(xml_attr(spot, "MIN_INTENSITY_CH1")),
            MAX_INTENSITY_CH1 = as.numeric(xml_attr(spot, "MAX_INTENSITY_CH1")),
            TOTAL_INTENSITY_CH1 = as.numeric(xml_attr(spot, "TOTAL_INTENSITY_CH1")),
            STD_INTENSITY_CH1 = as.numeric(xml_attr(spot, "STD_INTENSITY_CH1")),
            CONTRAST_CH1 = as.numeric(xml_attr(spot, "CONTRAST_CH1")),
            SNR_CH1 = as.numeric(xml_attr(spot, "SNR_CH1"))
          )
        )
        rm(spot)
      }
      
      spot_data <- spot_data %>% 
        mutate(
          TRACK_ID_index = ID %in% track_info$SPOT
        ) %>% 
        filter(
          TRACK_ID_index == TRUE
        ) %>% 
        rowwise() %>%  # Operate on each row individually
        mutate(
          TRACK_ID_index = list(which(ID == track_info$SPOT)),  # Use list to ensure the result is treated correctly
          TRACK_NAME = track_info$TRACK_NAME[TRACK_ID_index],
          TRACK_ID = track_info$TRACK_ID[TRACK_ID_index],
          TRACK_DURATION = track_info$TRACK_DURATION[TRACK_ID_index]
        ) %>% 
        select(
          -c(TRACK_ID_index)
        )
      
      # Write csv
      write.csv(spot_data, Output_Directory)
      
    } else{
      # Initialize an empty data frame with the specified column names
      spot_data <- data.frame(
        LABEL = character(),
        ID = integer(),
        QUALITY = numeric(),
        POSITION_X = numeric(),
        POSITION_Y = numeric(),
        POSITION_Z = numeric(),
        POSITION_T = numeric(),
        FRAME = integer(),
        RADIUS = numeric(),
        VISIBILITY = integer(),
        MEAN_INTENSITY_CH1 = numeric(),
        MEDIAN_INTENSITY_CH1 = numeric(),
        MIN_INTENSITY_CH1 = numeric(),
        MAX_INTENSITY_CH1 = numeric(),
        TOTAL_INTENSITY_CH1 = numeric(),
        STD_INTENSITY_CH1 = numeric(),
        CONTRAST_CH1 = numeric(),
        SNR_CH1 = numeric(),
        TRACK_NAME = character(),
        TRACK_ID = character(),
        TRACK_DURATION = character(),
        stringsAsFactors = FALSE
      )
      
      # Write csv
      write.csv(spot_data, Output_Directory)
    }
  } 
  else {
    # Initialize an empty data frame with the specified column names
    spot_data <- data.frame(
      LABEL = character(),
      ID = integer(),
      QUALITY = numeric(),
      POSITION_X = numeric(),
      POSITION_Y = numeric(),
      POSITION_Z = numeric(),
      POSITION_T = numeric(),
      FRAME = integer(),
      RADIUS = numeric(),
      VISIBILITY = integer(),
      MEAN_INTENSITY_CH1 = numeric(),
      MEDIAN_INTENSITY_CH1 = numeric(),
      MIN_INTENSITY_CH1 = numeric(),
      MAX_INTENSITY_CH1 = numeric(),
      TOTAL_INTENSITY_CH1 = numeric(),
      STD_INTENSITY_CH1 = numeric(),
      CONTRAST_CH1 = numeric(),
      SNR_CH1 = numeric(),
      TRACK_NAME = character(),
      TRACK_ID = character(),
      TRACK_DURATION = character(),
      stringsAsFactors = FALSE
    )
    
    # Write csv
    write.csv(spot_data, Output_Directory)
  }
  
  Output_Directory <- as.data.table(Output_Directory)
  return(Output_Directory)
}  


# Define the number of cores to use
numCores <- detectCores() - 3  # Reserve two core for system processes
# Use mclapply to parallelize XMLtoCSV function over rows of Table_Final
Trackmate_csv_Path <- mclapply(1:nrow(Table_Final), XMLtoCSV, mc.cores = numCores)
Trackmate_csv_Path <- rbindlist(Trackmate_csv_Path)
Trackmate_csv_Path <- Trackmate_csv_Path %>% as.data.table

# Define the output directory for the final CSV file
Output_Path <- strsplit(CSV_PATH, split = "/")[[1]]
Output_Directory <- c("/")
for(index in 2:(length(Output_Path)-1)){
  Output_Directory <- file.path(Output_Directory, Output_Path[index])
}
Output_Directory <- file.path(Output_Directory, "Trackmate_List.csv")

# Write the final CSV file
write.csv(Trackmate_csv_Path, Output_Directory)

# Clean up the environment by removing all objects
rm(list = ls())
gc()
