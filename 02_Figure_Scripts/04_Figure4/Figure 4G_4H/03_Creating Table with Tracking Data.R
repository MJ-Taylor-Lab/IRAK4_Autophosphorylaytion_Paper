# Load necessary libraries using pacman for easy package management
library(pacman)
pacman::p_load(dplyr, data.table, stringr)

# Load the cell list CSV file into a data.table called 'Table'
Trackmate_List_Directory <- "/Users/u_niranjan/Desktop/Image_Analysis_Grid/IRAK2KinaseDomainOnly/Trackmate_List.csv"
Table <- fread(Trackmate_List_Directory, header = TRUE)
Table <- Table %>% 
  select(
    Output_Directory
  )

Cell_List <- dirname(Trackmate_List_Directory)
Cell_List <- file.path(Cell_List, "Cell_List.csv")

Cell_List <- fread(Cell_List, header = TRUE)

# Define the folder path where the final table will be saved
Table_Save_Folder <- dirname(Trackmate_List_Directory)

# Function to process each Trackmate file and count spots
Spot_Counter <- function(TrackMate_Path){
  # Extracting Information from Final Table
  Table_information <- Table$Output_Directory[TrackMate_Path]
  
  # Read the TrackMate file
  Trackmate_File <-  fread(Table_information, header = TRUE)
  Trackmate_File <- Trackmate_File %>% 
    select(
      -c(1)
    )
  
  # Check if Track list is empty if so it will be set as 0
  if(nrow(Trackmate_File) != 0) {
    Trackmate_File <- Trackmate_File %>% 
      arrange(
        TRACK_ID, 
        FRAME
      ) %>% 
      group_by(
        TRACK_ID
      ) %>% 
      mutate(
        POSITION_T = as.numeric(POSITION_T)
      ) %>% 
      filter(
        !is.na(POSITION_T)
      ) %>% 
      group_by(
        TRACK_ID
      ) %>% 
      mutate(
        MAX_TIME = max(POSITION_T),
        MIN_TIME = min(POSITION_T)
      ) %>% 
      ungroup()
  } else {
    Trackmate_File <- data.table(
      LABEL = c(0),
      ID = c(0),
      QUALITY = c(0),
      POSITION_X = c(0),
      POSITION_Y = c(0),
      POSITION_Z = c(0),
      POSITION_T = c(0),
      FRAME = c(0),
      RADIUS = c(0),
      VISIBILITY = c(0),
      MEAN_INTENSITY_CH1 = c(0),
      MEDIAN_INTENSITY_CH1 = c(0),
      MIN_INTENSITY_CH1 = c(0),
      MAX_INTENSITY_CH1 = c(0),
      TOTAL_INTENSITY_CH1 = c(0),
      STD_INTENSITY_CH1 = c(0),
      CONTRAST_CH1 = c(0),
      SNR_CH1 = c(0),
      TRACK_NAME = "Track_0",
      TRACK_ID = c(0),
      TRACK_DURATION = c(0),
      MAX_TIME = c(0),
      MIN_TIME = c(0)
    )
  }
  
# Find Number of Distinct Events
Unique_Events = length(unique(Trackmate_File$TRACK_ID))   

# Adjusting Distinct Event COunt in case trackmate had 0 tracks
if(Unique_Events != 1){
  Unique_Events = Unique_Events
} else if(Unique_Events == 1 && Trackmate_File$MAX_TIME[1] == 0){
  Unique_Events = 0
} else {
  Unique_Events = Unique_Events
}

    
Trackmate_File <- Trackmate_File %>%     
    mutate(
      TIME_ADJUSTED = POSITION_T - MIN_TIME,
      LIFETIME = MAX_TIME - MIN_TIME,
      NUMBER_OF_EVENTS = Unique_Events
    ) %>% 
    group_by(
      TRACK_ID
    ) %>%
    mutate(
      NEW_ID = row_number(),
      CSV_PATH = Table_information,
      CSV_DIRECTORY = dirname(CSV_PATH),
      CELL_NUMBER = CSV_DIRECTORY,
      CELL_NUMBER = as.integer(sub(".*Cell(\\d+).*", "\\1", CELL_NUMBER)),
      CSV_DIRECTORY = dirname(CSV_DIRECTORY)
    ) %>% 
    select(
      -c(
        VISIBILITY,
        MEAN_INTENSITY_CH1,
        MEDIAN_INTENSITY_CH1,
        MIN_INTENSITY_CH1,
        MAX_INTENSITY_CH1,
        TOTAL_INTENSITY_CH1,
        STD_INTENSITY_CH1,
        CONTRAST_CH1,
        SNR_CH1
      )
    ) %>% 
  rowwise() %>% 
  mutate(
    INDEX = which(Cell_List$PATH == CSV_DIRECTORY),
    IMAGE = Cell_List$IMAGE[INDEX],
    CELL_LINE_NUMBER = Cell_List$CELL_LINE_NUMBER[INDEX],
    CELL_LINE = Cell_List$CELL_LINE[INDEX],
    COHORT = Cell_List$COHORT[INDEX],
    NUMBER_OF_CELLS_FOR_IMAGE = Cell_List$NUMBER_OF_CELLS[INDEX],
    UNIVERSAL_TRACK_ID = paste0(IMAGE, "...", CELL_NUMBER, "...", TRACK_ID),
    UNIVERSAL_SPOT_ID = paste0(UNIVERSAL_TRACK_ID, "...", NEW_ID)
  ) %>% 
  select(
    -c("INDEX")
  ) %>% 
  as.data.table()


  
  # Reorder columns to bring certain identifiers to the beginning
  col_name <- names(Trackmate_File)
  # Get the last 8 columns
  tail_name <- head(col_name, -8)
  
  new_order <- c(
    "IMAGE",
    "CELL_LINE_NUMBER",
    "CELL_NUMBER",
    "CELL_LINE",
    "COHORT",
    "UNIVERSAL_TRACK_ID",
    "UNIVERSAL_SPOT_ID",
    "NUMBER_OF_CELLS_FOR_IMAGE",
    tail_name
  )
  
  # Reorder the columns
  Trackmate_File <- Trackmate_File[, ..new_order]
  
  return(Trackmate_File)
}

# Define the number of cores to use
numCores <- detectCores() - 2  # Reserve two core for system processes
# Use mclapply to parallelize XMLtoCSV function over rows of Table_Final
Spot_Table <- mclapply(1:nrow(Table), Spot_Counter, mc.cores = numCores)
Spot_Table <- rbindlist(Spot_Table)

# Arrange Spot_Table by COHORT and UNIVERSAL_SPOT_ID
Spot_Table <- Spot_Table %>% 
  arrange(
    COHORT,
    UNIVERSAL_SPOT_ID
  ) %>% 
  filter(
    !is.na(NUMBER_OF_EVENTS)
  ) %>% 
  group_by(
    IMAGE
  ) %>% 
  mutate(
    REPLICATE = sub("_\\d+$", "", IMAGE)
  )
# Save the Spot_Table to a CSV file
Table_Path <- file.path(Table_Save_Folder, "TrackMate_Information.csv")
fwrite(Spot_Table, Table_Path)

# Cleanup the environment and run garbage collection
rm(list = ls())
gc()