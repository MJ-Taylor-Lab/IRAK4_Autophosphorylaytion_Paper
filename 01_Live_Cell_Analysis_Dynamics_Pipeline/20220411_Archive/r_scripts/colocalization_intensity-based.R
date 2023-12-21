#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

parameters_path = args[1]

# # CLEAN ENVIRONMENT----
# remove(list = ls())
# gc(reset = TRUE)
# pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
# 

# parameters_path = "/Users/u_deliz/Desktop/NewPipeline/Input/parameter_tables/"
# parameters_path = "/raven/u/deliz/new_pipeline/Input/parameter_tables"

new_image_ending = "_intensity_ref.tif"
results_table_name = "_intensity.csv.gz"


if("pacman" %in% rownames(installed.packages()) == FALSE)
{install.packages("pacman")}

pacman::p_load(dplyr, igraph, parallel, tidyr, data.table, ff, changepoint)
setDTthreads(parallel::detectCores(logical = F))

# Data frame must have x, y and id columns
ChangepointFx <- function(df){
  tryCatch({
    # Expand time list
    x_sequence <- min(df$x):max(df$x)
    # List of times to remove after running changepoint
    remove_x <- setdiff(x_sequence, df$x)
    remove_x <- !x_sequence %in%remove_x
    
    # Adjust df table
    df <- 
      df %>% 
      complete(
        # Expand times
        x = seq(min(x), max(x))
      ) %>% 
      fill(
        # Add missing y with previous value
        y,
        .direction = "down"
      ) %>% 
      mutate(
        # Log transform y
        y = log(y + 1)
      ) %>% 
      distinct()
    
    # Run changepoint
    m.pelt <- cpt.mean(df$y, method = "PELT", penalty = "Manual", pen.value = 0.3, minseglen = 3)
    # Get y
    y <- param.est(m.pelt)$mean
    # Transform y back to linear
    y <- exp(y)-1
    # Expand y to match time
    len <- seg.len(m.pelt)
    y <- rep(y, times = len)
    # Get ID
    id <- df$id[1]
    # Get times
    x <- df$x
    # Make new table
    new_df <- suppressMessages(bind_cols(x, y, id))
    names(new_df) <- c("x", "y", "id")
    # Keep only original times
    new_df <- new_df[remove_x, ]
    return(new_df)
  }, error = function(e) {print("Error with ChangepointFx")})}

# Get directories
directories_list = file.path(parameters_path, "directories.csv")
directories_list = fread(directories_list)
processing_path = directories_list$path[directories_list$contains == "processing"]
extraction_path = file.path(processing_path, "05_IntensityExtraction")
colocalized_path = file.path(processing_path, "06_Colocalization")
input_path = directories_list$path[directories_list$contains == "input"]
summary_path = file.path(input_path, "summary.csv")
# Get image list
file_list = fread(summary_path)
# Get images table
image_list = NULL
image_list$table = paste0(file_list$protein_relative_path, results_table_name)
image_list$table = file.path(extraction_path, image_list$table)
# Get list of files that need colocalization
colocalization_list <- as_tibble(image_list)
# Get colocalization list
colocalization_list <-
  colocalization_list %>% 
  mutate(
    cell = dirname(table),
    image = dirname(cell),
    cohort = dirname(image),
    cohort = basename(cohort),
    cohort = file.path(colocalized_path, cohort)
  ) %>% 
  filter(
    file.exists(table)
  ) %>% 
  group_by(
    cell
  ) %>% 
  mutate(
    n = n()
  )

# Create re-colocalized path if it doesn't exist
if(!file.exists(colocalized_path)){
  dir.create(colocalized_path)
}

# Get list of files that need colocalization
AllColocalizationNeeded <-
  colocalization_list %>% 
  filter(
    n > 1
  )

Images <- unique(AllColocalizationNeeded$image)
ColocalizeImage <- function(ImageX){
  tryCatch({
    ColocalizationNeeded <-
      AllColocalizationNeeded %>% 
      filter(
        image == Images[ImageX]
      )
    
    # Screen which frames have multiple proteins
    Cells <- unique(ColocalizationNeeded$cell)
    ScreeningFx <- function(CellX){
      tryCatch({
        # Get tables
        TableList <-
          ColocalizationNeeded %>% 
          filter(
            cell == Cells[CellX]
          )
        Tables <- lapply(TableList$table, fread)
        Tables <- rbindlist(Tables, fill = T)
        # Get only frames with multiple puncta
        Tables <-
          Tables %>% 
          arrange(
            PROTEIN,
            TRACK_ID,
            FRAME
          ) %>% 
          ungroup() %>% 
          group_by(
            FRAME
          ) %>% 
          mutate(
            N = NROW(unique(PROTEIN)),
          ) %>% 
          filter(
            N > 1
          )
        return(Tables)
      }, error = function(e){print(paste("ERROR with ScreeningFx CellX =", CellX))})
    }
    ScreenedTables <- mclapply(1:NROW(Cells), ScreeningFx, mc.cores = detectCores(logical = F))
    ScreenedTables <- ScreenedTables[(which(sapply(ScreenedTables,is.list), arr.ind=TRUE))]
    ScreenedTables <- rbindlist(ScreenedTables, fill=TRUE)
    
    SplitCellTables <-
      ScreenedTables %>% 
      mutate(
        CELL_PATH = dirname(RELATIVE_PATH),
        CELL_PATH = file.path(extraction_path, CELL_PATH),
        IMAGE_PATH = paste0(RELATIVE_PATH, new_image_ending),
        IMAGE_PATH = file.path(extraction_path, IMAGE_PATH)
      ) %>% 
      group_split(
        CELL_PATH
      )
    
    remove(ScreenedTables)
    
    GetImageTablePairs <- function(CellTable){
      tryCatch({
        # Get combinations
        Proteins = unique(CellTable$PROTEIN)
        Proteins = expand.grid(Proteins, Proteins)
        names(Proteins) = c("IMAGE", "TABLE")
        Proteins <-
          Proteins %>%
          filter(
            IMAGE != TABLE
          ) %>% 
          mutate(
            COMPLEMENTARY_PROTEIN = IMAGE,
            SOURCE_TABLE = TABLE,
            
            IMAGE = paste0(IMAGE, new_image_ending),
            IMAGE = file.path(CellTable$CELL_PATH[1], IMAGE),
            
            TABLE = paste0(TABLE, results_table_name),
            TABLE = file.path(CellTable$CELL_PATH[1], TABLE)
          ) %>% 
          filter(
            file.exists(IMAGE),
            file.exists(TABLE)
          ) %>% 
          group_by(
            SOURCE_TABLE
          ) %>% 
          mutate(
            COLUMN_NAME = 1:n(),
          )
        return(Proteins)
      }, error = function(e){print(paste("ERROR with GetImageTablePairs"))})
    }
    PairsList <- mclapply(SplitCellTables, GetImageTablePairs, mc.cores = detectCores(logical = F))
    PairsList <- PairsList[(which(sapply(PairsList,is.list), arr.ind=TRUE))]
    PairsList <- rbindlist(PairsList, fill = TRUE)
    
    GetIntensities <- function(PairX){
      tryCatch({
      print(paste("ImageX =", ImageX, "of", NROW(Images), "+ PairX =", PairX, "of", NROW(PairsList)))
      # Get parameters
      image_path <- PairsList$IMAGE[PairX]
      table_path <- PairsList$TABLE[PairX]
      CellTable <- fread(table_path)
      puncta_radius = CellTable$PUNCTA_DIAMETER[1]/2
      
      #Import images
      img = ijtiff::read_tif(image_path, frames = unique(CellTable$FRAME), msg = FALSE)
      GetFrameTable <- function(FrameX){
        tryCatch({
          img_frame = which(unique(CellTable$FRAME)==FrameX)
          Z = img[,,,img_frame]
          X = NROW(Z)
          X = rep(1:X, NCOL(Z))
          
          Y = NCOL(Z)
          Y = rep(1:Y, each = NROW(Z))
          
          Z = as.vector(Z)
          # Flipped on purpose. ijtiff flips coordinates
          frame_img = cbind(Y, X, Z)
          frame_img = as.data.frame(frame_img)
          names(frame_img) = c("x", "y", "z")
          frame_img$t = FrameX
          
          FrameCellTable <- CellTable[CellTable$FRAME == FrameX,]
          
          Table <- NULL
          Table[[1]] = frame_img
          Table[[2]] = as.data.frame(FrameCellTable)
          
          return(Table)
        }, error = function(e){print(paste("     ERROR with GetFrameTable FrameX =", FrameX))})
      }
      N_FRAMES = unique(CellTable$FRAME)
      N_FRAMES = sort(N_FRAMES)
      Tables <- mclapply(N_FRAMES, GetFrameTable, mc.cores = detectCores(logical = F))
      # Get frames
      FrameFx <- function(TableX){
        tryCatch({
          frame_img = TableX[[1]]
          FrameShortFullTable = TableX[[2]]
          # Get pixel intensity
          SubPixelLocalization <- function(SpotX){
            # Get coordinates
            x.coordinate = FrameShortFullTable$POSITION_X[SpotX]
            x.coordinate = as.numeric(x.coordinate)
            y.coordinate = FrameShortFullTable$POSITION_Y[SpotX]
            y.coordinate = as.numeric(y.coordinate)
            # Filter out spot
            TotalIntensity <- frame_img[frame_img$x >= x.coordinate - puncta_radius - .5 &
                                          frame_img$x <= x.coordinate + puncta_radius + .5 &
                                          frame_img$y >= y.coordinate - puncta_radius - .5 &
                                          frame_img$y <= y.coordinate + puncta_radius + .5,]
            
            TotalIntensity$x.low <- TotalIntensity$x - .5
            TotalIntensity$x.high <- TotalIntensity$x + .5
            TotalIntensity$x.coord.low <- x.coordinate - puncta_radius
            TotalIntensity$x.coord.high <- x.coordinate + puncta_radius
            TotalIntensity$x.high.capture <-
              data.table::fifelse(
                TotalIntensity$x.coord.high > TotalIntensity$x.high,
                TotalIntensity$x.high,
                TotalIntensity$x.coord.high
              )
            TotalIntensity$x.low.capture <-
              data.table::fifelse(
                x.coordinate - puncta_radius < TotalIntensity$x.low,
                TotalIntensity$x.low,
                TotalIntensity$x.coord.low
              )
            TotalIntensity$x.capture <- TotalIntensity$x.high.capture - TotalIntensity$x.low.capture
            
            TotalIntensity$y.low <- TotalIntensity$y - .5
            TotalIntensity$y.high <- TotalIntensity$y + .5
            TotalIntensity$y.coord.low <- y.coordinate - puncta_radius
            TotalIntensity$y.coord.high <- y.coordinate + puncta_radius
            TotalIntensity$y.high.capture <-
              data.table::fifelse(
                TotalIntensity$y.coord.high > TotalIntensity$y.high,
                TotalIntensity$y.high,
                TotalIntensity$y.coord.high
              )
            TotalIntensity$y.low.capture <-
              data.table::fifelse(
                TotalIntensity$y.coord.low < TotalIntensity$y.low,
                TotalIntensity$y.low,
                TotalIntensity$y.coord.low
              )
            
            TotalIntensity$y.capture <- TotalIntensity$y.high.capture - TotalIntensity$y.low.capture
            
            TotalIntensity$SPOT_AREA <- TotalIntensity$x.capture * TotalIntensity$y.capture
            TotalIntensity$TOTAL_INTENSITY <- TotalIntensity$SPOT_AREA*TotalIntensity$z
            
            TotalIntensity = TotalIntensity[TotalIntensity$x.capture > 0,]
            TotalIntensity = TotalIntensity[TotalIntensity$y.capture > 0,]
            
            Result <- NULL
            Result$UNIVERSAL_SPOT_ID = FrameShortFullTable$UNIVERSAL_SPOT_ID[SpotX]
            Result$UNIVERSAL_TRACK_ID = FrameShortFullTable$UNIVERSAL_TRACK_ID[SpotX]
            Result$FRAME = FrameShortFullTable$FRAME[SpotX]
            Result$TOTAL_INTENSITY = sum(TotalIntensity$TOTAL_INTENSITY)
            Result$STANDARD_DEVIATION = sd(TotalIntensity$TOTAL_INTENSITY)
            
            return(Result)
          }
          Intensities <- lapply(1:NROW(FrameShortFullTable), SubPixelLocalization)
          Intensities <- data.table::rbindlist(Intensities)
          return(Intensities)
        }, error = function(e){print(paste("     ERROR with FrameFx"))})
      }
      Intensities <- mclapply(Tables, FrameFx, mc.cores = detectCores(logical = F))
      Intensities <- Intensities[(which(sapply(Intensities,is.list), arr.ind=TRUE))]
      Intensities <- data.table::rbindlist(Intensities)
      # Add complementary protein name
      Intensities$COMPLEMENTARY_PROTEIN = PairsList$COMPLEMENTARY_PROTEIN[PairX]
      COLUMN_NUMBER <- PairsList$COLUMN_NAME[PairX]
      names(Intensities) <-c(
        "UNIVERSAL_SPOT_ID",
        "UNIVERSAL_TRACK_ID",
        "FRAME",
        paste0("COMPLEMENTARY_TOTAL_INTENSITY_", COLUMN_NUMBER),
        paste0("COMPLEMENTARY_STANDARD_DEVIATION_", COLUMN_NUMBER),
        paste0("COMPLEMENTARY_PROTEIN_", COLUMN_NUMBER)
      )
      
      SplitTracks <- Intensities[,2:4]
      names(SplitTracks) <- c(
        "id",
        "x",
        "y"
      )
      
      SplitTracks <- 
        SplitTracks %>%
        group_by(
          id
        ) %>% 
        mutate(
          LIFETIME = n()
        ) %>% 
        filter(
          LIFETIME >= 6
        ) %>% 
        arrange(
          LIFETIME,
          id,
          x
        ) %>% 
        ungroup() %>% 
        select(
          id,
          x,
          y
        ) %>% 
        group_split(
          id
        )
      
      # Run changepoint
      CPResults <- mclapply(SplitTracks, ChangepointFx, mc.cores = detectCores(logical = F))
      CPResults <- CPResults[(which(sapply(CPResults,is.list), arr.ind=TRUE))]
      CPResults <- rbindlist(CPResults)
      
      # Rename fields
      CPResults$id <- paste(CPResults$id, CPResults$x, sep = "...")
      CPResults$x <- NULL
      names(CPResults) <- c(paste0("COMPLEMENTARY_CHANGEPOINT_RAW_INTENSITY_", COLUMN_NUMBER), "UNIVERSAL_SPOT_ID")
      # Add changepoint to table
      Intensities <- merge(Intensities, CPResults, by = "UNIVERSAL_SPOT_ID", fill = T)
      
      TableName <- paste(PairsList$SOURCE_TABLE[PairX],  PairsList$COMPLEMENTARY_PROTEIN[PairX], "colocalization_intensity.csv.gz", sep = "_")
      TableName <- file.path(dirname(PairsList$IMAGE[PairX]), TableName)
      file.remove(TableName, showWarnings = F)
      fwrite(Intensities, TableName, row.names = F, na = "")
      
      return(TableName)
      }, error = function(e){print(paste("ERROR with GetIntensities. PairX =", PairX))})
    }
    
    lapply(1:NROW(PairsList), GetIntensities)
    
  }, error = function(e){print(paste("ERROR with ColocalizeImage ImageX =", ImageX))})
}
lapply(1:NROW(Images), ColocalizeImage)
print('colocalization-intensity based is now complete')
