#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

parameters_path = args[1]

# CLEAN ENVIRONMENT----
# remove(list = ls())
# gc(reset = TRUE)
# pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
# 
# parameters_path = "/Users/u_deliz/Desktop/NewPipeline/Input/parameter_tables/"

# Ending of files
intensity_ending = "_intensity.csv.gz"
colocalization_intensity_ending = "_colocalization_intensity.csv.gz"
colocalization_cost_filename = "colocalization_cost.csv.gz"
colocalization_group_filename = "colocalization_groups.csv.gz"
changepoint_ending = "_changepoint.csv.gz"


if("pacman" %in% rownames(installed.packages()) == FALSE)
{install.packages("pacman")}

pacman::p_load(dplyr, stringr, parallel, tidyr, data.table, ff)
setDTthreads(parallel::detectCores(logical = F))

# Get directories
directories_list = file.path(parameters_path, "directories.csv")
directories_list = fread(directories_list)
processing_path = directories_list$path[directories_list$contains == "processing"]
changepoint_path = file.path(processing_path, "07_Changepoint")
input_path = directories_list$path[directories_list$contains == "input"]
output_path = directories_list$path[directories_list$contains == "output"]
summary_path = file.path(input_path, "summary.csv")
# Get image list
file_list = fread(summary_path)
# Prep up
file_list <-
  file_list %>% 
  mutate(
    image = dirname(dirname(protein_relative_path)),
    image = basename(image),
    date = substr(image, 0, 8),
    date =  as.Date(date, format = '%Y%m%d'),
    exposure = ifelse(word(exposure, -1) == "ms", as.numeric(word(exposure, 1))/1000, ifelse(word(exposure, -1) == "s", as.numeric(word(exposure, 1)), exposure)),
    direction = direction*pi/180,
    angle = angle*pi/180
  )

# Separate calibrations and the rest
calibration_list = file_list[file_list$cohort=="Calibrations",]
Path <- calibration_list$protein_relative_path
Path <- dirname(Path)
Path <- file.path(changepoint_path, Path)
Path <- file.exists(Path)
calibration_list <- calibration_list[Path]

image_list = file_list[file_list$cohort!="Calibrations",]
Path <- image_list$protein_relative_path
Path <- dirname(Path)
Path <- file.path(changepoint_path, Path)
Path <- file.exists(Path)
image_list <- image_list[Path]

# Pairing
GetCalibrationImages <- function(ImageX){
  tryCatch({
    # Get calibration date
    IMAGE_DATE = image_list$date[ImageX]
    CALIBRATION_DATES = calibration_list$date
    DATE_DIFFERENCE = abs(CALIBRATION_DATES - IMAGE_DATE)
    CALIBRATION_DATE = CALIBRATION_DATES[which.min(DATE_DIFFERENCE)]
    
    # Calibration File
    CalibrationFile <-
      calibration_list %>% 
      filter(
        channel == image_list$channel[ImageX]
      )
    
    if(NROW(CalibrationFile)>0){
      # Get calibration
      # Nearest laser
      PowerFiltered = CalibrationFile[which.min(abs(CalibrationFile$power - image_list$power[ImageX])),]
      # Nearest exposure
      ExposureFiltered = PowerFiltered[which.min(abs(PowerFiltered$exposure - image_list$exposure[ImageX])),]
      # Nearest direction
      DirectionFiltered = ExposureFiltered[which.min(min(abs(ExposureFiltered$direction-image_list$direction[ImageX]), 360-abs(ExposureFiltered$direction-image_list$direction[ImageX]))),]
      # Nearest angle
      AngleFiltered = DirectionFiltered[which.min(min(abs(DirectionFiltered$angle-image_list$angle[ImageX]), 360-abs(DirectionFiltered$angle-image_list$angle[ImageX]))),]
      # Nearest date
      DateFiltered = AngleFiltered[which.min(abs(AngleFiltered$date - AngleFiltered$date[ImageX])),]
      # Select latest date
      CalibrationImage = DateFiltered$image[NROW(DateFiltered)]
      CalibrationProtein =  AngleFiltered$protein_name[1]
    } else{
      CalibrationImage = NA
      CalibrationProtein = NA
    }
    # Create export table
    ExportTable <- image_list[ImageX,]
    ExportTable$CALIBRATION_IMAGE = CalibrationImage
    ExportTable$FLUOROPHORE = CalibrationProtein
    
    return(ExportTable)
  }, error = function(e){print(paste("ERROR with MoveNoColocalizationNeededFx ImageX =", ImageX))})
}
PairedList <- mclapply(1:NROW(image_list), GetCalibrationImages, mc.cores = detectCores(logical = F))
PairedList <- rbindlist(PairedList, fill = TRUE, use.names = TRUE)
PairedList <- PairedList %>% distinct()

CalibrationImages <-
  PairedList %>% 
  select(
    CALIBRATION_IMAGE,
    FLUOROPHORE
  ) %>%
  distinct() %>% 
  drop_na()

GetCalibrationIntensity <- function(ImageX){
  tryCatch({
    CALIBRATION_IMAGE = CalibrationImages$CALIBRATION_IMAGE[ImageX]
    FLUOROPHORE = CalibrationImages$FLUOROPHORE[ImageX]
    CALIBRATION_IMAGE = file.path(changepoint_path, "Calibrations", CALIBRATION_IMAGE, "Cell_1", paste0(FLUOROPHORE, "_intensity.csv.gz"))
    
    IntensitiesTable <- fread(CALIBRATION_IMAGE)
    
    IntensitiesTable <-
      IntensitiesTable %>% 
      mutate(
        CALIBRATION_IMAGE = IMAGE
      ) %>% 
      filter(
        SPOT_AREA == PUNCTA_DIAMETER^2
      ) %>% 
      group_by(
        CALIBRATION_IMAGE
      ) %>% 
      summarize(
        CALIBRATION_STANDARD_DEVIATION = sd(TOTAL_INTENSITY),
        CALIBRATION_TOTAL_INTENSITY = median(TOTAL_INTENSITY)
      )
    return(IntensitiesTable)
  }, error = function(e){print(paste("ERROR with GetCalibrationIntensity ImageX =", ImageX))})
}
CalibrationIntensities <- mclapply(1:NROW(CalibrationImages), GetCalibrationIntensity, mc.cores = detectCores(logical = F))
CalibrationIntensities <- CalibrationIntensities[(which(sapply(CalibrationIntensities,is.list), arr.ind=TRUE))]
CalibrationIntensities <- rbindlist(CalibrationIntensities, fill = TRUE, use.names = TRUE)
# Pair image and calibration
PairedList <- as.data.frame(PairedList)
CalibrationIntensities <- as.data.frame(CalibrationIntensities)
PairedCalibrations <- merge(PairedList, CalibrationIntensities, by = "CALIBRATION_IMAGE", all = TRUE)
names(PairedCalibrations) <- toupper(names(PairedCalibrations))

PairedCalibrations <-
  PairedCalibrations %>%
  mutate(
    CELL_PATH = dirname(PROTEIN_RELATIVE_PATH),
    CALIBRATION_TOTAL_INTENSITY = ifelse(is.na(CALIBRATION_TOTAL_INTENSITY), 1, CALIBRATION_TOTAL_INTENSITY),
    FLUOROPHORE = ifelse(is.na(FLUOROPHORE), CHANNEL, FLUOROPHORE)
  ) %>%
  distinct()

Cells <- unique(PairedCalibrations$CELL_PATH)
CombineCellTables <- function(CellX){
  tryCatch({
    print(paste("CombineCellTables - CellX =", CellX))
    # Get tables
    CellPairedCalibrations <-
      PairedCalibrations %>% filter(
        CELL_PATH ==  Cells[CellX]
      ) %>% 
      mutate(
        PROTEIN = PROTEIN_NAME 
      ) %>% 
      select(
        PROTEIN,
        CALIBRATION_IMAGE,
        CALIBRATION_TOTAL_INTENSITY,
        CALIBRATION_STANDARD_DEVIATION
      )
    
    # Get cell path
    CellPath <- file.path(changepoint_path, Cells[CellX])
    
    # Get protein combinations
    Proteins = unique(CellPairedCalibrations$PROTEIN)
    Proteins = expand.grid(Proteins, Proteins)
    names(Proteins) = c("IMAGE", "TABLE")
    Proteins <-
      Proteins %>%
      filter(
        IMAGE != TABLE
      ) %>% 
      mutate(
        TABLE = paste0(IMAGE, "_", TABLE, colocalization_intensity_ending)
      )
    
    # Get tables
    IntensityTables <- paste0(CellPairedCalibrations$PROTEIN, intensity_ending)
    IntensityTables <- file.path(CellPath, IntensityTables)
    IntensityTables <- IntensityTables[file.exists(IntensityTables)]
    if(NROW(IntensityTables)>0){
      
      IntensityTables <- lapply(IntensityTables, fread)
      IntensityTables <- rbindlist(IntensityTables, fill = TRUE, use.names = TRUE)
      # Pair calibration
      IntensityTables <- as.data.frame(IntensityTables)
      CellPairedCalibrations <- as.data.frame(CellPairedCalibrations)
      IntensityTables <- merge(IntensityTables, CellPairedCalibrations, by = "PROTEIN", all = TRUE)
      
      # Add colocalization intensities (image and coordinates flipped)
      ColocalizationIntensityTables <- file.path(CellPath, Proteins$TABLE)
      ColocalizationIntensityTables <- ColocalizationIntensityTables[file.exists(ColocalizationIntensityTables)]
      if(NROW(ColocalizationIntensityTables)>0){
        ColocalizationIntensityTables <- lapply(ColocalizationIntensityTables, fread)
        ColocalizationIntensityTables <- rbindlist(ColocalizationIntensityTables, fill = TRUE, use.names = TRUE)
        
        # Normalize colocalization intensity
        ComplementaryProteins <- names(ColocalizationIntensityTables)
        ComplementaryProteinsIndex <- substr(ComplementaryProteins, 0, nchar("COMPLEMENTARY_PROTEIN_"))
        ComplementaryProteinsIndex <- which(ComplementaryProteinsIndex == "COMPLEMENTARY_PROTEIN_")
        ComplementaryProteins <- ComplementaryProteins[ComplementaryProteinsIndex]
        ComplementaryProteins <- grep("COMPLEMENTARY_PROTEIN_", ComplementaryProteins)
        
        # Get names of columns to fill
        FillNames <- c("COMPLEMENTARY_TOTAL_INTENSITY_",
                       "COMPLEMENTARY_STANDARD_DEVIATION_",
                       "COMPLEMENTARY_PROTEIN_",
                       "COMPLEMENTARY_CHANGEPOINT_RAW_INTENSITY_")
        FillNames <- expand.grid(FillNames, ComplementaryProteins)
        FillNames <- paste0(FillNames$Var1, FillNames$Var2)
        # Fill data
        ColocalizationIntensityTables <-
          ColocalizationIntensityTables %>% 
          group_by(
            UNIVERSAL_SPOT_ID
          ) %>% 
          fill(
            all_of(FillNames),
            .direction = "updown"
          ) %>% 
          distinct()
        
        ColocalizationIntensityNormalization <- function(ProteinX){
          tryCatch({
            # Get protein name
            ProteinColumnName <- paste0("COMPLEMENTARY_PROTEIN_", ProteinX)
            ProteinName = which(names(ColocalizationIntensityTables) == ProteinColumnName)
            ProteinName <- ColocalizationIntensityTables[, ProteinName]
            # Get calibration intensity
            ProteinName <- as.data.frame(ProteinName)
            CellPairedCalibrations <- as.data.frame(CellPairedCalibrations)
            ProteinCalibration <- merge(ProteinName, CellPairedCalibrations, by.x = ProteinColumnName, by.y = "PROTEIN", all = T, fill = T)
            
            # Get intensity
            IntensityColumnName <- paste0("COMPLEMENTARY_TOTAL_INTENSITY_", ProteinX)
            ProteinIntensity = which(names(ColocalizationIntensityTables) == IntensityColumnName)
            ProteinIntensity <- ColocalizationIntensityTables[, ProteinIntensity]
            ProteinIntensity <- ProteinIntensity/ProteinCalibration$CALIBRATION_TOTAL_INTENSITY
            names(ProteinIntensity) <- paste0("COMPLEMENTARY_NORMALIZED_INTENSITY_", ProteinX)
            RegProteinIntensity <- ProteinIntensity
            
            IntensityColumnName <- paste0("COMPLEMENTARY_CHANGEPOINT_RAW_INTENSITY_", ProteinX)
            ProteinIntensity = which(names(ColocalizationIntensityTables) == IntensityColumnName)
            ProteinIntensity <- ColocalizationIntensityTables[, ProteinIntensity]
            ProteinIntensity <- ProteinIntensity/ProteinCalibration$CALIBRATION_TOTAL_INTENSITY
            names(ProteinIntensity) <- paste0("COMPLEMENTARY_CHANGEPOINT_NORMALIZED_INTENSITY_", ProteinX)
            
            ProteinIntensity <- bind_cols(ProteinIntensity, RegProteinIntensity)
            
            ProteinIntensity$UNIVERSAL_SPOT_ID <- ColocalizationIntensityTables$UNIVERSAL_SPOT_ID
            return(ProteinIntensity)
          }, error = function(e){print(paste("ERROR with ColocalizationIntensityNormalization ProteinX =", ProteinX))})
        }
        ColocalizationNormalizedIntensities <- lapply(ComplementaryProteins, ColocalizationIntensityNormalization)
        ColocalizationNormalizedIntensities <- ColocalizationNormalizedIntensities[(which(sapply(ColocalizationNormalizedIntensities,is.list), arr.ind=TRUE))]
        ColocalizationNormalizedIntensities <- do.call(cbind, ColocalizationNormalizedIntensities)
        ColocalizationNormalizedIntensities <- ColocalizationNormalizedIntensities[!duplicated(as.list(ColocalizationNormalizedIntensities))]
        ColocalizationIntensityTables <- as.data.frame(ColocalizationIntensityTables)
        ColocalizationNormalizedIntensities <- as.data.frame(ColocalizationNormalizedIntensities)
        ColocalizationNormalizedIntensities <- merge(ColocalizationIntensityTables, ColocalizationNormalizedIntensities, by = "UNIVERSAL_SPOT_ID", all = TRUE)
        
        ColocalizationNormalizedIntensities$FRAME <- NULL
        ColocalizationNormalizedIntensities$UNIVERSAL_TRACK_ID <- NULL
        
        IntensityTables <- as.data.frame(IntensityTables)
        ColocalizationNormalizedIntensities <- as.data.frame(ColocalizationNormalizedIntensities)
        IntensityTables <- merge(IntensityTables, ColocalizationNormalizedIntensities, by = "UNIVERSAL_SPOT_ID", all = TRUE)
        remove(ColocalizationNormalizedIntensities, ColocalizationIntensityTables)
      }
      
      # Add colocalization data (coordinates paired, cost)
      ColocalizationCost <- file.path(CellPath, colocalization_cost_filename)
      ColocalizationCost <- ColocalizationCost[file.exists(ColocalizationCost)]
      ColocalizationCost <- if(NROW(ColocalizationCost)>0){
        ColocalizationCost <- fread(ColocalizationCost)
        IntensityTables <- as.data.frame(IntensityTables)
        ColocalizationCost <- as.data.frame(ColocalizationCost)
        IntensityTables <- merge(IntensityTables, ColocalizationCost, by = "UNIVERSAL_SPOT_ID", all = TRUE)
        remove(ColocalizationCost)
      }
      # Add colocalization data (coordinates paired, edge)
      ColocalizationGroups <- file.path(CellPath, colocalization_group_filename)
      ColocalizationGroups <- ColocalizationGroups[file.exists(ColocalizationGroups)]
      ColocalizationGroups <- if(NROW(ColocalizationGroups)>0){
        ColocalizationGroups <- fread(ColocalizationGroups)
        IntensityTables <- as.data.frame(IntensityTables)
        ColocalizationGroups <- as.data.frame(ColocalizationGroups)
        IntensityTables <- merge(IntensityTables, ColocalizationGroups, by = "UNIVERSAL_TRACK_ID", all = TRUE)
        remove(ColocalizationGroups)
      }
      
      # Add changepoint analysis
      ChangepointTables <- paste0(CellPairedCalibrations$PROTEIN, changepoint_ending)
      ChangepointTables <- file.path(CellPath, ChangepointTables)
      ChangepointTables <- ChangepointTables[file.exists(ChangepointTables)]
      ChangepointTables <- if(NROW(ChangepointTables)>0){
        ChangepointTables <- lapply(ChangepointTables, fread)
        ChangepointTables <- ChangepointTables[(which(sapply(ChangepointTables,is.list), arr.ind=TRUE))]
        ChangepointTables <- rbindlist(ChangepointTables, fill = TRUE, use.names = TRUE)
        IntensityTables <- as.data.frame(IntensityTables)
        IntensityTables <- as.data.frame(IntensityTables)
        IntensityTables <- merge(IntensityTables, ChangepointTables, by = "UNIVERSAL_SPOT_ID", all = TRUE)
        IntensityTables <- as.data.table(IntensityTables)
        IntensityTables[, NORMALIZED_CHANGEPOINT_INTENSITY := CHANGEPOINT_TOTAL_INTENSITY/CALIBRATION_TOTAL_INTENSITY, by = PROTEIN]
        remove(ChangepointTables)
      }
      
      # Determine nearest neighbor to spot of the same protein
      NearestNeighborTable <-
        IntensityTables %>% 
        filter(
          !is.na(POSITION_X)
        ) %>% 
        group_by(
          PROTEIN,
          FRAME
        ) %>% 
        mutate(
          N = n()
        ) %>% 
        filter(
          N > 1
        ) %>% 
        select(
          PROTEIN,
          FRAME,
          PUNCTA_DIAMETER,
          UNIVERSAL_SPOT_ID,
          POSITION_X,
          POSITION_Y
        ) %>% 
        ungroup() %>% 
        group_split(
          PROTEIN,
          FRAME
        )
      # Nearest neighbor971
      RADIUS <- max(IntensityTables$PUNCTA_DIAMETER, na.rm = T) * 3
      NearestNeighborSearch <- function(ProteinFrameX){
        # Get nearest spot
        CoordiantesTable <- ProteinFrameX %>% select(POSITION_X, POSITION_Y)
        DistanceTable <- RANN::nn2(CoordiantesTable, k = 2)
        DistanceTable <- DistanceTable$nn.dists[,2]
        # Get number of spots nearby
        ClusterTable <- RANN::nn2(CoordiantesTable, searchtype = c("radius"), radius = RADIUS, k = NROW(ProteinFrameX))
        ClusterTable <- ClusterTable$nn.dists
        ClusterTable <- ClusterTable > 0 & ClusterTable <  1e+153 
        ClusterTable <- rowSums(ClusterTable)
        # Put results together
        DistanceResults <- ProteinFrameX %>% select(UNIVERSAL_SPOT_ID)
        DistanceResults$NEAREST_SPOT <- DistanceTable
        DistanceResults$SPOTS_WITHIN_RADIUS = ClusterTable
        DistanceResults$SPOT_RADIUS_LIMIT = RADIUS
        return(DistanceResults)
      }
      NeighborResults <- lapply(NearestNeighborTable, NearestNeighborSearch)
      NeighborResults <- NeighborResults[(which(sapply(NeighborResults,is.list), arr.ind=TRUE))]
      NeighborResults <- rbindlist(NeighborResults, fill = TRUE, use.names = TRUE)
      IntensityTables <- as.data.frame(IntensityTables)
      NeighborResults <- as.data.frame(NeighborResults)
      IntensityTables <- merge(IntensityTables, NeighborResults, by = "UNIVERSAL_SPOT_ID", all = TRUE)
      # Add missing values
      MissingIndex <- which(is.na(IntensityTables$NEAREST_SPOT))
      IntensityTables$NEAREST_SPOT[MissingIndex] <- Inf
      MissingIndex <- which(is.na(IntensityTables$SPOTS_WITHIN_RADIUS))
      IntensityTables$SPOTS_WITHIN_RADIUS[MissingIndex] <- 0
      MissingIndex <- which(is.na(IntensityTables$SPOT_RADIUS_LIMIT))
      IntensityTables$SPOT_RADIUS_LIMIT[MissingIndex] <- RADIUS
      
      # Get landing frame
      LANDING_FRAME = min(IntensityTables$FRAME)
      
      # Add actual time
      TimeTable <- file.path(dirname(CellPath), "timesteps.csv")
      if(file.exists(TimeTable)){
        TimeTable <- fread(TimeTable, header = FALSE)
        names(TimeTable) <- "TIME"
        TimeTable <- TimeTable %>% mutate(FRAME = 1:n(), TIME = TIME/1000)
        IntensityTables <- as.data.frame(IntensityTables)
        TimeTable <- as.data.frame(TimeTable)
        IntensityTables <- merge(IntensityTables, TimeTable, by = "FRAME")
      } else{
        IntensityTables <-
          IntensityTables %>% 
          mutate(
            TIME = FRAME*FRAME_RATE
          )
      }
      # Calculate basic parameters and save
      IntensityTables <- as.data.table(IntensityTables)
      IntensityTables[, MIN_X := min(POSITION_X/CALIBRATION_UM), by = UNIVERSAL_TRACK_ID]
      IntensityTables[, MIN_X := MIN_X > PUNCTA_DIAMETER*2]
      IntensityTables <- IntensityTables[IntensityTables$MIN_X]
      IntensityTables$MIN_X <- NULL
      
      IntensityTables[, MAX_X := max(ABSOLUTE_POSITION_X/CALIBRATION_UM), by = UNIVERSAL_TRACK_ID]
      IntensityTables[, MAX_X := MAX_X < WIDTH - PUNCTA_DIAMETER*2]
      IntensityTables <- IntensityTables[IntensityTables$MAX_X]
      IntensityTables$MAX_X <- NULL
      
      IntensityTables[, MIN_Y := min(POSITION_Y/CALIBRATION_UM), by = UNIVERSAL_TRACK_ID]
      IntensityTables[, MIN_Y := MIN_Y > PUNCTA_DIAMETER*2]
      IntensityTables <- IntensityTables[IntensityTables$MIN_Y]
      IntensityTables$MIN_Y <- NULL
      
      IntensityTables[, MAX_Y := max(ABSOLUTE_POSITION_Y/CALIBRATION_UM), by = UNIVERSAL_TRACK_ID]
      IntensityTables[, MAX_Y := MAX_Y < HEIGHT - PUNCTA_DIAMETER*2]
      IntensityTables <- IntensityTables[IntensityTables$MAX_Y]
      IntensityTables$MAX_Y <- NULL
      
      # Remove small spots
      IntensityTables <- IntensityTables[IntensityTables$SPOT_AREA >= (IntensityTables$PUNCTA_DIAMETER^2)]
      
      # Normalize spots
      IntensityTables[, NORMALIZED_INTENSITY := TOTAL_INTENSITY/CALIBRATION_TOTAL_INTENSITY, by = PROTEIN]
      
      # If infinite distance, put NA
      is.infinite(IntensityTables$NEAREST_SPOT)
      
      IntensityTables <- data.table(IntensityTables)
      invisible(lapply(names(IntensityTables),function(.name) set(IntensityTables, which(is.infinite(IntensityTables[[.name]])), j = .name,value =NA)))
      IntensityTables <- data.frame(IntensityTables)
      
      IntensityTables <-
        IntensityTables %>% 
        group_by(
          UNIVERSAL_TRACK_ID
        ) %>% 
        mutate(
          # Ligand density every 0.5 steps using log base 10
          LIGAND_DENSITY_CAT =
            # Round to nearest half log
            round(
              # Classifies ligands based on log using base 3.162 (10^.5)
              log(
                LIGAND_DENSITY,
                base = (10^.5)),
              digits = 0
            ) * 0.5,
          # Convert back to linear
          LIGAND_DENSITY_CAT = signif(10^LIGAND_DENSITY_CAT, 2),
          # Get frame number
          FRAMES_ADJUSTED = FRAME - min(FRAME),
          # Time adjusted
          TIME_ADJUSTED = TIME - min(TIME),
          # Get puncta lifetime
          LIFETIME = max(TIME) - min(TIME),
          # Time since first spot
          TIME_SINCE_LANDING = (FRAME - LANDING_FRAME)/FRAME_RATE,
          # Max total intensity of track
          MAX_NORMALIZED_INTENSITY = max(NORMALIZED_INTENSITY, na.rm = TRUE),
          # To be used later for overall delta and for categorizing de-novo and disassembling
          STARTING_NORMALIZED_INTENSITY = NORMALIZED_INTENSITY[1],
          # Ending intensity
          ENDING_NORMALIZED_INTENSITY = NORMALIZED_INTENSITY[n()],
          # Overall change in intensity from start to max
          START_TO_MAX_INTENSITY = MAX_NORMALIZED_INTENSITY - STARTING_NORMALIZED_INTENSITY,
          # For pointing out which frame contains the max intensity
          MAX_INTENSITY_TIME = ifelse(MAX_NORMALIZED_INTENSITY == NORMALIZED_INTENSITY, TIME, NA),
          # Get first frame to reach track max intensity in case of duplicates
          MAX_INTENSITY_TIME = min(MAX_INTENSITY_TIME, na.rm = TRUE)
        )
      
      CellData <-
        IntensityTables %>% 
        ungroup() %>% 
        select(
          one_of(
            "RELATIVE_PATH",
            
            "LIGAND",
            "LIGAND_DENSITY",
            
            "CHANNEL",
            "POWER",
            "EXCITATION",
            "EMMISION",
            "EXPOSURE",
            "ANGLE",
            "DIRECTION",
            "FOCUS",
            "OBJECTIVE",
            
            "WIDTH",
            "HEIGHT",
            "CALIBRATION_UM",
            
            "CELL_DIAMETER",
            "PUNCTA_DIAMETER",
            
            "SEGMENT_WITH",
            
            "TRACKMATE_THRESHOLD",
            "TRACKMATE_FRAME_GAP",
            "TRACKMATE_GAP_LINK_DISTANCE",
            "TRACKMATE_MAX_LINK_DISTANCE",
            
            "DISTANCE_THRESHOLD",
            "SPOT_RADIUS_LIMIT",
            "ASSOCIATION_TIME_THRESHOLD",
            
            "CELL_POSITION_X",
            "CELL_POSITION_Y",
            
            "TIME_START",
            "FRAME_RATE",
            
            "CALIBRATION_IMAGE",
            "CALIBRATION_TOTAL_INTENSITY",
            "CALIBRATION_STANDARD_DEVIATION"
          )
        ) %>% 
        distinct() %>% 
        fill(
          one_of(
            "DISTANCE_THRESHOLD",
            "ASSOCIATION_TIME_THRESHOLD"
          ),
          .direction = "updown"
        ) %>% 
        distinct()
      
      # Write file
      DestinationPath <- file.path(CellPath, "Parameters.csv.gz")
      file.remove(DestinationPath, showWarnings = F)
      fwrite(CellData, DestinationPath, row.names = F, na = "")
      remove(CellData)
      
      Parameterless <-
        IntensityTables %>% 
        ungroup() %>% 
        select(
          !one_of(
            "LIGAND",
            "LIGAND_DENSITY",
            
            "CHANNEL",
            "POWER",
            "EXCITATION",
            "EMMISION",
            "EXPOSURE",
            "ANGLE",
            "DIRECTION",
            "FOCUS",
            "OBJECTIVE",
            
            "WIDTH",
            "HEIGHT",
            "CALIBRATION_UM",
            
            "CELL_DIAMETER",
            "PUNCTA_DIAMETER",
            
            "SEGMENT_WITH",
            
            "TRACKMATE_THRESHOLD",
            "TRACKMATE_FRAME_GAP",
            "TRACKMATE_GAP_LINK_DISTANCE",
            "TRACKMATE_MAX_LINK_DISTANCE",
            
            "DISTANCE_THRESHOLD",
            "SPOT_RADIUS_LIMIT",
            "ASSOCIATION_TIME_THRESHOLD",
            
            "CELL_POSITION_X",
            "CELL_POSITION_Y",
            
            "TIME_START",
            "FRAME_RATE",
            
            "CALIBRATION_IMAGE",
            "CALIBRATION_TOTAL_INTENSITY",
            "CALIBRATION_STANDARD_DEVIATION"
          )
        )
      
      # Write file
      DestinationPath <- file.path(CellPath, "Analysis.csv.gz")
      file.remove(DestinationPath, showWarnings = F)
      fwrite(Parameterless, DestinationPath, row.names = F, na = "")
      remove(Parameterless)
      
      
      
      # Get names of colocalization
      if(exists("ComplementaryProteins")==FALSE){
        ColocalizationNames = "ColocalizationNames"
      } else{
        ColocalizationNames <- c(
          "COMPLEMENTARY_NORMALIZED_INTENSITY_","COMPLEMENTARY_CHANGEPOINT_NORMALIZED_INTENSITY_", "COMPLEMENTARY_PROTEIN_",
          "COLOCALIZATION_PROTEIN_", "COLOCALIZATION_SPOT_", "COLOCALIZED_DISTANCE_", "COLOCALIZATION_INCLUDED_")
        ColocalizationNames <- expand.grid(ColocalizationNames, ComplementaryProteins)
        ColocalizationNames <- paste0(ColocalizationNames$Var1, ColocalizationNames$Var2)
      }

      Essential <-
        IntensityTables %>% 
        ungroup() %>% 
        select(
          one_of(
            c(
              "RELATIVE_PATH",
              "PROTEIN",
              "COHORT",
              "LIGAND_DENSITY_CAT",
              "IMAGE",
              "CELL",
              "UNIVERSAL_TRACK_ID",
              "UNIVERSAL_SPOT_ID",
              
              "FRAME",
              "TIME",
              "FRAMES_ADJUSTED",
              "TIME_ADJUSTED",
              "TIME_SINCE_LANDING",
              "LIFETIME",
              
              "NORMALIZED_INTENSITY",
              "NORMALIZED_CHANGEPOINT_INTENSITY",
              "STARTING_NORMALIZED_INTENSITY",
              "MAX_NORMALIZED_INTENSITY",
              "START_TO_MAX_INTENSITY",
              "MAX_INTENSITY_TIME",
              
              "CELL_AREA",
              "ABSOLUTE_POSITION_X",
              "ABSOLUTE_POSITION_Y",
              
              "NEAREST_SPOT",
              "SPOTS_WITHIN_RADIUS",
              "COLOCALIZATION_GROUP",
              
              ColocalizationNames
            )
          )
        )
      
      # Write file
      DestinationPath <- file.path(CellPath, "Essential.csv.gz")
      file.remove(DestinationPath, showWarnings = F)
      fwrite(Essential, DestinationPath, row.names = F, na = "")
      
      return(CellPath)
    }
  }, error = function(e){print(paste("ERROR with CombineCellTables CellX =", CellX))})
}
CellAnalysis <- mclapply(1:NROW(Cells), CombineCellTables, mc.cores = detectCores(logical = F))
CellAnalysis <- unlist(CellAnalysis)
CellAnalysis <- CellAnalysis[file.exists(CellAnalysis)]

# Combine all cell tables
Images <- dirname(CellAnalysis)
Images <- unique(Images)

CombineImageTables <- function(ImageX){
  tryCatch({
    print(paste("ImageX =", ImageX, "of", nImages))
    # Get cell table paths
    CellsList <- CellAnalysis[dirname(CellAnalysis) == Images[ImageX]]
    Outputs <- c("Parameters.csv.gz", "Analysis.csv.gz", "Essential.csv.gz")
    # Function to combine tables
    TableByOutput <- function(OutputX){
      print(paste("Combining",basename(Images[ImageX]), " - ", OutputX))
      # Get path
      Path <- file.path(CellsList, OutputX)
      Path <- Path[file.exists(Path)]
      # Pull tables
      Table <- lapply(Path, fread)
      Table <- rbindlist(Table, fill = TRUE, use.names = TRUE)
      # Save
      DestinationPath <- file.path(Images[ImageX], OutputX)
      file.remove(DestinationPath, showWarnings = F)
      fwrite(Table, DestinationPath, row.names = F, na = "")
      if(OutputX =="Essential.csv.gz"){
        return(Table)
      }
    }
    CellsList <- lapply(Outputs, TableByOutput)
    CellsList <- rbindlist(CellsList)
    # Make image summary for future analysis
    CellSummary <-
      CellsList %>%
      filter(
        FRAMES_ADJUSTED == 0
      ) %>%
      group_by(
        LIGAND_DENSITY_CAT,
        COHORT,
        IMAGE,
        CELL,
        PROTEIN,
        CELL_AREA
      ) %>%
      summarize(
        SPOTS = n(),
        LIFETIME = mean(LIFETIME, na.rm = T),
        STARTING_NORMALIZED_INTENSITY = mean(STARTING_NORMALIZED_INTENSITY, na.rm = T),
        MAX_NORMALIZED_INTENSITY = mean(MAX_NORMALIZED_INTENSITY, na.rm = T),
        START_TO_MAX_INTENSITY = mean(START_TO_MAX_INTENSITY, na.rm = T)
      )
    
    DestinationPath <- file.path(Images[ImageX], "CellSummary.csv.gz")
    file.remove(DestinationPath, showWarnings = F)
    fwrite(CellSummary, DestinationPath, row.names = F, na = "")
    
    ProteinSummary <-
      CellSummary %>%
      group_by(
        LIGAND_DENSITY_CAT,
        COHORT,
        IMAGE,
        PROTEIN
      ) %>%
      summarize(
        CELLS = n(),
        SPOTS = sum(SPOTS),
        LIFETIME = mean(LIFETIME, na.rm = T),
        STARTING_NORMALIZED_INTENSITY = mean(STARTING_NORMALIZED_INTENSITY, na.rm = T),
        MAX_NORMALIZED_INTENSITY = mean(MAX_NORMALIZED_INTENSITY, na.rm = T),
        START_TO_MAX_INTENSITY = mean(START_TO_MAX_INTENSITY, na.rm = T)
      )
    
    DestinationPath <- file.path(Images[ImageX], "ProteinSummary.csv.gz")
    file.remove(DestinationPath, showWarnings = F)
    fwrite(ProteinSummary, DestinationPath, row.names = F, na = "")
    return(dirname(DestinationPath))
  }, error = function(e){print(paste("ERROR with CombineImageTables ImageX =", ImageX))})
}
nImages <- NROW(Images)
ProcessedImages <- lapply(1:nImages, CombineImageTables)

# Move processed images
for(Image in ProcessedImages){
  old_path = Image
  Cohort = basename(dirname(Image))
  dir.create(file.path(output_path, Cohort))
  Image = basename(Image)
  new_path = file.path(output_path, Cohort, Image)
  file.move(old_path, new_path)
}

# Move calibrations
for(CalibrationX in 1:NROW(calibration_list)){
  Image = calibration_list$image[CalibrationX]
  Cohort = "Calibrations"
  old_path = file.path(changepoint_path, Cohort, Image)
  new_path = file.path(output_path, Cohort, Image)
  dir.create(file.path(output_path, Cohort))
  file.move(old_path, new_path)
}
