#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

parameters_path = args[1]

# # CLEAN ENVIRONMENT----
# remove(list = ls())
# gc(reset = TRUE)
# pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
#
# 
# parameters_path = "/Users/u_deliz/Desktop/NewPipeline/Input/parameter_tables/"
# parameters_path = '/raven/u/deliz/new_pipeline/Input/parameter_tables'
new_image_ending = "_intensity_ref.tif"
results_table_name = "_intensity.csv.gz"

if("pacman" %in% rownames(installed.packages()) == FALSE)
{install.packages("pacman")}

pacman::p_load(dplyr, igraph, parallel, tidyr, data.table, ff, RANN)
setDTthreads(parallel::detectCores(logical = F))

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
# Get cohort list
Cohorts <- unique(colocalization_list$cohort)
# Create cohort path if it doesn't exist
for(Cohort in Cohorts){
  if(!file.exists(Cohort)){
    dir.create(Cohort)
  }
}

# Get list of no colocalization needed
NoColocalizationNeeded <-
  colocalization_list %>% 
  filter(
    n == 1
  )

# Get list of files that need colocalization
ColocalizationNeeded <-
  colocalization_list %>% 
  filter(
    n > 1
  )

# Get list of cells
Cells <- unique(ColocalizationNeeded$cell)
# Compute cost (i.e., find nearest neighbors for each puncta)
CellColocalizationFx <- function(CellX){
  tryCatch({
    print(paste("CellColocalizationFx - CellX =", CellX))
    # Import table list
    CellTable <-
      ColocalizationNeeded %>% 
      filter(
        cell == Cells[CellX]
      )
    
    # Read tables
    CellTable <- lapply(CellTable$table, fread)
    CellTable <- rbindlist(CellTable, fill = T)
    
    # Get parameters
    Radius <- max(CellTable$PUNCTA_DIAMETER)
    # Get minimum association time
    AssocTime <- max(CellTable$TRACKMATE_FRAME_GAP, na.rm = T) + 1
    
    # Filter out edge puncta
    CellTable[, MIN_X := min(POSITION_X/CALIBRATION_UM), by = UNIVERSAL_TRACK_ID]
    CellTable[, MIN_X := (MIN_X > PUNCTA_DIAMETER*2)]
    CellTable <- CellTable[CellTable$MIN_X]
    CellTable$MIN_X <- NULL
    
    CellTable[, MAX_X := max(ABSOLUTE_POSITION_X/CALIBRATION_UM), by = UNIVERSAL_TRACK_ID]
    CellTable[, MAX_X := MAX_X < (WIDTH - PUNCTA_DIAMETER*2)]
    CellTable <- CellTable[CellTable$MAX_X]
    CellTable$MAX_X <- NULL
    
    CellTable[, MIN_Y := min(POSITION_Y/CALIBRATION_UM), by = UNIVERSAL_TRACK_ID]
    CellTable[, MIN_Y := MIN_Y > (PUNCTA_DIAMETER*2)]
    CellTable <- CellTable[CellTable$MIN_Y]
    CellTable$MIN_Y <- NULL
    
    CellTable[, MAX_Y := max(ABSOLUTE_POSITION_Y/CALIBRATION_UM), by = UNIVERSAL_TRACK_ID]
    CellTable[, MAX_Y := MAX_Y < (HEIGHT - PUNCTA_DIAMETER*2)]
    CellTable <- CellTable[CellTable$MAX_Y]
    CellTable$MAX_Y <- NULL
    
    # Filter out small puncta
    CellTable <- CellTable[CellTable$SPOT_AREA >= (CellTable$PUNCTA_DIAMETER^2)]

    # Get only frames with multiple puncta
    CellTable <-
      CellTable %>%
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
      ) %>% 
      arrange(
        N
      ) %>% 
      group_split(
        IMAGE,
        CELL,
        FRAME
      ) 
    
    if(NROW(CellTable) != 0){
      
      # Perform cost (distance) function
      CostFx <- function(FrameTable){
        tryCatch({
          # Get number of proteins
          Proteins <- unique(FrameTable$PROTEIN)
          # Get combinations
          Combinations <- combn(Proteins, 2)
          Tests <- 1:NCOL(Combinations)
          
          FindNeighbors <- function(PairX){
            tryCatch({
              # Proteins to run
              ReferenceProtein <- Combinations[1,PairX]
              QueryProtein <- Combinations[2,PairX]
              # Get tables
              ReferenceProtein <- FrameTable %>% filter(PROTEIN == ReferenceProtein)
              QueryProtein <- FrameTable %>% filter(PROTEIN == QueryProtein)
              ReferenceProteinCoordinates <- ReferenceProtein %>% select(POSITION_X, POSITION_Y)
              QueryProteinCoordinates <- QueryProtein %>% select(POSITION_X, POSITION_Y)
              
              Distance <- RANN::nn2(QueryProteinCoordinates,  ReferenceProteinCoordinates, searchtype = c("radius"), radius = Radius,k = 1)
              
              GetSpotAndDistanceFx <- function(DistanceX){
                tryCatch({
                  QuerySpotName <- Distance$nn.idx[DistanceX]
                  if(QuerySpotName != 0){
                    Result <- NULL
                    Result$REFERENCE_SPOT <- ReferenceProtein$UNIVERSAL_SPOT_ID[DistanceX]
                    Result$REFERENCE_TRACK <- ReferenceProtein$UNIVERSAL_TRACK_ID[DistanceX]
                    Result$QUERY_SPOT <- QueryProtein$UNIVERSAL_SPOT_ID[QuerySpotName]
                    Result$QUERY_TRACK <- QueryProtein$UNIVERSAL_TRACK_ID[QuerySpotName]
                    
                    Result$COLOCALIZATION_DISTANCE <- Distance$nn.dists[DistanceX]
                    Result$ASSOCIATION_TIME_THRESHOLD <- AssocTime
                    Result$REFERENCE_PROTEIN = ReferenceProtein$PROTEIN[1]
                    Result$QUERY_PROTEIN = QueryProtein$PROTEIN[1]
                    # Result$UNIVERSAL_SPOT_ID <- c(QueryProtein$UNIVERSAL_SPOT_ID[QuerySpotName], ReferenceProtein$UNIVERSAL_SPOT_ID[DistanceX])
                    # Result$COLOCALIZATION_SPOT <- c(ReferenceProtein$UNIVERSAL_SPOT_ID[DistanceX], QueryProtein$UNIVERSAL_SPOT_ID[QuerySpotName])
                    # Result$DISTANCE <- rep(Distance$nn.dists[DistanceX], 2)
                    Result <- as.data.frame(Result)
                    return(Result)
                  }
                }, error = function(e){print(paste("ERROR with GetSpotAndDistanceFx. DistanceX =", DistanceX))})
              }
              Distance <- lapply(1:NROW(Distance$nn.idx), GetSpotAndDistanceFx)
              Distance <- Distance[(which(sapply(Distance,is.list), arr.ind=TRUE))]
              Distance <- rbindlist(Distance)
              
              return(Distance)
            }, error = function(e){print(paste("ERROR with FindNeighbors. PairX = ", PairX))})
          }
          TableDistances <- lapply(Tests, FindNeighbors)
          TableDistances <- rbindlist(TableDistances, fill = TRUE)
          
          # Add cell identifier
          if(NROW(TableDistances) > 0){
            TableDistances$IMAGE_CELL = paste(FrameTable$IMAGE[1], FrameTable$CELL[1], sep = "...")
            TableDistances$RELATIVE_PATH = FrameTable$RELATIVE_PATH[1]
            TableDistances$DISTANCE_THRESHOLD <- Radius
            TableDistances$FRAME <- FrameTable$FRAME[1]
          }
          return(TableDistances)
        }, error = function(e){print(paste("ERROR with CostFx"))})
      }
      CostTable <- lapply(CellTable, CostFx)
      CostTable <- CostTable[(which(sapply(CostTable,is.list), arr.ind=TRUE))]
      CostTable <- rbindlist(CostTable, fill = TRUE)
      
      if(NROW(CostTable) != 0){
        
        # Remove spots shared between reference/query
        CostTable <-
          CostTable %>% 
          group_by(
            QUERY_PROTEIN,
            REFERENCE_SPOT
          ) %>% 
          filter(
            COLOCALIZATION_DISTANCE == min(COLOCALIZATION_DISTANCE)
          ) %>% 
          group_by(
            REFERENCE_PROTEIN,
            QUERY_SPOT
          ) %>% 
          filter(
            COLOCALIZATION_DISTANCE == min(COLOCALIZATION_DISTANCE)
          ) %>% 
          ungroup() %>% 
          mutate(
            PAIR_ID = 1:n()
          )
        
        Edges <-
          CostTable %>% 
          arrange(
            REFERENCE_TRACK,
            QUERY_TRACK,
            FRAME
          ) %>% 
          group_by(
            REFERENCE_TRACK,
            QUERY_TRACK
          ) %>%
          mutate(
            ASSOCIATION_TIME = lead(FRAME)-FRAME
          ) %>% 
          filter(
            ASSOCIATION_TIME == 1
          ) %>%
          mutate(
            ASSOCIATION_TIME = ifelse(lead(FRAME)-FRAME == 1, cumsum(ASSOCIATION_TIME), 0),
            ASSOCIATION_TIME = max(ASSOCIATION_TIME, na.rm = T),
          ) %>% 
          filter(
            ASSOCIATION_TIME > (ASSOCIATION_TIME_THRESHOLD-1),
          ) %>% 
          ungroup()
        
        EdgePairs <-
          Edges %>% 
          select(
            PAIR_ID
          ) %>% 
          mutate(
            COLOCALIZATION_INCLUDED = TRUE
          )
        
        CostTable <- as.data.frame(CostTable)
        EdgePairs <- as.data.frame(EdgePairs)
        CostTable <- merge(CostTable, EdgePairs, by = "PAIR_ID", all = TRUE, fill = TRUE)
        CostTable$PAIR_ID <- NULL
        CostTable <- CostTable %>% mutate(COLOCALIZATION_INCLUDED = ifelse(is.na(COLOCALIZATION_INCLUDED), FALSE, TRUE))
        
        # Get protein names
        Proteins <- unique(c(CostTable$REFERENCE_PROTEIN, CostTable$QUERY_PROTEIN))
        Proteins <- expand.grid(Proteins, Proteins)
        names(Proteins) <- c("PROTEIN", "COLOCALIZATION_PROTEIN")
        Proteins <- Proteins[Proteins$PROTEIN != Proteins$COLOCALIZATION_PROTEIN,]
        # Get pair ID
        Proteins <-
          Proteins %>% 
          arrange(
            PROTEIN,
            COLOCALIZATION_PROTEIN
          ) %>% 
          group_by(
            PROTEIN
          ) %>% 
          mutate(
            PAIR = paste(PROTEIN, COLOCALIZATION_PROTEIN, sep = "..."),
            COLUMN_NUMBER = 1:n()
          ) %>% 
          ungroup()
        
        # Get reference spots
        CostTableReferenceExport <- 
          CostTable %>% 
          mutate(
            UNIVERSAL_SPOT_ID = REFERENCE_SPOT,
            COLOCALIZATION_SPOT = QUERY_SPOT,
            PROTEIN = REFERENCE_PROTEIN,
            COLOCALIZATION_PROTEIN = QUERY_PROTEIN
          ) %>% 
          select(
            UNIVERSAL_SPOT_ID,
            PROTEIN,
            COLOCALIZATION_PROTEIN,
            UNIVERSAL_SPOT_ID,
            COLOCALIZATION_SPOT,
            COLOCALIZATION_DISTANCE,
            COLOCALIZATION_INCLUDED
          )
        
        # Get query spots
        CostTableQueryExport <- 
          CostTable %>% 
          mutate(
            UNIVERSAL_SPOT_ID = QUERY_SPOT,
            COLOCALIZATION_SPOT = QUERY_SPOT,
            PROTEIN = QUERY_PROTEIN,
            COLOCALIZATION_PROTEIN = REFERENCE_PROTEIN
          ) %>% 
          select(
            UNIVERSAL_SPOT_ID,
            PROTEIN,
            COLOCALIZATION_PROTEIN,
            COLOCALIZATION_SPOT,
            COLOCALIZATION_DISTANCE,
            COLOCALIZATION_INCLUDED
          )
        
        CostTableExport <- rbind(CostTableReferenceExport, CostTableQueryExport)
        remove(CostTableReferenceExport, CostTableQueryExport)
        # Add column number
        SplitNames <- function(PairX){
          tryCatch({
            # Get column number
            TempProteins <-
              Proteins %>% 
              filter(
                PAIR == Proteins$PAIR[PairX]
              )
            ColumnNumber <- TempProteins$COLUMN_NUMBER[1]
            # Get pair table
            TempCostTable <-
              CostTableExport %>% 
              filter(
                PROTEIN == Proteins$PROTEIN[PairX],
                COLOCALIZATION_PROTEIN == Proteins$COLOCALIZATION_PROTEIN[PairX]
              ) %>% 
              select(-c(
                PROTEIN
              ))
            # Rename columns
            names(TempCostTable) <- 
              c(
                "UNIVERSAL_SPOT_ID",
                paste0("COLOCALIZATION_PROTEIN_", ColumnNumber),
                paste0("COLOCALIZATION_SPOT_", ColumnNumber),
                paste0("COLOCALIZED_DISTANCE_", ColumnNumber),
                paste0("COLOCALIZATION_INCLUDED_", ColumnNumber)
              )
            return(TempCostTable)
          }, error = function(e){print(paste("ERROR with SplitNames"))})
        }
        CostTableExport <- lapply(1:NROW(Proteins), SplitNames)
        CostTableExport <- CostTableExport[(which(sapply(CostTableExport,is.list), arr.ind=TRUE))]
        CostTableExport <- rbindlist(CostTableExport, fill = TRUE)
        
        # Get names of columns to fill
        FillNames <- c("COLOCALIZATION_PROTEIN_",
                       "COLOCALIZATION_SPOT_",
                       "COLOCALIZED_DISTANCE_",
                       "COLOCALIZATION_INCLUDED_")
        
        FillNames <- expand.grid(FillNames, unique(Proteins$COLUMN_NUMBER))
        FillNames <- paste0(FillNames$Var1, FillNames$Var2)
        
        # Put table together
        CostTableExport <-
          CostTableExport %>% 
          group_by(
            UNIVERSAL_SPOT_ID
          ) %>% 
          fill(
            all_of(FillNames),
            .direction = "updown"
          ) %>% 
          distinct() %>% 
          mutate(
            DISTANCE_THRESHOLD = Radius,
            ASSOCIATION_TIME_THRESHOLD = AssocTime
          )
        
        # Save colocalized spots
        save_path = file.path(Cells[CellX], "colocalization_cost.csv.gz")
        file.remove(save_path, showWarnings = F)
        fwrite(CostTableExport, save_path, row.names = F, na = "")
        remove(CostTableExport)
        
        # Get edges (i.e., the colocalized spot pair)
        Edges <-
          Edges %>% 
          select(
            REFERENCE_TRACK,
            QUERY_TRACK
          ) %>% 
          distinct()
        
        # Get group numbers 
        names(Edges) <- c("from", "to")
        SubGroups <- igraph::graph_from_data_frame(Edges, directed = F)
        SubGroups <- igraph::clusters(SubGroups)$membership
        SubGroups <- as.data.frame(SubGroups)
        
        # Create new table from result
        Groups <- NULL
        Groups$UNIVERSAL_TRACK_ID <- rownames(SubGroups)
        Groups$COLOCALIZATION_GROUP <- SubGroups$SubGroups
        Groups <- as_tibble(Groups)
        
        # Save colocalized group
        save_path = file.path(Cells[CellX], "colocalization_groups.csv.gz")
        file.remove(save_path, showWarnings = F)
        fwrite(Groups, save_path, row.names = F, na = "")
        
        return(save_path)
      }
    }
  }, error = function(e){print(paste("ERROR with CellColocalizationFx CellX =", CellX))})
}
CostTablePaths <- mclapply(1:NROW(Cells), CellColocalizationFx, mc.cores = detectCores(logical = F))
CostTablePaths <- unlist(CostTablePaths)
CostTablePaths <- CostTablePaths[file.exists(CostTablePaths)]
print('Colocalization by coordinates complete')
# Move images that needed colocalization
ImagePath <- ColocalizationNeeded$image
for(Image in ImagePath){
  old_path = Image
  Cohort = basename(dirname(Image))
  Image = basename(Image)
  new_path = file.path(colocalized_path, Cohort, Image)
  file.move(old_path, new_path)
}

# Move images that don't need colocalization
ImagePath <- NoColocalizationNeeded$image
for(Image in ImagePath){
  old_path = Image
  Cohort = basename(dirname(Image))
  Image = basename(Image)
  new_path = file.path(colocalized_path, Cohort, Image)
  file.move(old_path, new_path)
}
