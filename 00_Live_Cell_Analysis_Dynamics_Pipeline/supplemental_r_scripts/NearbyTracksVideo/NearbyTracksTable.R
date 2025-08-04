setwd(RETURN_TO_DIRECTORY)

# Loop over tracks
DistanceTrackFx <- function(TrackX){
  tryCatch({
    # Tracks to work on
    FOCUS_TRACK = SubTracks[TrackX]
    # Generate table for reference track
    Reference <-
      ExpTracks %>%
      filter(
        UNIVERSAL_TRACK_ID == FOCUS_TRACK
      ) %>%
      select(
        UNIVERSAL_TRACK_ID,
        FRAME,
        POSITION_X,
        POSITION_Y
      )
    
    # Generate table for query track
    Query <-
      ExpTracks %>%
      filter(
        UNIVERSAL_TRACK_ID != FOCUS_TRACK
      ) %>%
      select(
        UNIVERSAL_TRACK_ID,
        FRAME,
        POSITION_X,
        POSITION_Y
      )
    # Generate frame range
    START <- min(Reference$FRAME)
    END <- max(Reference$FRAME)
    FRAMES <- START:END
    # Find nearest spot ± 3 frames
    DistanceFrameFx <- function(FrameX){
      tryCatch({
        FOCUS_FRAME = FRAMES[FrameX]
        # Reference spot
        FrameReference <-
          Reference %>%
          filter(
            FRAME == FOCUS_FRAME
          ) %>%
          select(-c(
            UNIVERSAL_TRACK_ID,
            FRAME
          ))
        
        # Query spots (±0 frames)
        FrameQuery <-
          Query %>%
          filter(
            FRAME >= FOCUS_FRAME - 0 &
              FRAME <= FOCUS_FRAME + 0
          ) %>%
          select(-c(
            FRAME
          ))
        
        TRACK_ORDER <- FrameQuery$UNIVERSAL_TRACK_ID
        FrameQuery$UNIVERSAL_TRACK_ID <- NULL
        
        if(NROW(FrameReference) == 0|NROW(FrameQuery) == 0){
          FrameReference$CLOSEST_DISTANCE = NA
          FrameReference$CLOSEST_TRACK = NA
          FrameReference$UNIVERSAL_TRACK_ID = NULL
          FrameReference
        } else {
          #Find nearest neighbor
          Pair <- as.data.frame(RANN::nn2(FrameReference, FrameQuery, k=1))
          # Make it a table
          CLOSEST_TRACK <- which(Pair$nn.dists==min(Pair$nn.dists))
          CLOSEST_TRACK <- TRACK_ORDER[CLOSEST_TRACK]
          CLOSEST_DISTANCE <- min(Pair$nn.dists)
          
          FrameReference$CLOSEST_DISTANCE <- CLOSEST_DISTANCE
          FrameReference$CLOSEST_TRACK <- CLOSEST_TRACK
          
          FrameReference$UNIVERSAL_TRACK_ID <- NULL
        }
        FrameReference
        
        
      }, error = function(e){print(paste("     ERROR with DistanceFrameFx. FrameX =", FrameX))})
    }
    Distances <- mclapply(1:NROW(FRAMES), DistanceFrameFx)
    Distances <- Distances[(which(sapply(Distances,is.list), arr.ind=TRUE))]
    Distances <- data.table::rbindlist(Distances)
    
    if(is.infinite(min(Distances$CLOSEST_DISTANCE, na.rm = T))){
      Distances <-
        Distances %>%
        mutate(
          UNIVERSAL_TRACK_ID = FOCUS_TRACK
        ) %>%
        select(
          UNIVERSAL_TRACK_ID,
          CLOSEST_DISTANCE,
          CLOSEST_TRACK
        )
      Distances <- Distances[1]
    } else {
      # Only nearest point for entire track
      Distances <-
        Distances %>%
        select(
          CLOSEST_DISTANCE,
          CLOSEST_TRACK
        ) %>%
        filter(
          CLOSEST_DISTANCE == min(Distances$CLOSEST_DISTANCE, na.rm = T)
        ) %>%
        mutate(
          CLOSEST_DISTANCE = CLOSEST_DISTANCE/PIXEL_SIZE,
          CLOSEST_DISTANCE = round(CLOSEST_DISTANCE),
          UNIVERSAL_TRACK_ID = FOCUS_TRACK
        )
    }
    
    Distances <- Distances %>% select(UNIVERSAL_TRACK_ID, CLOSEST_DISTANCE, CLOSEST_TRACK)
    Distances
  }, error = function(e){print(paste("     ERROR with DistanceTrackFx. TrackX =", TrackX))})
}
Distances <- mclapply(1:NROW(SubTracks), DistanceTrackFx)
Distances <- Distances[(which(sapply(Distances,is.list), arr.ind=TRUE))]
Distances <- data.table::rbindlist(Distances)
Distances <- Distances %>% distinct()
