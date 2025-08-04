setwd(RETURN_TO_DIRECTORY)

PlotFx <- function(TrackX){
  # tryCatch({
  FOCUS_TRACK = SubTracks[TrackX]
  
  TRACK_ID = strsplit(FOCUS_TRACK, '...', fixed = T)[[1]][4]
  
  TRACK_FOLDER <- file.path(SUBFOLDER, TRACK_ID)
  dir.create(TRACK_FOLDER)
  setwd(TRACK_FOLDER)
  # Get table of track coordinates
  plot_spot <-
    ExpTracks %>%
    filter(
      UNIVERSAL_TRACK_ID == FOCUS_TRACK
    ) %>%
    mutate(
      t = FRAME
    ) %>%
    select(
      t,
      NORMALIZED_INTENSITY,
      POSITION_X,
      POSITION_Y
    )
  
  # Get intensity of five spots back
  plot_spot_back <- NULL
  plot_spot_back$t <- (min(plot_spot$t)-5):(min(plot_spot$t)-1)
  plot_spot_back <- as_tibble(plot_spot_back)
  plot_spot_back <-
    plot_spot_back %>%
    mutate(
      POSITION_X =  plot_spot$POSITION_X[1],
      POSITION_Y = plot_spot$POSITION_Y[1]
    ) %>%
    filter(
      t >= 1
    )
  
  # Combine tables
  plot_spot <- bind_rows(plot_spot_back, plot_spot)
  
  # Range of frames to show
  RANGE <- min(plot_spot$t):max(plot_spot$t)
  
  # Add missing spots
  # and middle point between coordinates
  if(NROW(plot_spot) != NROW(RANGE)){
    
    MISSING_T <- RANGE[which(!RANGE %in% unique(plot_spot$t))]
    
    plot_spot_missing <- NULL
    plot_spot_missing$t <- MISSING_T
    plot_spot_missing <- as_tibble(plot_spot_missing)
    
    plot_spot <- bind_rows(plot_spot_missing, plot_spot)
    
    plot_spot <-
      plot_spot %>%
      arrange(
        t
      ) %>%
      mutate(
        POSITION_X = ifelse(is.na(POSITION_X), (lag(POSITION_X) + lead(POSITION_X))/2, POSITION_X),
        POSITION_Y = ifelse(is.na(POSITION_Y), (lag(POSITION_Y) + lead(POSITION_Y))/2, POSITION_Y)
      )
  }
  
  plot_img_spot <- merge(plot_spot, plot_img, by = "t")
  
  plot_img_spot <-
    plot_img_spot %>%
    group_by(
      t
    ) %>%
    mutate(
      x = x - POSITION_X,
      x = x / PIXEL_SIZE,
      x = round(x),
      y = y - POSITION_Y,
      y = y / PIXEL_SIZE,
      y = round(y),
      tadj = t-min(t),
      t = t
    ) %>%
    filter(
      x >= -20,
      y >= -20,
      x <= 20,
      y <= 20
    )
  
  COLUMNS <- ceiling(NROW(RANGE)/5)
  ROWS <- ceiling(NROW(RANGE)/COLUMNS)
  ROWS <- ifelse(ROWS>COLUMNS, COLUMNS, ROWS)
  
  plot_all_img_spot <-
    plot_img_spot %>%
    filter(
      x >= -12,
      y >= -12,
      x <= 12,
      y <= 12
    )
  
  RETURN_TO_DIRECTORY <- getwd()
  setwd(NEARBY_TRACKS_SCRIPT)
  print("Running script PlotAllSpots.R")
  source("PlotAllSpots.R", local = T)
  
  PlotFrameFx <- function(TimeX) {
    tryCatch({
      FOCUS_TIME = RANGE[TimeX]
      
      plot_tracks <-
        ExpTracks %>%
        group_by(
          UNIVERSAL_TRACK_ID
        ) %>%
        mutate(
          KEEP = ifelse(FRAME==FOCUS_TIME, T, F),
          KEEP = sum(KEEP)
        ) %>%
        filter(
          KEEP > 0
        )
      
      RETURN_TO_DIRECTORY <- getwd()
      setwd(NEARBY_TRACKS_SCRIPT)
      print("Running script PlotCellView.R")
      source("PlotCellView.R", local = T)
      
      RETURN_TO_DIRECTORY <- getwd()
      setwd(NEARBY_TRACKS_SCRIPT)
      print("Running script PlotLimitedTracks.R")
      source("PlotLimitedTracks.R", local = T)
      
      RETURN_TO_DIRECTORY <- getwd()
      setwd(NEARBY_TRACKS_SCRIPT)
      print("Running script PlotTracks.R")
      source("PlotTracks.R", local = T)
      
      RETURN_TO_DIRECTORY <- getwd()
      setwd(NEARBY_TRACKS_SCRIPT)
      print("Running script PlotOtherSpots.R")
      source("PlotOtherSpots.R", local = T)
      
      RETURN_TO_DIRECTORY <- getwd()
      setwd(NEARBY_TRACKS_SCRIPT)
      print("Running script PlotSingleSpot.R")
      source("PlotSingleSpot.R", local = T)
      
      RETURN_TO_DIRECTORY <- getwd()
      setwd(NEARBY_TRACKS_SCRIPT)
      print("Running script PlotLimitedSingleSpot.R")
      source("PlotLimitedSingleSpot.R", local = T)

      AllPlots <-
        arrangeGrob(
          CellView, LimitedSingleSpot, SingleSpot, LimitedTracksPlot, TracksPlot, AllSpots,
          layout_matrix =
            rbind(c(1,2,4),
                  c(1,2,4),
                  c(1,3,5),
                  c(1,3,5),
                  c(6,6,6),
                  c(6,6,6),
                  c(6,6,6)
            )
        )

      setwd(TRACK_FOLDER)
      
      cowplot::ggdraw(AllPlots) +
        theme(plot.background = element_rect(fill="black", color = NA)) +
        ggsave(
          paste0(FOCUS_TIME, ".png"),
          height = 9*0.8,
          width = 16*0.8
        )
      pdf(NULL)
    }, error = function(e){print(paste("     ERROR with PlotFrameFx TimeX =", TimeX))})
  }
  
  ObjectList <- c(
    ObjectList,
    ls(envir=globalenv()),
    ls(envir=environment())
  )

  cl <- parallel::makeCluster(parallel::detectCores())
  doParallel::registerDoParallel(cl)
  foreach(
    TimeX = 1:NROW(RANGE),
    .packages = (.packages()),
    .export = ObjectList
  ) %dopar% {
    PlotFrameFx(TimeX)
    return(print(TimeX))
  }
  parallel::stopCluster(cl)
  
  FILE_LIST <-  paste0(RANGE, ".png")
  FILE_LIST <- file.path(TRACK_FOLDER, FILE_LIST)
  FILE_LIST <- FILE_LIST[file.exists(FILE_LIST)]
  
  setwd(SUBFOLDER)
  # Output name
  DISTANCE = which(ExpTracks$UNIVERSAL_TRACK_ID == FOCUS_TRACK)
  DISTANCE = min(DISTANCE)
  DISTANCE = ExpTracks$CLOSEST_DISTANCE[DISTANCE]
  DISTANCE = round(DISTANCE)
  # DISTANCE = formatC(DISTANCE, format = "f", digits = 1)
  DISTANCE = paste(DISTANCE, "px")
  
  OUTPUT_NAME = paste(DISTANCE, FOCUS_TRACK, sep = " | ")
  OUTPUT_NAME = paste0(OUTPUT_NAME, ".mp4")
  
  av::av_encode_video(
    FILE_LIST,
    framerate = VIDEO_FRAME_RATE,
    output = file.path(OUTPUT_FOLDER, OUTPUT_NAME)
  ) 
  
  # }, error = function(e){print(paste("     ERROR with PlotFx TrackX =", TrackX))})
}

for(TrackX in 1:NROW(SubTracks)){
  print(paste("     ", TrackX))
  PlotFx(TrackX)
}
