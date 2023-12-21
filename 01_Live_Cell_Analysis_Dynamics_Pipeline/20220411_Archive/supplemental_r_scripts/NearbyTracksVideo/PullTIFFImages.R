setwd(RETURN_TO_DIRECTORY)

# Get image 
if(file.exists(paste0(lp.PROTEIN, "_MED.tif"))){
  str_name = paste0(lp.PROTEIN, "_MED.tif")
} else{
  str_name = paste0(lp.PROTEIN, ".tif")
}

#Import image
img = ijtiff::read_tif(str_name)
Frames = ijtiff::frames_count(str_name)[1]
FrameFx <- function(FrameX){
  tryCatch({
    
    Z = img[,,,FrameX]
    X = NROW(Z)
    X = rep(1:X, NCOL(Z))
    
    Y = NCOL(Z)
    Y = rep(1:Y, each = NROW(Z))
    
    Z = as.vector(Z)
    
    plot_img = cbind(X, Y, Z)
    plot_img = as_tibble(plot_img)
    
    plot_img <-
      plot_img %>%
      mutate(
        y = X*PIXEL_SIZE,
        x = Y*PIXEL_SIZE,
        t = FrameX,
        z = Z/ExpTracks$PROTEIN_INTENSITY[1]*9
      ) %>%
      select(-c(
        X,
        Y,
        Z
      ))
    plot_img
  }, error = function(e){print(paste("     ERROR with FrameFx FrameX =", FrameX))})
}
plot_img <- mclapply(1:Frames, FrameFx)
plot_img <- data.table::rbindlist(plot_img)

SUBFOLDER <- file.path(FOLDER, "Track_Videos")
dir.create(SUBFOLDER)
