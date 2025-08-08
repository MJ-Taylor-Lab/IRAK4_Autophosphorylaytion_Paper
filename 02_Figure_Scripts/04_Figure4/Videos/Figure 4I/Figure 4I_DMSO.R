library(pacman)
# Comment: `pacman` is an excellent package for managing other R packages.
# `p_load` checks if a package is installed and loads it; if not, it installs and then loads it.
# This makes the script more portable as it handles dependencies automatically.
pacman::p_load(EBImage)

Image_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/00_New_Figures/Figure_4/Video/IRAK3-KinaseDomain/DMSO/MyD88_crop.tif"

image <- readImage(Image_Path)

# Replicate the image 400 times
stack <- replicate(301, image, simplify = FALSE)
# Combine into a stack
stack <- combine(stack)


Image_Path <- dirname(Image_Path)
Save_Path <- file.path(Image_Path, "MyD88_box_crop_stack.tif")
writeImage(stack, Save_Path, type = "tiff")



# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()

