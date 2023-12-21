#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# Import variable
ImageListPath = args[1]

# Import libraries
library(dplyr)
library(data.table)
library(R.utils)
library(parallel)

# Script for combining image tables
# To have a table of everything, merge Analysis and Parameters by RELATIVE_PATH
TablesToCombine <- c("Parameters.csv.gz", "Essential.csv.gz", "Analysis.csv.gz")

# Read table
TablesList <- fread(ImageListPath)

# Remove extension if added
TablesList <-
  TablesList %>% 
  mutate(
    IMAGE = gsub(".nd2", "", IMAGE)
  )

# Get path to images
ImagesPath <- dirname(ImageListPath)
print(paste("Images Path:", ImagesPath))
setwd(ImagesPath)

# Iterate over table names to combine
CombineFx <- function(OutputX){
  print(TablesToCombine[OutputX])
  
  # Pull all tables and merge
  TempTablesList <- TablesList
  TempTablesList$PATH <-
    file.path(ImagesPath, TablesList$COHORT, TablesList$IMAGE, TablesToCombine[OutputX])
  
  Files <- TempTablesList$PATH
  Files <- Files[file.exists(Files)]
  print(paste("Total Images =", NROW(Files)))

  ReadTable <- function(ImageX){
    tryCatch({
      print(paste("Importing ImageX =" , ImageX))
      
      if(grepl("macOS", osVersion)){
        TempTable <- fread(Files[ImageX], fill=TRUE)
      } else{
        # Import
        TempPath <- gsub(".gz", "", Files[ImageX])
        if(file.exists(TempPath)){
          file.remove(TempPath)
        }
        gunzip(Files[ImageX], remove = FALSE)
        TempTable <- fread(file = file.path(TempPath))
        # Delete file
        file.remove(TempPath)
      }
      return(TempTable)
    }, error = function(e){print(paste("Error with ImageX =", ImageX))})
  }
  GrandTable <- mclapply(1:NROW(Files), ReadTable)
  GrandTable <- GrandTable[(which(sapply(GrandTable,is.list), arr.ind=TRUE))]
  GrandTable <- rbindlist(GrandTable, use.names = TRUE, fill = TRUE)
  # Add notes to the images
  TempTablesList <- TempTablesList %>% dplyr::select(IMAGE, NOTE)
  GrandTable <- as.data.frame(GrandTable)
  TempTablesList <- as.data.frame(TempTablesList)
  tryCatch({
    GrandTable <- merge(GrandTable, TempTablesList, by = "IMAGE")
  }, error = function(e){print("Cannot add notes")})
  # Save
  print(paste("Saving", TablesToCombine[OutputX]))
  file.remove(file.path(ImagesPath, TablesToCombine[OutputX]))
  fwrite(GrandTable, file.path(ImagesPath, TablesToCombine[OutputX]))
}
lapply(1:NROW(TablesToCombine), CombineFx)

print("Completed GrandTable")
