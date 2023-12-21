library(data.table)
library(parallel)
library(tidyverse)

Experimental <- fread("/Users/u_deliz/Desktop/PhasePortraitAnalysis/Essential.csv.gz")
Calibrations <- fread("/Users/u_deliz/Desktop/PhasePortraitAnalysis/Fluorophores.csv.gz")

Experimental <-
  Experimental %>% 
  mutate(
    FLUOROPHORE = ifelse(PROTEIN == "MyD88", "GFP", "mScarlet")
    DATE =  as.Date(substr(IMAGE, 1, 8), format = '%Y%m%d')
  )

CalibrationsSummary <- 
  Calibrations %>% 
  filter(
    SPOT_AREA == 25
  ) %>% 
  group_by(
    IMAGE,
    PROTEIN
  ) %>% 
  summarize(
    DATE = substr(IMAGE, 1, 8),
    DATE =  as.Date(DATE, format = '%Y%m%d'),
    N = n(),
    TOTAL_INTENSITY = median(TOTAL_INTENSITY)
  ) %>% 
  distinct() %>% 
  group_by(
    DATE,
    PROTEIN
  ) %>% 
  mutate(
    RANK = rank(-N)
  ) %>% 
  filter(
    RANK == 1
  ) %>% 
  select(-c(
    RANK
  ))

ImgToCalibrate <-
  Experimental %>% 
  select(
    IMAGE,
    FLUOROPHORE,
    DATE
  ) %>% 
  distinct()

NormalizeIntensitites <- function(ImageX){
  
  TempExperimental <-
    Experimental %>% 
    filter(
      IMAGE == ImgToCalibrate$IMAGE[ImageX],
      FLUOROPHORE == ImgToCalibrate$FLUOROPHORE[ImageX],
      DATE == ImgToCalibrate$DATE[ImageX]
    )
  
  TempCalibrationsSummary <-
    CalibrationsSummary %>% 
    filter(
      PROTEIN == TempExperimental$FLUOROPHORE[1]
    )
  
  DateFiltered = TempCalibrationsSummary[which.min(abs(TempCalibrationsSummary$DATE - TempExperimental$DATE[1])),]
  
  
  OtherChannel <-
    Experimental %>% 
    filter(
      IMAGE == ImgToCalibrate$IMAGE[ImageX],
      FLUOROPHORE != ImgToCalibrate$FLUOROPHORE[ImageX],
      DATE == ImgToCalibrate$DATE[ImageX]
    ) %>% 
    select(
      FLUOROPHORE,
      DATE
    ) %>% 
    distinct()
  
  TempCalibrationsSummary <-
    CalibrationsSummary %>% 
    filter(
      PROTEIN == OtherChannel$FLUOROPHORE[1]
    )
  
  OtherDateFiltered = TempCalibrationsSummary[which.min(abs(TempCalibrationsSummary$DATE - OtherChannel$DATE)),]
  
  
  TempExperimental <-
    TempExperimental %>% 
    mutate(
      NORMALIZED_INTENSITY = NORMALIZED_INTENSITY/DateFiltered$TOTAL_INTENSITY,
      COMPLEMENTARY_NORMALIZED_INTENSITY_1 = COMPLEMENTARY_NORMALIZED_INTENSITY_1/OtherDateFiltered$TOTAL_INTENSITY
    )
  
  return(TempExperimental)
  
}

Table <- mclapply(1:NROW(ImgToCalibrate), NormalizeIntensitites)
Table <- rbindlist(Table)

fwrite(Table, "Essential.csv.gz")

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