# find ./ -type d | less >> Folders.txt

library(data.table)
library(tidyverse)

ImgList <- read.delim("Folders.txt")
names(ImgList) <- "Path"


ImgList <-
  ImgList %>% 
  mutate(
    IMAGE = basename(dirname(Path)),
    COHORT = dirname(dirname(Path))
    
    # IMAGE = basename((Path)),
    # COHORT = "Calibrations"
  ) %>% 
  filter(
    IMAGE != ".",
    COHORT != "."
  ) %>% 
  select(-c(
    Path
  )) %>% 
  mutate(
    IMAGE = gsub("./", "", IMAGE),
    COHORT = gsub("./", "", COHORT)
  ) %>% 
  distinct() %>% 
  mutate(
    NOTE = ""
  )

write.csv(ImgList, "CombineImages.csv")


Imgs1 <- file.path("/raven/u/deliz/new_pipeline/Output4/Calibrations", ImgList$IMAGE, "Cell_1/GFP_intensity.csv.gz")
Imgs2 <- file.path("/raven/u/deliz/new_pipeline/Output4/Calibrations", ImgList$IMAGE, "Cell_1/mScarlet_intensity.csv.gz")
Imgs <- bind_rows(Imgs1, Imgs2)


Imgs <- lapply(Imgs, fread)


