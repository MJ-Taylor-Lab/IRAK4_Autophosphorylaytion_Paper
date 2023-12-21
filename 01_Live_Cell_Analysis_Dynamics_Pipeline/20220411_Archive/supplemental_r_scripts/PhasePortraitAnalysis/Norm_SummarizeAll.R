setwd(OUTPUT_DIRECTORY)

# Table Name
SaveName <-
  paste0(
    "Normalized StatTable - ",
    "LeadLag ", LEAD_LAG, " - ",
    "StepSize ", round(STEP_SIZE, 2),
    ".csv.gz"
  )
SaveName <- gsub(" - \\.", "\\.", SaveName)

StatTable <- fread(file.path(OUTPUT_DIRECTORY, SaveName))

# Separate by landing time
LandingTimeFx <- function(ThresholdX){
  TempTable <- 
    StatTable %>%
    filter(
      FRAMES_SINCE_LANDING_CAT <= ThresholdX
    ) %>% 
    mutate(
      FRAMES_SINCE_LANDING_CAT = ThresholdX
    ) %>% 
    as.data.table()
  
  return(TempTable)
}
LandingTimes <- unique(StatTable$FRAMES_SINCE_LANDING_CAT)
StatTable <- mclapply(LandingTimes, LandingTimeFx, mc.cores = detectCores(logical = F))
StatTable <- rbindlist(StatTable)

SplitNormalization <-
  StatTable %>% 
  mutate(
    NORMALIZATION_GROUP = paste(COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN)
  ) %>% 
  as_tibble() %>% 
  group_split(
    NORMALIZATION_GROUP
  )

SplitNormalization <- SplitNormalization[order(sapply(SplitNormalization,nrow))]
SplitNormalization <- SplitNormalization[c(NROW(SplitNormalization):1)]

NormalizeFx <- function(TableX){
  
  Normalization <-
    TableX %>%
    group_by(
      IMAGE
    ) %>% 
    summarize(
      RATIO = mean(MAX_REFERENCE/MAX_QUERY),
      MIN_REFERENCE = mean(MIN_REFERENCE),
      MIN_QUERY = mean(MIN_QUERY),
      MAX_REFERENCE = mean(MAX_REFERENCE),
      MAX_QUERY = mean(MAX_QUERY)
    ) %>% 
    ungroup() %>% 
    mutate(
      FINAL_RATIO = median(RATIO),
    ) %>% 
    mutate(
      TEST = which.min(abs(RATIO - FINAL_RATIO)),
      ID = 1:n()
    ) %>% 
    filter(
      ID == TEST
    ) %>% 
    as_tibble()
  
  NormalizedTable <-
    TableX %>%
    mutate(
      MIN_REFERENCE = Normalization$MIN_REFERENCE[1],
      MAX_REFERENCE = Normalization$MAX_REFERENCE[1],
      MIN_QUERY = Normalization$MIN_QUERY[1],
      MAX_QUERY = Normalization$MAX_QUERY[1],
    ) %>%
    mutate(
      REFERENCE_TOTAL_INTENSITY = REFERENCE_TOTAL_INTENSITY*MAX_REFERENCE+MIN_REFERENCE,
      ROUNDED_REFERENCE_TOTAL_INTENSITY = round(REFERENCE_TOTAL_INTENSITY/STEP_SIZE)*STEP_SIZE,
      DELTA_REFERENCE_TOTAL_INTENSITY = DELTA_REFERENCE_TOTAL_INTENSITY*MAX_REFERENCE,
      ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY = ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY*MAX_REFERENCE,
      
      QUERY_TOTAL_INTENSITY = QUERY_TOTAL_INTENSITY*MAX_QUERY+MIN_QUERY,
      ROUNDED_QUERY_TOTAL_INTENSITY = round(QUERY_TOTAL_INTENSITY/STEP_SIZE)*STEP_SIZE,
      DELTA_QUERY_TOTAL_INTENSITY = DELTA_QUERY_TOTAL_INTENSITY*MAX_QUERY,
      ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY = ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY*MAX_QUERY
    ) %>%
    as.data.table()
  
  NormalizedTable$NORMALIZATION_GROUP = NULL
  
  return(NormalizedTable)
}
NormStatTable <- mclapply(SplitNormalization, NormalizeFx)
NormStatTable <- rbindlist(NormStatTable)
remove(SplitNormalization)

SplitStatTable <-
  NormStatTable %>% 
  arrange(
    COHORT, LIGAND_DENSITY_CAT
  ) %>% 
  group_by(
    COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN
  ) %>% 
  mutate(
    FPS = round(FPS, 2),
    FACET = paste(COHORT, LIGAND_DENSITY_CAT, FPS, FRAMES_SINCE_LANDING_CAT, REFERENCE_PROTEIN, QUERY_PROTEIN),
    PLOT_FACETS = paste(COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN)
  ) %>% 
  as_tibble() %>% 
  group_split(
    FACET
  )

SplitStatTable <- SplitStatTable[order(sapply(SplitStatTable,nrow))]
SplitStatTable <- SplitStatTable[c(NROW(SplitStatTable):1)]

SplitSummaryFx <- function(TableX){

  # Summarize results dealing with either reference OR query increasing
  PhasePortrait <-
    TableX %>% 
    filter(
      !is.infinite(DELTA_REFERENCE_TOTAL_INTENSITY),
      !is.infinite(DELTA_QUERY_TOTAL_INTENSITY),
      
      !is.infinite(ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY),
      !is.infinite(ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY)
    ) %>% 
    as_tibble() %>%
    group_by(
      COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN,
      FACET, PLOT_FACETS,
      ROUNDED_REFERENCE_TOTAL_INTENSITY, ROUNDED_QUERY_TOTAL_INTENSITY,
      FRAMES_SINCE_LANDING_CAT
    ) %>% 
    summarize(
      # Get number of spots in bin so that bins with few spots can be filtered out
      N = n(),
      FRAMES_ADJUSTED = median(FRAMES_ADJUSTED),
      TIME_ADJUSTED = median(TIME_ADJUSTED),
      
      # Deviation from median
      MAD_DELTA_REFERENCE_TOTAL_INTENSITY = mad(DELTA_REFERENCE_TOTAL_INTENSITY),
      MAD_DELTA_QUERY_TOTAL_INTENSITY = mad(DELTA_QUERY_TOTAL_INTENSITY),
      MAD_ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY = mad(ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY),
      MAD_ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY = mad(ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY),
      MAD_DELTA_RATIO = mad(DELTA_REFERENCE_TOTAL_INTENSITY/DELTA_QUERY_TOTAL_INTENSITY),
      MAD_ADJUSTED_DELTA_RATIO = mad(ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY/ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY),

      DELTA_REFERENCE_Q1 = quantile(DELTA_REFERENCE_TOTAL_INTENSITY, .25, na.rm = T),
      DELTA_REFERENCE_Q3 = quantile(DELTA_REFERENCE_TOTAL_INTENSITY, .75, na.rm = T),
      DELTA_QUERY_Q1 = quantile(DELTA_QUERY_TOTAL_INTENSITY, .25, na.rm = T),
      DELTA_QUERY_Q3 = quantile(DELTA_QUERY_TOTAL_INTENSITY, .75, na.rm = T),
      ADJUSTED_DELTA_REFERENCE_Q1 = quantile(ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY, .25, na.rm = T),
      ADJUSTED_DELTA_REFERENCE_Q3 = quantile(ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY, .75, na.rm = T),
      ADJUSTED_DELTA_QUERY_Q1 = quantile(ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY, .25, na.rm = T),
      ADJUSTED_DELTA_QUERY_Q3 = quantile(ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY, .75, na.rm = T),

      # Get median of intensity changes. mean could skew the data, thus not used
      DELTA_REFERENCE_TOTAL_INTENSITY = median(DELTA_REFERENCE_TOTAL_INTENSITY, na.rm = T),
      DELTA_QUERY_TOTAL_INTENSITY = median(DELTA_QUERY_TOTAL_INTENSITY, na.rm = T),
      
      ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY = median(ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY, na.rm = T),
      ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY = median(ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY, na.rm = T)
    ) %>% 
    group_by(
      COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN,
      FACET, PLOT_FACETS,
      ROUNDED_REFERENCE_TOTAL_INTENSITY, ROUNDED_QUERY_TOTAL_INTENSITY,
      FRAMES_SINCE_LANDING_CAT
    ) %>% 
    mutate(
      THETA = atan2(DELTA_QUERY_TOTAL_INTENSITY, DELTA_REFERENCE_TOTAL_INTENSITY),
      ADJUSTED_THETA = atan2(ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY, ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY),
      
      # Get angle to account for negative magnitudes
      THETA = atan2(-DELTA_REFERENCE_TOTAL_INTENSITY, -DELTA_QUERY_TOTAL_INTENSITY),
      ADJUSTED_THETA = atan2(-ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY, -ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY),
      
      DELTA_RATIO = DELTA_REFERENCE_TOTAL_INTENSITY/DELTA_QUERY_TOTAL_INTENSITY,
      ADJUSTED_DELTA_RATIO = ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY/ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY,
    ) %>% 
    mutate(
      # Make angle
      THETA = (THETA*180/pi+180)/360,
      ADJUSTED_THETA = (ADJUSTED_THETA*180/pi+180)/360,
      
      DELTA_RATIO = ifelse(DELTA_RATIO<1, -1/DELTA_RATIO, DELTA_RATIO),
      ADJUSTED_DELTA_RATIO = ifelse(ADJUSTED_DELTA_RATIO<1, -1/ADJUSTED_DELTA_RATIO, ADJUSTED_DELTA_RATIO)
    ) %>% 
    filter(
      N >= 5
    ) %>% 
    group_by(
      COHORT, LIGAND_DENSITY_CAT, FPS, REFERENCE_PROTEIN, QUERY_PROTEIN,
      FRAMES_SINCE_LANDING_CAT,
      FACET, PLOT_FACETS
    ) %>% 
    mutate(
      # Keep only N's above the 15th percentile
      N_TEST = ifelse(N >= 50, T, F),
      
      # # Get absolute magnitude so that the log can be taken later
      # MAD_MAGNITUDE = abs(MAD_DELTA_REFERENCE_TOTAL_INTENSITY) + abs(MAD_DELTA_QUERY_TOTAL_INTENSITY),
      # # Get angle to account for negative magnitudes
      # MAD_RAD_ANGLE = atan2(-MAD_DELTA_REFERENCE_TOTAL_INTENSITY, -MAD_DELTA_QUERY_TOTAL_INTENSITY),
      # # Make anglez
      # MAD_DEG_ANGLE = MAD_RAD_ANGLE*180/pi+180,
      # MAD_DEG_ANGLE = MAD_DEG_ANGLE/360,
      # 
      # # Get absolute magnitude so that the log can be taken later
      # MAD_ADJUSTED_MAGNITUDE = abs(MAD_ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY) + abs(MAD_ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY),
      # # Get angle to account for negative magnitudes
      # MAD_ADJUSTED_RAD_ANGLE = atan2(-MAD_ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY, -MAD_ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY),
      # # Make angle
      # MAD_ADJUSTED_DEG_ANGLE = MAD_ADJUSTED_RAD_ANGLE*180/pi+180,
      # MAD_ADJUSTED_DEG_ANGLE = MAD_ADJUSTED_DEG_ANGLE/360,
      # 
      # # # Get absolute magnitude so that the log can be taken later
      # MAGNITUDE = abs(DELTA_REFERENCE_TOTAL_INTENSITY) + abs(DELTA_QUERY_TOTAL_INTENSITY),
      # # Get angle to account for negative magnitudes
      # RAD_ANGLE = atan2(-DELTA_REFERENCE_TOTAL_INTENSITY, -DELTA_QUERY_TOTAL_INTENSITY),
      # # Make angle
      # DEG_ANGLE = RAD_ANGLE*180/pi+180,
      # DEG_ANGLE = DEG_ANGLE/360,
      # 
      # # Get absolute magnitude so that the log can be taken later
      # ADJUSTED_MAGNITUDE = abs(ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY) + abs(ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY),
      # # Get angle to account for negative magnitudes
      # ADJUSTED_RAD_ANGLE = atan2(-ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY, -ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY),
      # # Make angle
      # ADJUSTED_DEG_ANGLE = ADJUSTED_RAD_ANGLE*180/pi+180,
      # ADJUSTED_DEG_ANGLE = ADJUSTED_DEG_ANGLE/360,
    ) %>%
    as.data.table()
    
  PhasePortrait$QUERY_PROTEIN = factor(PhasePortrait$QUERY_PROTEIN, levels = PROTEIN_ORDER)
  return(PhasePortrait)
}
LinearPhasePortrait <- mclapply(SplitStatTable, SplitSummaryFx, mc.cores = detectCores(logical = F))
LinearPhasePortrait <- rbindlist(LinearPhasePortrait)

SaveName <-
  paste0(
    "Normalized All PhasePortrait - ",
    "LeadLag ", LEAD_LAG, " - ",
    "StepSize ", STEP_SIZE,
    ".csv.gz"
  )

SaveName <- gsub(" - \\.", "\\.", SaveName)

file.remove(file.path(OUTPUT_DIRECTORY, SaveName))
fwrite(LinearPhasePortrait, file.path(OUTPUT_DIRECTORY, SaveName))
