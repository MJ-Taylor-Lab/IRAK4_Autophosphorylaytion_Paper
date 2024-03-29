#!/usr/bin/env Rscript

remove(list = ls())
gc(reset = TRUE)
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Import variables
args = commandArgs(trailingOnly=TRUE)
# args = strsplit("/raven/u/deliz/new_pipeline/Output4 /raven/u/deliz/dynamics_pipeline/supplemental_r_scripts/PhasePortraitAnalysis /raven/u/deliz/new_pipeline/Output4/PhasePortraitAnalysis 5 0.05 100 200 F", " ")
if(grepl("macOS", osVersion)){
  args = strsplit("/Users/u_deliz/Desktop/PhasePortraitAnalysis /Users/u_deliz/dynamics_pipeline/supplemental_r_scripts/PhasePortraitAnalysis /Users/u_deliz/Desktop/PhasePortraitAnalysis 5 2 100 200 0.01 0.95 0", " ")
}
args = unlist(args)

# Table location
TABLE_PATH = args[1]
print(paste("TABLE_PATH =", TABLE_PATH))

# Script location
SCRIPTS_DIRECTORY = args[2]
print(paste("SCRIPTS_DIRECTORY =", SCRIPTS_DIRECTORY))

# Output folder
OUTPUT_DIRECTORY = args[3]
print(paste("OUTPUT_DIRECTORY =", OUTPUT_DIRECTORY))

# Lead/Lag factor
LEAD_LAG = args[4] #frame(s) before and __ frame(s) after
LEAD_LAG = as.numeric(LEAD_LAG)
print(paste("LEAD_LAG =", LEAD_LAG))

# Grid size for binning and plots
STEP_SIZE = args[5] # molecules
STEP_SIZE = as.numeric(STEP_SIZE)
print(paste("STEP_SIZE =", STEP_SIZE))

MAX_FRAMES = args[6] #FALSE for linear scale
MAX_FRAMES = as.numeric(MAX_FRAMES)
print(paste("MAX_FRAMES =", MAX_FRAMES))

FRAMES_SINCE_LANDING_THRESHOLD = args[7] #FALSE for linear scale
FRAMES_SINCE_LANDING_THRESHOLD = as.numeric(FRAMES_SINCE_LANDING_THRESHOLD)
print(paste("FRAMES_SINCE_LANDING_THRESHOLD =", FRAMES_SINCE_LANDING_THRESHOLD))

LOWER_BOUNDARY = args[8] #FALSE for linear scale
LOWER_BOUNDARY = as.numeric(LOWER_BOUNDARY)
print(paste("LOWER_BOUNDARY =", LOWER_BOUNDARY))

UPPER_BOUNDARY = args[9] #FALSE for linear scale
UPPER_BOUNDARY = as.numeric(UPPER_BOUNDARY)
print(paste("UPPER_BOUNDARY =", UPPER_BOUNDARY))

MOVING_AVERAGE_WINDOW = args[10] #FALSE for linear scale
MOVING_AVERAGE_WINDOW = as.numeric(MOVING_AVERAGE_WINDOW)
print(paste("MOVING_AVERAGE_WINDOW =", MOVING_AVERAGE_WINDOW))


# Import libraries
pacman::p_load(dtplyr, lemon, ggquiver, ggplot2, ggdark, scales, parallel, ggforce, data.table, metR, viridis, RcppRoll, tidyr, dplyr, compiler)
filter <- dplyr::filter

if(
  !file.exists(OUTPUT_DIRECTORY)){
  dir.create(OUTPUT_DIRECTORY)
}

# OTHER
# Reference protein
USE_REFERENCE_PROTEIN = "MyD88"
# Order of protein appearance
PROTEIN_ORDER = c("MyD88", "IRAK4", "IRAK1", "TRAF6", "Pellino", "HOIL1", "NEMO", "A20")

# Run Setup----
tryCatch({
    print(":::::::::::::::::::: Setup.R ::::::::::::::::::::")
    setwd(SCRIPTS_DIRECTORY)
    print("Running script Norm_Setup.R")
    source("Norm_Setup.R", local = T)
}, error = function(e) {print("Error with Norm_Setup.R")})

# Run Analysis----
tryCatch({
  print(":::::::::::::::::::: Analysis.R ::::::::::::::::::::")
  setwd(SCRIPTS_DIRECTORY)
  print("Running script Norm_Analysis.R")
  source("Norm_Analysis.R", local = T)
}, error = function(e) {print("Error with Norm_Analysis.R")})

# Run Summary----
tryCatch({
  print(":::::::::::::::::::: Summarize.R ::::::::::::::::::::")
  setwd(SCRIPTS_DIRECTORY)
  print("Running script Norm_SummarizeAll.R")
  source("Norm_SummarizeAll.R", local = T)
}, error = function(e) {print("Error with Norm_SummarizeAll.R")})

# Plot Phase Portrait Linear----
tryCatch({
  print(":::::::::::::::::::: PlotPhasePortraitLinear.R ::::::::::::::::::::")
  setwd(SCRIPTS_DIRECTORY)
  print("Running script Norm_PlotPhasePortraitLinearAll.R")
  source("Norm_PlotPhasePortraitLinearAll.R", local = T)
}, error = function(e) {print("Error with Norm_PlotPhasePortraitLinearAll.R")})

# Run Summary----
tryCatch({
  print(":::::::::::::::::::: Norm_Summarize.R ::::::::::::::::::::")
  setwd(SCRIPTS_DIRECTORY)
  print("Running script Norm_Summarize.R")
  source("Norm_Summarize.R", local = T)
}, error = function(e) {print("Error with Norm_Summarize.R")})

#
# # Plot Phase Portrait Linear----
tryCatch({
  print(":::::::::::::::::::: PlotPhasePortraitLinear.R ::::::::::::::::::::")
  setwd(SCRIPTS_DIRECTORY)
  print("Running script Norm_PlotPhasePortraitLinear.R")
  source("Norm_PlotPhasePortraitLinear.R", local = T)
}, error = function(e) {print("Error with PlotPhasePortraitLinear.R")})

