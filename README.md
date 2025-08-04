# IRAK4 Phosphorylation - Scripts

This repository contains scripts needed to analyse data for the paper. 
All the way from analysing the '.nd2' files to plots in figures.
Repository contains 4 subfolders. Description of the function of each subfolder is litsed below.
For those who wish to plot the data themselves from the ‘Primary_Data.xlsx’ file please skip directly to '03_PlotFromPrimaryData'

Scripts On Github:
00_Live_Cell_Analysis_Dynamics_Pipeline
    Images acquired in .nd2 format are supplied to the script ‘00_Live_Cell_Analysis_Dynamics_Pipeline’. Running these scripts generates ‘Analysis.csv.gz’ tables that are used for further processing. To install and run this script please refer to the READ.me file within ‘00_Live_Cell_Analysis_Dynamics_Pipeline’.



01_Cohort_Table_Creation
    Usually, images supplied to ‘00_Live_Cell_Analysis_Dynamics_Pipeline’ are based on experimental day. This results in separate output files for each of the days. ‘01_Cohort_Table_Creation’ takes the different ‘Analysis.csv.gz’ files from different output days and combines them together for each cohort for connivence. It also provides an additional column called ‘SHORT_LABEL’ that makes data. Processing and plot labelling easier in future steps.



02_Figure_Scripts
    This folder contains the scripts needed to produce the plots seen in the figures. The subfolders are labelled based on Plot position in the Figures. ‘02_Figure_Scripts’ use the necessary tables from ‘01_Cohort_Table_Creation’ to create plots. Most Figures have a “00_Figurexx_CentralScript.R” script. By providing the Folder in which this script is stored and running multiple figures based on the same cohort can be generated. ‘00_Figurexx_CentralScript.R’ calls on other scripts in the same subfolder which in turn calls scripts in the folder ‘00_Functions’. Some user input is required to operate scripts like path to ‘00_Functions’, Path to subfolder where the script that is to be run is and Path to store the plots generated. Apart from that some scripts require user input to adjust graph axis or options. Currently parameters are set based on plots in figures.



03_PlotFromPrimaryData
    For readers that want to analyse the data and generate their own plots we have generated Scripts in ‘03_PlotFromPrimaryData’. Just download the ‘Primary_Data.xlsx’. And the script of the figure you wish to analyse. Provide the path of ‘Primary_Data.xlsx’, path to the ‘Functions’ file within ‘03_PlotFromPrimaryData’ stored on your computer and a place to store the generated plots. 

