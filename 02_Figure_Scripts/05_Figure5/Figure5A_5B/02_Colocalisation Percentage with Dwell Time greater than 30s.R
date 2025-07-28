# Colocalised Track Information--------------------------------------------------------
# Get A list of Universal tracks with Dwell Times
Function_script <- file.path(Function_folder, "CellSummarybyImage_ForColocalisationPercentage.R")
source(Function_script)
Cell_Summary_by_Track <- Cell_Summary_by_Track_for_ColocalisationPercentage(Table)


# Dwell Cycle -------------------------------------------------------------
# We can create a list of dwell time/ frame filters we can choose from but for publication we only need 10 frames or 30s
Dwell_Test_List <- c(10) #Number of Dwell Frames
Dwell_Test_List <- Dwell_Test_List %>%  as.data.table()
colnames(Dwell_Test_List) <- c("DWELL_FRAMES")
Dwell_Test_List <- Dwell_Test_List %>% 
  mutate(
    ROW_NUMBER = row_number(),
    ROW_NUMBER = sprintf("%02d", ROW_NUMBER),
    SAVE_NAME = paste0(ROW_NUMBER, "_1 - Colocalisation Violin Dwell Frames greater than equal to ", DWELL_FRAMES,  " frames or Dwell time greater than or equal to ", round(DWELL_FRAMES*3), " s"),
    SAVE_NAME_2 = paste0(ROW_NUMBER, "_2 - Colocalisation Violin Dwell Frames greater than equal to ", DWELL_FRAMES,  " frames or Dwell time greater than or equal to ", round(DWELL_FRAMES*3), " s without statistics and axis label")
  )


# Colocalisation Function -------------------------------------------------
Colocalisation_Fx <- function(Count){

  # Table Creation ----------------------------------------------------------
  ### Finding cell means by applying Dwell Time Filter

  Colocalisation_Percentage_byCell <- Colocalisation_Percentage_byCell(
    Table, 
    Cell_Summary_by_Track, 
    Count
  )
  
  # Assigning the right table to the right variable from function
  Colocalisation_Percentage_byImage <- Colocalisation_Percentage_byCell$Colocalisation_Percentage_byImage
  Colocalisation_Percentage_byCOHORT <- Colocalisation_Percentage_byCell$Colocalisation_Percentage_byCOHORT
  Colocalisation_Percentage_byCell <- Colocalisation_Percentage_byCell$Colocalisation_Percentage_byCell
  
  
  Colocalisation_Percentage_byCell <- Colocalisation_Percentage_byCell %>% 
  mutate(
    VIOLIN_COLOR = case_when(
      SHORT_LABEL == "D329A/\nT342E/\nT345E/\nS346E" ~ "#F99108",  # Assign colors based on SHORT_LABEL
      SHORT_LABEL == "D329A" ~ "#F99108",
      SHORT_LABEL == "WT" ~ "#33AFCC"
    )
  )
  
  Colocalisation_Percentage_byCell$SHORT_LABEL <- 
    factor(
      Colocalisation_Percentage_byCell$SHORT_LABEL,
      levels = c("D329A/\nT342E/\nT345E/\nS346E", "D329A", "WT")
    )
  
  Colocalisation_Percentage_byImage <- Colocalisation_Percentage_byImage %>% 
    mutate(
      JITTER_COLOR = case_when(
        SHORT_LABEL == "D329A/\nT342E/\nT345E/\nS346E" ~ "#F99108",  # Assign colors based on SHORT_LABEL
        SHORT_LABEL == "D329A" ~ "#F99108",
        SHORT_LABEL == "WT" ~ "#33AFCC"
      )
    )
  
  #Arrange the COhort Table in order to plot mean values correctly
  Colocalisation_Percentage_byCOHORT <- Colocalisation_Percentage_byCOHORT %>% 
    mutate(
      ORDER = case_when(
        SHORT_LABEL == "D329A/\nT342E/\nT345E/\nS346E" ~ 1,  # Assign value based on y axis position in plot
        SHORT_LABEL == "D329A" ~ 2,
        SHORT_LABEL == "WT" ~ 3
      )
    ) %>% 
    arrange(
      ORDER
    )
  
  
    # t-test ------------------------------------------------------------------
  ### DMSO vs Kinase Inhibitor 20 uM
  p_value_test_confirmation = "Y"
  if(p_value_test_confirmation == "Y"){
    # Original list with 2-5 components
    Combination_Table <- unique(Colocalisation_Percentage_byImage$SHORT_LABEL)
    
    # Convert the list to a vector since combn works with vectors
    Combination_Table <- unlist(Combination_Table)
    Combination_Table <- sort(Combination_Table)
    
    # Generate unique combinations of 2 components
    Combination_Table <- combn(Combination_Table, 2, simplify = FALSE)
    
    Combination_Function <- function(index){
      Combination <- Combination_Table[[index]] 
      Combination <- Combination %>% as.data.table()
      colnames(Combination) <- c("SHORT_LABEL")
      Combination <- Combination %>% 
        mutate(
          SHORT_LABEL = as.character(SHORT_LABEL)
        )
      
      p_value_Result <- Colocalisation_Percentage_byImage %>% 
        rowwise() %>% 
        mutate(
          TEST = SHORT_LABEL %in% Combination$SHORT_LABEL
        ) %>% 
        filter(
          TEST == TRUE
        )
      
      p_value_Result_1 <- wilcox.test(
        data = p_value_Result,
        MEAN_COLOCLIZED_SPOT_TEST ~ SHORT_LABEL
      )$p.value
      p_value_Result_1 <- signif(p_value_Result_1, digits = 3)
      
      Temp <- data.table(
        COHORT1 = Combination$SHORT_LABEL[1],
        COHORT2 = Combination$SHORT_LABEL[2],
        p_value = p_value_Result_1
      )
      
      return(Temp)
    }
    
    p_value_Table <- lapply(1:length(Combination_Table), Combination_Function)
    p_value_Table <- rbindlist(p_value_Table)
    
    rm(Combination_Table)
  }
  
  # Write csv ---------------------------------------------------------------
  Source_Data_Path <- paste0(Dwell_Test_List$SAVE_NAME[Count], " ",  basename(Plot_Script_Directory), "_Colocalisation_Percentage_byCell.csv")
  Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Source_Data_Path)
  write.csv(Colocalisation_Percentage_byCell, Plot_Save_Path)
  
  Source_Data_Path <- paste0(Dwell_Test_List$SAVE_NAME[Count], " ",  basename(Plot_Script_Directory), "_Colocalisation_Percentage_byImage.csv")
  Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Source_Data_Path)
  write.csv(Colocalisation_Percentage_byImage, Plot_Save_Path)
  
  Source_Data_Path <- paste0(Dwell_Test_List$SAVE_NAME[Count], " ",  basename(Plot_Script_Directory), "_Colocalisation_Percentage_byCOHORT.csv")
  Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Source_Data_Path)
  write.csv(Colocalisation_Percentage_byCOHORT, Plot_Save_Path)

  # GGPLOT VIOLIN -----------------------------------------------------------
  # Generate a Violin Baseplot
  Function_script <- file.path(Function_folder, "ColocalisationPercentage_3Cohort.R")
  source(Function_script)
  Plot <- Plot_ColocalisationPercentage_grey(
    Colocalisation_Percentage_byCell, 
    Colocalisation_Percentage_byImage, 
    Colocalisation_Percentage_byCOHORT
  )
  
  
  Plot_pvalue <- Plot_ColocalisationPercentage_pvalues(
    Plot, 
    p_value_Table
  )
  
  # Save the reference plot with P-value as a PDF
  Plot_Save_Path_1 <- paste0(Dwell_Test_List$SAVE_NAME[Count],".pdf")
  Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Plot_Save_Path_1)
  ggsave(
    Plot_Save_Path,
    plot = Plot_pvalue,
    height = 3 * 3,
    width = 5 * 4
  )

  
  
  ### Plot without stats and axis labels
  Plot <- Plot_ColocalisationPercentage_forPublication(Plot)
  
  Plot_Save_Path_1 <- paste0(Dwell_Test_List$SAVE_NAME_2[Count],".pdf")
  Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Plot_Save_Path_1)
  ggsave(
    Plot_Save_Path,
    plot = Plot,
    height = 40,
    width = 32,
    units = "mm"
  )

}

lapply(1:nrow(Dwell_Test_List), Colocalisation_Fx)



# Cleanup -----------------------------------------------------------------
rm(
  Plot_Directory_Save_Path,
  Cell_Summary_by_Track,
  Dwell_Test_List,
  Colocalisation_Fx,
  Colocalisation_Percentage_byCell,
  Cell_Summary_by_Track_for_ColocalisationPercentage,
  Plot_ColocalisationPercentage_pvalues,
  Plot_ColocalisationPercentage_forPublication,
  Plot_ComplexLifetime_forPublication_grey,
  Plot_ComplexLifetime_grey,
  Plot_ComplexLifetime_withMeanSEMandpvalues_grey,
  Function_script
)

gc()
