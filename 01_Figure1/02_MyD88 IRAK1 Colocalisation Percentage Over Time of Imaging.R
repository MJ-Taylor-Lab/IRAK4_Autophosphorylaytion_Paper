Plot_Directory_Save_Path_Sub_Dir <- file.path(Plot_Directory_Save_Path, "02_IRAK1 Colocalisation Number per cell Sweep")
if(!file.exists(Plot_Directory_Save_Path_Sub_Dir)){
  dir.create(Plot_Directory_Save_Path_Sub_Dir)
}


# Table FIlters -----------------------------------------------------------
Table_filtered <- Table %>% 
  filter(
    IMAGE == "20220921 plate002_well9B_5nM_cl082_MyD88_IRAK1KI_notreat_"
  )


# Parameter Set -----------------------------------------------------------
Seconds_per_Frame = 3
Dwell_Frame_List <- seq(3, 3, by = 1) ### How many frames count as colocalisation... add one to the list to know number of frames
Dwell_Frame_List <- Dwell_Frame_List %>%  as.data.table()
colnames(Dwell_Frame_List) <- c("DWELL_FRAMES")

Dwell_Test_List <- Dwell_Frame_List %>% 
  mutate(
    ROW_NUMBER = row_number(),
    ROW_NUMBER = sprintf("%02d", ROW_NUMBER),
    DWELL_TIME = DWELL_FRAMES*Seconds_per_Frame,
    SAVE_NAME = paste0(ROW_NUMBER, "_1 - Colocalisation Count per cell vs Time since cell landing for colocalisation Thresholdgreater than and equal to ", DWELL_FRAMES,  " frames or ", round(DWELL_FRAMES*3), " s.pdf"),
    SAVE_NAME_2 = paste0(ROW_NUMBER, "_2 - Colocalisation Count per cell vs Time since cell landing for colocalisation Threshold greater than and equal to ", DWELL_FRAMES,  " frames or ", round(DWELL_FRAMES*3), " s without statistics.pdf")
  )

rm(Dwell_Frame_List)

# Colocalisation Calculation ----------------------------------------------
for(index in 1:nrow(Dwell_Test_List)){
  Time_Step <- c(as.numeric(Dwell_Test_List$DWELL_FRAMES[index]*3), seq(30, 1080, by = 30)) #SECONDS_PER_FRAME*2
  
  Colocalisation_Fx <- function(Count){
    # Making a List of Tracks that Colocalize
    Cell_Summary_by_Track <- Table_filtered %>%
      filter(
        PROTEIN == "MyD88",
        NORMALIZED_INTENSITY >= 1.0,
        MAX_NORMALIZED_INTENSITY >= 0.75
      ) %>% 
      group_by(
        UNIVERSAL_TRACK_ID
      ) %>% 
      mutate(
        LIFETIME = max(TIME_ADJUSTED) - min(TIME_ADJUSTED)
      ) %>% 
      group_by(
        IMAGE,
        CELL
      ) %>% 
      mutate(
        TRACK_TIME_SINCE_LANDING = TIME_SINCE_LANDING - min(TIME_SINCE_LANDING)
      ) %>% 
      group_by(
        UNIVERSAL_TRACK_ID
      ) %>% 
      arrange(
        FRAME,
        .by_group = TRUE
      ) %>% 
      filter(
        TRACK_TIME_SINCE_LANDING <= Time_Step[Count]
      ) %>% 
      mutate(
        COLOCALIZATION = NORMALIZED_INTENSITY >= 0.75 #threshold at which recruitment is counted
      ) %>% 
      group_by(
        UNIVERSAL_TRACK_ID
      ) %>%
      mutate(
        STREAK = cumsum(!COLOCALIZATION), #new column which creates a grouping variable for continuous frames which are above threshold
        DWELL_TIME = lead(TIME_ADJUSTED, default = last(TIME_ADJUSTED)) -  TIME_ADJUSTED
      ) %>%
      filter(
        COLOCALIZATION == TRUE
      ) %>% 
      group_by(
        COHORT,
        SHORT_LABEL,
        IMAGE,
        UNIVERSAL_TRACK_ID,
        STREAK #group continuous frames above threshold together
      ) %>%
      summarise(
        DWELL_FRAMES = sum(COLOCALIZATION), #number of frames that complementary protein is above threshold in a continuous stretch
        DWELL_TIME = sum(DWELL_TIME),
        MAX_NORMALIZED_INTENSITY = max(COMPLEMENTARY_NORMALIZED_INTENSITY_1),
        LIFETIME = max(LIFETIME)
      )  %>% 
      filter(
        DWELL_FRAMES >= Dwell_Test_List$DWELL_FRAMES[index]-1 #I want to look at transient recruitment events so dwell frame of 2 implies 3 contionous frames is interesting for me, different cutoffs are prbly useful
      ) %>%
      mutate(
        DWELL_TIME_BIN = round(DWELL_TIME)
      ) %>%
      as.data.table()
    
    
    # Table Creation ----------------------------------------------------------
    # Creating Colocalised List
    Colocalisation_Percentage_List <- Cell_Summary_by_Track %>% 
      distinct(
        UNIVERSAL_TRACK_ID
      )
    
    ### On Cell By cell basis we need to calculate how many tracks are colocalised
    Colocalisation_Percentage_byCell <- Table_filtered %>%
      filter(
        PROTEIN == "MyD88",
        NORMALIZED_INTENSITY >= 1.0,
        MAX_NORMALIZED_INTENSITY >= 0.75
      ) %>% 
      group_by(
        IMAGE,
        CELL
      ) %>% 
      distinct(
        UNIVERSAL_TRACK_ID,
        .keep_all = TRUE
      ) %>% 
      mutate(
        COLOCLIZED_SPOT_TEST = UNIVERSAL_TRACK_ID %in% Colocalisation_Percentage_List$UNIVERSAL_TRACK_ID,
        COLOCLIZED_SPOT_TEST = case_when(
          COLOCLIZED_SPOT_TEST == "TRUE" ~ 1,
          COLOCLIZED_SPOT_TEST != "TRUE" ~ 0
        )
      ) %>% 
      arrange(
        UNIVERSAL_TRACK_ID
      ) %>%
      group_by(
        IMAGE,
        COHORT,
        SHORT_LABEL,
        CELL
      ) %>%
      summarise(
        NUMBER_OF_COLOCLIZED_TRACK = sum(COLOCLIZED_SPOT_TEST)
      ) %>% 
      arrange(
        CELL
      ) %>% 
      mutate(
        TIME_STEP = Time_Step[Count]
      ) %>% 
      as.data.table()
    
    Colocalisation_Percentage_byCell <- Colocalisation_Percentage_byCell[,-(1:3)]
    
    return(Colocalisation_Percentage_byCell)
  }
  
  Colocalisation_Table_Split <- lapply(1:length(Time_Step), Colocalisation_Fx)
  Colocalisation_Table_Split <- rbindlist(Colocalisation_Table_Split) %>% as.data.table()
  
  Colocalisation_Table_Split <- Colocalisation_Table_Split %>% 
    mutate(
      COLOR = case_when(
        CELL == 1 ~ "#fde725",
        CELL == 2 ~ "#efe51c",
        CELL == 3 ~ "#dde318",
        CELL == 4 ~ "#cae11f",
        CELL == 5 ~ "#b8de29",
        CELL == 6 ~ "#a5db36",
        CELL == 7 ~ "#93d741",
        CELL == 8 ~ "#81d34d",
        CELL == 9 ~ "#70cf57",
        CELL == 10 ~ "#60ca60",
        CELL == 11 ~ "#52c569",
        CELL == 12 ~ "#44bf70",
        CELL == 13 ~ "#38b977",
        CELL == 14 ~ "#2fb47c",
        CELL == 15 ~ "#27ad81",
        CELL == 16 ~ "#22a785",
        CELL == 17 ~ "#1fa188",
        CELL == 18 ~ "#1f9a8a",
        CELL == 19 ~ "#20938c",
        CELL == 20 ~ "#228d8d",
        CELL == 21 ~ "#24868e",
        CELL == 22 ~ "#27808e",
        CELL == 23 ~ "#29798e",
        CELL == 24 ~ "#2c728e",
        CELL == 25 ~ "#2f6c8e",
        CELL == 26 ~ "#31668e",
        CELL == 27 ~ "#355f8d",
        CELL == 28 ~ "#38588c",
        CELL == 29 ~ "#3c508b",
        CELL == 30 ~ "#3f4889",
        CELL == 31 ~ "#424086",
        CELL == 32 ~ "#453882",
        CELL == 33 ~ "#472f7d",
        CELL == 34 ~ "#482677",
        CELL == 35 ~ "#481d6f",
        CELL == 36 ~ "#481467",
        CELL == 37 ~ "#460a5d",
        TRUE ~ "#440154"
      )
    )
  

  # GGPLOT_Line Graph -------------------------------------------------------
  Plot <- ggplot(
    data = Colocalisation_Table_Split,
    aes(
      x = TIME_STEP,
      y = NUMBER_OF_COLOCLIZED_TRACK
    )
  ) +
    geom_line(
      aes(
        colour = Colocalisation_Table_Split$COLOR,
      ),
      linewidth = 0.3,
    ) +
    labs(
      x = "Time (s)",
      y = "Colocalised Tracks per Cell (a.u.)"
    ) +
    scale_x_continuous(
      limits = c(0,max(Colocalisation_Table_Split$TIME_STEP)+1),
      breaks = seq(0,max(Colocalisation_Table_Split$TIME_STEP)+1, by = 240),
      expand = c(0,0)
    ) +
    scale_y_continuous(
      n.breaks = 6
    ) +
    theme_classic()
  
  #### Plot with Significance and axis labels and facet labels
  Plot +
    theme(
      axis.title = element_text(size = 30),
      axis.text.y = element_text(size = 20),
      axis.text.x = element_text(size = 20)
    ) 
  
  Plot_Save_Path <- file.path(Plot_Directory_Save_Path_Sub_Dir, Dwell_Test_List$SAVE_NAME[index])
  ggsave(
    Plot_Save_Path,
    plot = last_plot(),
    height = 5*2,
    width = 5*4
  )
  
  #### Plot with y axis and axis labels
  Plot +
    theme(
      axis.title = element_text(size = 5),
      axis.text = element_text(size = 5),
      legend.position = "none",
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.background = element_blank(),
      legend.background = element_blank()
    )
  
  Plot_Save_Path <-   file.path(Plot_Directory_Save_Path_Sub_Dir, Dwell_Test_List$SAVE_NAME_2[index])
  ggsave(
    Plot_Save_Path,
    plot = last_plot(),
    height = 38,
    width = 36,
    units = "mm"
  )

  rm(
    Colocalisation_Table_Split,
    Plot_Save_Path,
    Plot
  )
}


# Cleanup -----------------------------------------------------------------
rm(
  Dwell_Test_List,
  Table_filtered,
  index,
  Seconds_per_Frame,
  Time_Step,
  Plot_Directory_Save_Path_Sub_Dir,
  Colocalisation_Fx
)


