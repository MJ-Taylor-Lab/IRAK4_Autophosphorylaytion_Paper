Plot_Directory_Save_Path_Sub_Dir <- file.path(Plot_Directory_Save_Path, "01_Colocalisation Violin Sweep")
if(!file.exists(Plot_Directory_Save_Path_Sub_Dir)){
  dir.create(Plot_Directory_Save_Path_Sub_Dir)
}


# Colocalised Track Information--------------------------------------------------------
Cell_Summary_by_Track<-
  Table %>%
  filter(
    PROTEIN == "MyD88",
    MAX_NORMALIZED_INTENSITY >= 1.5,
    NORMALIZED_INTENSITY >= 0.75
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  arrange(
    FRAME,
    .by_group = TRUE
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    LIFETIME = max(TIME_ADJUSTED) - min(TIME_ADJUSTED)
  ) %>% 
  mutate(
    COLOCALIZATION = COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 1.5 #threshold at which recruitment is counted
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
    ORDER_NUMBER,
    IMAGE,
    UNIVERSAL_TRACK_ID,
    STREAK #group continuous frames above threshold together
  ) %>%
  summarise(
    DWELL_FRAMES = sum(COLOCALIZATION), #number of frames that complementary protein is above threshold in a continuous stretch
    DWELL_TIME = sum(DWELL_TIME),
    MAX_NORMALIZED_INTENSITY = max(NORMALIZED_INTENSITY),
    LIFETIME = max(LIFETIME)
  ) %>%
  filter(
    DWELL_FRAMES >= 2 #I want to look at transient recruitment events so dwell frame of 2 implies 3 contionous frames is interesting for me, different cutoffs are prbly useful
  ) %>%
  mutate(
    DWELL_TIME_BIN = round(DWELL_TIME)
  ) %>% 
  group_by(
    IMAGE,
    COHORT,
    SHORT_LABEL
  ) %>% 
  mutate(
    DATE = strsplit(IMAGE, " ", fixed = T)[[1]][1],
    DATE = as.integer(DATE)
  ) %>% 
  arrange(
    DATE,
    COHORT,
    IMAGE
  ) %>% 
  transform(
    ID_DATE = as.numeric(factor(DATE))
  ) %>% 
  group_by(
    DATE
  ) %>% 
  mutate(
    COLOR = case_when(
      ID_DATE == 1 ~ "#7CAE00",
      ID_DATE == 2 ~ "#CD9600",
      ID_DATE == 3 ~ "#F8766D",
      ID_DATE == 4 ~ "#00BE67",
      ID_DATE == 5 ~ "#00BFC4",
      ID_DATE == 6 ~ "#00A9FF",
      ID_DATE == 7 ~ "#C77CFF",
      ID_DATE == 8 ~ "#FF61CC"
    )
  ) %>% 
  ungroup() %>% 
  arrange(
    ORDER_NUMBER,
    UNIVERSAL_TRACK_ID
  ) %>% 
  as.data.table()


Cell_Summary_by_Track$SHORT_LABEL <- 
  factor(
    Cell_Summary_by_Track$SHORT_LABEL,
    levels = c("DMSO", "Kinase Inhibitor 500 nM", "Kinase Inhibitor 20 uM")
  )


# Dwell Cycle -------------------------------------------------------------
Dwell_Test_List <- c(1:30) #Number of Dwell Frames
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
  # Creating Colocalised List
  Colocalisation_Percentage_List <- Cell_Summary_by_Track %>% 
    filter(
      DWELL_FRAMES >= Dwell_Test_List$DWELL_FRAMES[Count] #If we want to look at 3 contionous frames we need Dwell Frames = 2 as counting starts from 0
    ) %>%
    distinct(
      UNIVERSAL_TRACK_ID
    )
  
  ### Finding cell means
  Colocalisation_Percentage_byCell <- 
    Table %>% 
    filter(
      PROTEIN == "MyD88",
      MAX_NORMALIZED_INTENSITY >= 1.0,
      NORMALIZED_INTENSITY >= 0.75
    ) %>% 
    distinct(
      UNIVERSAL_TRACK_ID,
      .keep_all = TRUE
    ) %>% 
    mutate(
      COLOCLIZED_SPOT_TEST = UNIVERSAL_TRACK_ID %in% Colocalisation_Percentage_List$UNIVERSAL_TRACK_ID,
      COLOCLIZED_SPOT_TEST = case_when(
        COLOCLIZED_SPOT_TEST == "TRUE" ~ 100,
        COLOCLIZED_SPOT_TEST != "TRUE" ~ 0
      )
    ) %>% 
    arrange(
      UNIVERSAL_TRACK_ID,
      ORDER_NUMBER
    ) %>%
    group_by(
      IMAGE,
      COHORT,
      SHORT_LABEL,
      ORDER_NUMBER,
      CELL
    ) %>%
    summarise(
      MEAN_COLOCLIZED_SPOT_TEST = mean(COLOCLIZED_SPOT_TEST),                   ### We group all the tracks of a cell together and find the mean percentage of colocalisation per cell
    ) %>% 
    group_by(
      IMAGE,
      COHORT,
      SHORT_LABEL
    ) %>% 
    mutate(
      DATE = strsplit(IMAGE, " ", fixed = T)[[1]][1],
      DATE = as.integer(DATE)
    ) %>% 
    arrange(
      DATE,
      COHORT,
      IMAGE
    ) %>% 
    transform(
      ID_DATE = as.numeric(factor(DATE))
    ) %>% 
    group_by(
      DATE
    ) %>% 
    mutate(
      COLOR = case_when(
        ID_DATE == 1 ~ "#7CAE00",
        ID_DATE == 2 ~ "#CD9600",
        ID_DATE == 3 ~ "#F8766D",
        ID_DATE == 4 ~ "#00BE67",
        ID_DATE == 5 ~ "#00BFC4",
        ID_DATE == 6 ~ "#00A9FF",
        ID_DATE == 7 ~ "#C77CFF",
        ID_DATE == 8 ~ "#FF61CC"
      )
    ) %>% 
    ungroup() %>% 
    mutate(
      VIOLIN_COLOR = case_when(
        SHORT_LABEL == "DMSO" ~ "orange",
        SHORT_LABEL == "Kinase Inhibitor 20 uM" ~ "lightblue",
        SHORT_LABEL == "Kinase Inhibitor 500 nM" ~ "lightblue"
      ),
      VIOLIN_SHAPE = case_when(
        SHORT_LABEL == "DMSO" ~ 21,
        SHORT_LABEL == "Kinase Inhibitor 20 uM" ~ 22,
        SHORT_LABEL == "Kinase Inhibitor 500 nM" ~ 24
      )
    ) %>% 
    arrange(
      SHORT_LABEL
    ) %>% 
    as.data.table()
  
  Colocalisation_Percentage_byCell$SHORT_LABEL <- 
    factor(
      Colocalisation_Percentage_byCell$SHORT_LABEL,
      levels = c("DMSO", "Kinase Inhibitor 500 nM", "Kinase Inhibitor 20 uM")
    )
  
  Colocalisation_Percentage_byCell$VIOLIN_COLOR <- 
    factor(
      Colocalisation_Percentage_byCell$VIOLIN_COLOR,
      levels = c("orange", "lightblue")
    )
  
  # Colocalisation by Image
  Colocalisation_Percentage_byImage <- Colocalisation_Percentage_byCell %>% 
    group_by(
      IMAGE,
      COHORT,
      ORDER_NUMBER,
      SHORT_LABEL,
      VIOLIN_SHAPE,
      COLOR
    ) %>% 
    summarise(
      MEAN_COLOCLIZED_SPOT_TEST = mean(MEAN_COLOCLIZED_SPOT_TEST)
    ) %>% 
    arrange(
      ORDER_NUMBER
    ) %>% 
    as.data.table()
  
    # t-test ------------------------------------------------------------------
  ### DMSO vs Kinase Inhibitor 20 uM
  p_value_Result <- Colocalisation_Percentage_byImage %>% 
    filter(
      SHORT_LABEL != "Kinase Inhibitor 500 nM"
    )
  
  p_value_Result_MyD88 <- wilcox.test(
    data = p_value_Result,
    MEAN_COLOCLIZED_SPOT_TEST ~ SHORT_LABEL
  )$p.value
  
  p_value_Result <- signif(p_value_Result_MyD88, digits = 3)
  
  df_p_val_DMSOvsKI20uM <- data.frame(
    group1 = "DMSO",
    group2 = "Kinase Inhibitor 20 uM",
    label = p_value_Result,
    y.position = 23.5
  )
  
  ### Kinase Inhibitor 20 uM vs Kinase Inhibitor 500 nM
  p_value_Result <- Colocalisation_Percentage_byImage %>% 
    filter(
      SHORT_LABEL != "DMSO"
    )
  
  p_value_Result_MyD88 <- wilcox.test(
    data = p_value_Result,
    MEAN_COLOCLIZED_SPOT_TEST ~ SHORT_LABEL
  )$p.value
  
  p_value_Result <- signif(p_value_Result_MyD88, digits = 3)
  
  df_p_val_KI20uMvsKI500nM <- data.frame(
    group1 = "Kinase Inhibitor 20 uM",
    group2 = "Kinase Inhibitor 500 nM",
    label = p_value_Result,
    y.position = 20.5
  )
  
  ### DMSO vs Kinase Inhibitor 500 nM
  p_value_Result <- Colocalisation_Percentage_byImage %>% 
    filter(
      SHORT_LABEL != "Kinase Inhibitor 20 uM"
    )
  
  p_value_Result_MyD88 <- wilcox.test(
    data = p_value_Result,
    MEAN_COLOCLIZED_SPOT_TEST ~ SHORT_LABEL
  )$p.value
  
  p_value_Result <- signif(p_value_Result_MyD88, digits = 3)
  
  df_p_val_DMSOvsKI500nM <- data.frame(
    group1 = "DMSO",
    group2 = "Kinase Inhibitor 500 nM",
    label = p_value_Result,
    y.position = 20.5
  )
  
  rm(
    p_value_Result,
    p_value_Result_MyD88,
    Colocalisation_Percentage_List
  ) 
  

  # Cohort Means ------------------------------------------------------------
  ### Finding Cohort means
  Colocalisation_Percentage_byCOHORT <- Colocalisation_Percentage_byCell %>% 
    arrange(
      SHORT_LABEL,
      IMAGE
    ) %>% 
    transform(
      ID_COHORT = as.numeric(factor(SHORT_LABEL))
    ) %>% 
    group_by(
      ID_COHORT
    ) %>% 
    mutate(
      ID_IMAGE = as.numeric(factor(IMAGE)),
      TOTAL_NUMBER_OF_IMAGES = max(ID_IMAGE)
    ) %>% 
    ungroup() %>% 
    group_by(
      SHORT_LABEL,
      COHORT
    ) %>%
    summarise(
      TOTAL_NUMBER_OF_IMAGES = mean(TOTAL_NUMBER_OF_IMAGES),
      STANDARD_DEVIATION_OF_MEAN_COLOCLIZED_SPOT_TEST = sd(MEAN_COLOCLIZED_SPOT_TEST),
      STANDARD_ERROR_Of_MEAN_SPOT_TEST = STANDARD_DEVIATION_OF_MEAN_COLOCLIZED_SPOT_TEST/TOTAL_NUMBER_OF_IMAGES,
      MEAN_COLOCLIZED_SPOT_TEST = mean(MEAN_COLOCLIZED_SPOT_TEST),
      YMAX_COLOCLIZED_SPOT_TEST = MEAN_COLOCLIZED_SPOT_TEST + STANDARD_ERROR_Of_MEAN_SPOT_TEST,
      YMIN_COLOCLIZED_SPOT_TEST = MEAN_COLOCLIZED_SPOT_TEST - STANDARD_ERROR_Of_MEAN_SPOT_TEST,
    )
  

  # GGPLOT VIOLIN -----------------------------------------------------------
  #####Plot
  Plot <- ggplot(
    data = Colocalisation_Percentage_byCOHORT
  ) +
    geom_violin(
      data = Colocalisation_Percentage_byCell,
      aes(
        y = MEAN_COLOCLIZED_SPOT_TEST,
        x = SHORT_LABEL,
        fill = SHORT_LABEL
      ),
      scale = "width",
      alpha = 0.5
    ) +
    scale_fill_manual(
      values=c("orange", "lightblue", "lightblue")
    ) +
    geom_segment(
      aes(
        x = 0.9, 
        xend = 1.1,
        y = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[1],
        yend = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[1]
      ),
      size = 0.5
    ) +
    geom_segment(
      aes(
        x = 1, 
        xend = 1,
        y = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[1] - Colocalisation_Percentage_byCOHORT$STANDARD_ERROR_Of_MEAN_SPOT_TEST[1],
        yend = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[1] + Colocalisation_Percentage_byCOHORT$STANDARD_ERROR_Of_MEAN_SPOT_TEST[1]
      ),
      size = 0.5
    ) +
    geom_segment(
      aes(
        x = 0.9, 
        xend = 1.1,
        y = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[1] - Colocalisation_Percentage_byCOHORT$STANDARD_ERROR_Of_MEAN_SPOT_TEST[1],
        yend = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[1] - Colocalisation_Percentage_byCOHORT$STANDARD_ERROR_Of_MEAN_SPOT_TEST[1]
      ),
      size = 0.5
    ) +
    geom_segment(
      aes(
        x = 0.9, 
        xend = 1.1,
        y = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[1] + Colocalisation_Percentage_byCOHORT$STANDARD_ERROR_Of_MEAN_SPOT_TEST[1],
        yend = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[1] + Colocalisation_Percentage_byCOHORT$STANDARD_ERROR_Of_MEAN_SPOT_TEST[1]
      ),
      size = 0.5
    ) +
    geom_segment(
      aes(
        x = 1.9, 
        xend = 2.1,
        y = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[2],
        yend = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[2]
      ),
      size = 0.5
    ) +
    geom_segment(
      aes(
        x = 2, 
        xend = 2,
        y = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[2] - Colocalisation_Percentage_byCOHORT$STANDARD_ERROR_Of_MEAN_SPOT_TEST[2],
        yend = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[2] + Colocalisation_Percentage_byCOHORT$STANDARD_ERROR_Of_MEAN_SPOT_TEST[2]
      ),
      size = 0.5
    ) +
    geom_segment(
      aes(
        x = 1.9, 
        xend = 2.1,
        y = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[2] - Colocalisation_Percentage_byCOHORT$STANDARD_ERROR_Of_MEAN_SPOT_TEST[2],
        yend = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[2] - Colocalisation_Percentage_byCOHORT$STANDARD_ERROR_Of_MEAN_SPOT_TEST[2]
      ),
      size = 0.5
    ) +
    geom_segment(
      aes(
        x = 1.9, 
        xend = 2.1,
        y = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[2] + Colocalisation_Percentage_byCOHORT$STANDARD_ERROR_Of_MEAN_SPOT_TEST[2],
        yend = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[2] + Colocalisation_Percentage_byCOHORT$STANDARD_ERROR_Of_MEAN_SPOT_TEST[2]
      ),
      size = 0.5
    ) +
    geom_segment(
      aes(
        x = 2.9, 
        xend = 3.1,
        y = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[3],
        yend = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[3]
      ),
      size = 0.5
    ) +
    geom_segment(
      aes(
        x = 3, 
        xend = 3,
        y = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[3] - Colocalisation_Percentage_byCOHORT$STANDARD_ERROR_Of_MEAN_SPOT_TEST[3],
        yend = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[3] + Colocalisation_Percentage_byCOHORT$STANDARD_ERROR_Of_MEAN_SPOT_TEST[3]
      ),
      size = 0.5
    ) +
    geom_segment(
      aes(
        x = 2.9, 
        xend = 3.1,
        y = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[3] + Colocalisation_Percentage_byCOHORT$STANDARD_ERROR_Of_MEAN_SPOT_TEST[3],
        yend = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[3] + Colocalisation_Percentage_byCOHORT$STANDARD_ERROR_Of_MEAN_SPOT_TEST[3]
      ),
      size = 0.5
    ) +
    geom_segment(
      aes(
        x = 2.9, 
        xend = 3.1,
        y = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[3] - Colocalisation_Percentage_byCOHORT$STANDARD_ERROR_Of_MEAN_SPOT_TEST[3],
        yend = Colocalisation_Percentage_byCOHORT$MEAN_COLOCLIZED_SPOT_TEST[3] - Colocalisation_Percentage_byCOHORT$STANDARD_ERROR_Of_MEAN_SPOT_TEST[3]
      ),
      size = 0.5
    ) +
    scale_y_continuous(
      limits = c(0,26),#c(0,100),
      breaks = seq(0,26, by = 5)#seq(0,100, by = 10)
    ) +
    labs(
      x = "Treatment",
      y = "Colocalised Tracks per Cell (%)"
    ) +
    theme_classic()
  
  ##### Plot with p-value and axis
  Plot +
    geom_jitter(
      data =  Colocalisation_Percentage_byImage,
      aes(
        x = SHORT_LABEL,
        y = MEAN_COLOCLIZED_SPOT_TEST
      ),
      shape = Colocalisation_Percentage_byImage$VIOLIN_SHAPE,
      color = "red", #Colocalisation_Percentage_byImage$COLOR,
      size = 6,
      position=position_jitter(0.06)
    ) +
    add_pvalue(
      df_p_val_DMSOvsKI20uM,
      xmin = "group1",
      xmax = "group2",
      label = "p = {label}",
      y.position = "y.position",
      label.size = 9
    ) +
    add_pvalue(
      df_p_val_KI20uMvsKI500nM,
      xmin = "group1",
      xmax = "group2",
      label = "p = {label}",
      y.position = "y.position",
      label.size = 9
    ) +
    add_pvalue(
      df_p_val_DMSOvsKI500nM,
      xmin = "group1",
      xmax = "group2",
      label = "p = {label}",
      y.position = "y.position",
      label.size = 9
    ) +
    theme(
      axis.title = element_text(size = 30),
      axis.text.y = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      legend.position = "none"
    )
  
  Plot_Save_Path_1 <- paste0(Dwell_Test_List$SAVE_NAME[Count],".pdf")
  Plot_Save_Path <- file.path(Plot_Directory_Save_Path_Sub_Dir, Plot_Save_Path_1)
  ggsave(
    Plot_Save_Path,
    plot = last_plot(),
    height = 5*2,
    width = 5*3
  )
  
  ### Plot without stats and axis labels
  Plot +
    geom_jitter(
      data =  Colocalisation_Percentage_byImage,
      aes(
        x = SHORT_LABEL,
        y = MEAN_COLOCLIZED_SPOT_TEST
      ),
      shape = Colocalisation_Percentage_byImage$VIOLIN_SHAPE,
      fill = "#888888",
      size = 1.2,
      position=position_jitter(0.08),
      alpha = 0.8
    ) +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_text(size = 5),
      axis.text.x = element_text(size = 5),
      legend.position = "none",
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
  
  Plot_Save_Path_1 <- paste0(Dwell_Test_List$SAVE_NAME_2[Count],".pdf")
  Plot_Save_Path <- file.path(Plot_Directory_Save_Path_Sub_Dir, Plot_Save_Path_1)
  ggsave(
    Plot_Save_Path,
    plot = last_plot(),
    height = 45,
    width = 35,
    units = "mm"
  )
  
  rm(
    Colocalisation_Percentage_byCell,
    Colocalisation_Percentage_byImage,
    Colocalisation_Percentage_byCOHORT,
    df_p_val_DMSOvsKI20uM,
    df_p_val_KI20uMvsKI500nM,
    df_p_val_DMSOvsKI500nM,
    Plot,
    Plot_Path,
    Plot_Save_Path_1,
    Plot_Save_Path
  )
}

lapply(1:nrow(Dwell_Test_List), Colocalisation_Fx)


rm(
  Plot_Directory_Save_Path_Sub_Dir,
  Cell_Summary_by_Track,
  Dwell_Test_List,
  Colocalisation_Fx
)
