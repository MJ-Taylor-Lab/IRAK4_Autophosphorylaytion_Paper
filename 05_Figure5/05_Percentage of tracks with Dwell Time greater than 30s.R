# Dwell Time Image average
Plot_Directory_Save_Path_Sub_Dir <- file.path(Plot_Directory_Save_Path, "02_Percentage of tracks with Dwell time above 30s")
if(!file.exists(Plot_Directory_Save_Path_Sub_Dir)){
  dir.create(Plot_Directory_Save_Path_Sub_Dir)
}

# Dwell Time List-------------------------------------------------------
Cell_Summary_by_Track_List <-
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
    DWELL_FRAMES >= 10 #I want to look at transient recruitment events so dwell frame of 2 implies 3 contionous frames is interesting for me, different cutoffs are prbly useful
  ) %>%
  mutate(
    DWELL_TIME_BIN = round(DWELL_TIME)
  ) %>% 
  as.data.table()

Cell_Summary_by_Track_List <- Cell_Summary_by_Track_List %>% 
  distinct(
    UNIVERSAL_TRACK_ID
  )


# Cell Summary Percentages ------------------------------------------------
Cell_Summary_by_Track <- Table %>% 
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
    RECRUITING_TRACK_TEST = UNIVERSAL_TRACK_ID %in% Cell_Summary_by_Track_List$UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    DWELL_TIME_SPOT_TEST = UNIVERSAL_TRACK_ID %in% Cell_Summary_by_Track_List$UNIVERSAL_TRACK_ID,
    DWELL_TIME_SPOT_TEST = case_when(
      DWELL_TIME_SPOT_TEST == "TRUE" ~ 100,
      DWELL_TIME_SPOT_TEST != "TRUE" ~ 0
    )
  ) %>% 
  arrange(
    UNIVERSAL_TRACK_ID
  ) %>%
  group_by(
    IMAGE,
    COHORT,
    SHORT_LABEL,
    ORDER_NUMBER,
    CELL
  ) %>%
  summarise(
    MEAN_DWELL_TIME_TEST = mean(DWELL_TIME_SPOT_TEST),                   ### We group all the tracks of a cell together and find the mean percentage of colocalisation per cell
  ) %>% 
  mutate(
    COLOR = case_when(
      SHORT_LABEL == "DMSO" ~ "orange",
      TRUE ~ "lightblue"
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
    ORDER_NUMBER,
    IMAGE
  ) %>% 
  as.data.table()


Cell_Summary_by_Track$SHORT_LABEL <- 
  factor(
    Cell_Summary_by_Track$SHORT_LABEL,
    levels = c("DMSO", "Kinase Inhibitor 500 nM", "Kinase Inhibitor 20 uM")
  )

# Summary by Images -------------------------------------------------------
#### Protien-Image Summary
Cell_Summary_by_Image <- Cell_Summary_by_Track %>% 
  group_by(
    IMAGE,
    COHORT,
    SHORT_LABEL,
    ORDER_NUMBER,
    VIOLIN_SHAPE
  ) %>% 
  summarise(
    DWELL_TIME_MEDIAN = median(MEAN_DWELL_TIME_TEST),
    DWELL_TIME_MEAN = mean(MEAN_DWELL_TIME_TEST)
  ) %>% 
  arrange(
    ORDER_NUMBER
  )

Cell_Summary_by_Image$SHORT_LABEL <- 
  factor(
    Cell_Summary_by_Image$SHORT_LABEL,
    levels = c("DMSO", "Kinase Inhibitor 500 nM", "Kinase Inhibitor 20 uM")
  )


# t-test ------------------------------------------------------------------
### DMSO vs Kinase Inhibitor 20 uM
p_value_Result <- Cell_Summary_by_Image %>% 
  filter(
    SHORT_LABEL != "Kinase Inhibitor 500 nM"
  )

p_value_Result_MyD88 <- wilcox.test(
  data = p_value_Result,
  DWELL_TIME_MEAN ~ SHORT_LABEL
)$p.value

p_value_Result <- signif(p_value_Result_MyD88, digits = 3)

df_p_val_DMSOvsKI20uM <- data.frame(
  group1 = "DMSO",
  group2 = "Kinase Inhibitor 20 uM",
  label = p_value_Result,
  y.position = 27.5
)

### Kinase Inhibitor 20 uM vs Kinase Inhibitor 500 nM
p_value_Result <- Cell_Summary_by_Image %>% 
  filter(
    SHORT_LABEL != "DMSO"
  )

p_value_Result_MyD88 <- wilcox.test(
  data = p_value_Result,
  DWELL_TIME_MEAN ~ SHORT_LABEL
)$p.value

p_value_Result <- signif(p_value_Result_MyD88, digits = 3)

df_p_val_KI20uMvsKI500nM <- data.frame(
  group1 = "Kinase Inhibitor 20 uM",
  group2 = "Kinase Inhibitor 500 nM",
  label = p_value_Result,
  y.position = 25.5
)

### DMSO vs Kinase Inhibitor 500 nM
p_value_Result <- Cell_Summary_by_Image %>% 
  filter(
    SHORT_LABEL != "Kinase Inhibitor 20 uM"
  )

p_value_Result_MyD88 <- wilcox.test(
  data = p_value_Result,
  DWELL_TIME_MEAN ~ SHORT_LABEL
)$p.value

p_value_Result <- signif(p_value_Result_MyD88, digits = 3)

df_p_val_DMSOvsKI500nM <- data.frame(
  group1 = "DMSO",
  group2 = "Kinase Inhibitor 500 nM",
  label = p_value_Result,
  y.position = 25.5
)

rm(
  p_value_Result,
  p_value_Result_MyD88,
  Cell_Summary_by_Track_List
) 

# Cohort Summary ----------------------------------------------------------
Plot_Table_by_Cohort <- Cell_Summary_by_Track %>% 
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
    COHORT,
    ORDER_NUMBER
  ) %>%
  summarise(
    TOTAL_NUMBER_OF_IMAGES = mean(TOTAL_NUMBER_OF_IMAGES),
    STANDARD_DEVIATION_OF_DWELL_TIME_TEST = sd(MEAN_DWELL_TIME_TEST),
    STANDARD_ERROR_OF_DWELL_TIME = STANDARD_DEVIATION_OF_DWELL_TIME_TEST/TOTAL_NUMBER_OF_IMAGES,
    MEAN_DWELL_TIME = mean(MEAN_DWELL_TIME_TEST),
    YMAX_DWELL_TIME_TEST = MEAN_DWELL_TIME + STANDARD_ERROR_OF_DWELL_TIME,
    YMIN_DWELL_TIME_TEST = MEAN_DWELL_TIME - STANDARD_ERROR_OF_DWELL_TIME,
  ) %>% 
  arrange(
    ORDER_NUMBER
  )


# GGPLOT Violin -----------------------------------------------------------
#####Plot
Plot <- ggplot(
  data = Plot_Table_by_Cohort
) +
  geom_violin(
    data = Cell_Summary_by_Track,
    aes(
      y = MEAN_DWELL_TIME_TEST,
      x = SHORT_LABEL,
      fill = SHORT_LABEL
    ),
    scale = "width",
    alpha = 0.5
  ) +
  scale_fill_manual(
    values = c("orange", "lightblue", "lightblue")
  ) +
  # geom_sina(
  #   data = Cell_Summary_by_Track,
  #   aes(
  #     y = MEAN_DWELL_TIME_TEST,
  #     x = SHORT_LABEL,
  #     fill = SHORT_LABEL
  #   ),
  #   scale = "width",
  #   alpha = 0.3,
  #   size = 0.1,
  #   method = "density",
  #   maxwidth = 0.8,
  #   shape = 27
  # ) +
  geom_segment(
    aes(
      x = 0.9, 
      xend = 1.1,
      y = Plot_Table_by_Cohort$MEAN_DWELL_TIME[1],
      yend = Plot_Table_by_Cohort$MEAN_DWELL_TIME[1]
    ),
    size = 1
  ) +
  geom_segment(
    aes(
      x = 1, 
      xend = 1,
      y = Plot_Table_by_Cohort$MEAN_DWELL_TIME[1] - Plot_Table_by_Cohort$STANDARD_ERROR_OF_DWELL_TIME[1],
      yend = Plot_Table_by_Cohort$MEAN_DWELL_TIME[1] + Plot_Table_by_Cohort$STANDARD_ERROR_OF_DWELL_TIME[1]
    ),
    size = 1
  ) +
  geom_segment(
    aes(
      x = 0.9, 
      xend = 1.1,
      y = Plot_Table_by_Cohort$MEAN_DWELL_TIME[1] - Plot_Table_by_Cohort$STANDARD_ERROR_OF_DWELL_TIME[1],
      yend = Plot_Table_by_Cohort$MEAN_DWELL_TIME[1] - Plot_Table_by_Cohort$STANDARD_ERROR_OF_DWELL_TIME[1]
    ),
    size = 1
  ) +
  geom_segment(
    aes(
      x = 0.9, 
      xend = 1.1,
      y = Plot_Table_by_Cohort$MEAN_DWELL_TIME[1] + Plot_Table_by_Cohort$STANDARD_ERROR_OF_DWELL_TIME[1],
      yend = Plot_Table_by_Cohort$MEAN_DWELL_TIME[1] + Plot_Table_by_Cohort$STANDARD_ERROR_OF_DWELL_TIME[1]
    ),
    size = 1
  ) +
  geom_segment(
    aes(
      x = 1.9, 
      xend = 2.1,
      y = Plot_Table_by_Cohort$MEAN_DWELL_TIME[2],
      yend = Plot_Table_by_Cohort$MEAN_DWELL_TIME[2]
    ),
    size = 1
  ) +
  geom_segment(
    aes(
      x = 2, 
      xend = 2,
      y = Plot_Table_by_Cohort$MEAN_DWELL_TIME[2] - Plot_Table_by_Cohort$STANDARD_ERROR_OF_DWELL_TIME[2],
      yend = Plot_Table_by_Cohort$MEAN_DWELL_TIME[2] + Plot_Table_by_Cohort$STANDARD_ERROR_OF_DWELL_TIME[2]
    ),
    size = 1
  ) +
  geom_segment(
    aes(
      x = 1.9, 
      xend = 2.1,
      y = Plot_Table_by_Cohort$MEAN_DWELL_TIME[2] - Plot_Table_by_Cohort$STANDARD_ERROR_OF_DWELL_TIME[2],
      yend = Plot_Table_by_Cohort$MEAN_DWELL_TIME[2] - Plot_Table_by_Cohort$STANDARD_ERROR_OF_DWELL_TIME[2]
    ),
    size = 1
  ) +
  geom_segment(
    aes(
      x = 1.9, 
      xend = 2.1,
      y = Plot_Table_by_Cohort$MEAN_DWELL_TIME[2] + Plot_Table_by_Cohort$STANDARD_ERROR_OF_DWELL_TIME[2],
      yend = Plot_Table_by_Cohort$MEAN_DWELL_TIME[2] + Plot_Table_by_Cohort$STANDARD_ERROR_OF_DWELL_TIME[2]
    ),
    size = 1
  ) +
  geom_segment(
    aes(
      x = 2.9, 
      xend = 3.1,
      y = Plot_Table_by_Cohort$MEAN_DWELL_TIME[3],
      yend = Plot_Table_by_Cohort$MEAN_DWELL_TIME[3]
    ),
    size = 1
  ) +
  geom_segment(
    aes(
      x = 3, 
      xend = 3,
      y = Plot_Table_by_Cohort$MEAN_DWELL_TIME[3] - Plot_Table_by_Cohort$STANDARD_ERROR_OF_DWELL_TIME[3],
      yend = Plot_Table_by_Cohort$MEAN_DWELL_TIME[3] + Plot_Table_by_Cohort$STANDARD_ERROR_OF_DWELL_TIME[3]
    ),
    size = 1
  ) +
  geom_segment(
    aes(
      x = 2.9, 
      xend = 3.1,
      y = Plot_Table_by_Cohort$MEAN_DWELL_TIME[3] + Plot_Table_by_Cohort$STANDARD_ERROR_OF_DWELL_TIME[3],
      yend = Plot_Table_by_Cohort$MEAN_DWELL_TIME[3] + Plot_Table_by_Cohort$STANDARD_ERROR_OF_DWELL_TIME[3]
    ),
    size = 1
  ) +
  geom_segment(
    aes(
      x = 2.9, 
      xend = 3.1,
      y = Plot_Table_by_Cohort$MEAN_DWELL_TIME[3] - Plot_Table_by_Cohort$STANDARD_ERROR_OF_DWELL_TIME[3],
      yend = Plot_Table_by_Cohort$MEAN_DWELL_TIME[3] - Plot_Table_by_Cohort$STANDARD_ERROR_OF_DWELL_TIME[3]
    ),
    size = 1
  ) +
  scale_y_continuous(
    limits = c(0,31),#c(0,100),
    breaks = seq(0,31, by = 5)#seq(0,100, by = 10)
  ) +
  labs(
    x = "Treatment",
    y = "Colocalised Tracks per Cell (%)"
  ) +
  theme_classic()

##### Plot with p-value and axis
Plot +
  geom_jitter(
    data =  Cell_Summary_by_Image,
    aes(
      x = SHORT_LABEL,
      y = DWELL_TIME_MEAN
    ),
    shape = Cell_Summary_by_Image$VIOLIN_SHAPE,
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

Plot_Save_Path_1 <- "01_Percentage of track over 30s Dwell Time with statistics.pdf"
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
    data =  Cell_Summary_by_Image,
    aes(
      x = SHORT_LABEL,
      y = DWELL_TIME_MEAN
    ),
    shape = Cell_Summary_by_Image$VIOLIN_SHAPE,
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

Plot_Save_Path_1 <-  "02_Percentage of track over 30s Dwell Time without statistics.pdf"
Plot_Save_Path <- file.path(Plot_Directory_Save_Path_Sub_Dir, Plot_Save_Path_1)
ggsave(
  Plot_Save_Path,
  plot = last_plot(),
  height = 45,
  width = 35,
  units = "mm"
)


# Cleanup -----------------------------------------------------------------
rm(
  Cell_Summary_by_Image,
  Cell_Summary_by_Track,
  Plot_Table_by_Cohort,
  df_p_val_DMSOvsKI20uM,
  df_p_val_KI20uMvsKI500nM,
  df_p_val_DMSOvsKI500nM,
  Plot,
  Plot_Path,
  Plot_Save_Path_1,
  Plot_Save_Path,
  Plot_Directory_Save_Path_Sub_Dir
)

