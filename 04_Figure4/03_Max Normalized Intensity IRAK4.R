# Creating Folder to save Max Normalised Intenisty Plots
Plot_Directory_Save_Path_Sub_Dir <- file.path(Plot_Directory_Save_Path, "02_Max Normalised Intensity IRAK4 Lifetime Greater Than 50s")
if(!file.exists(Plot_Directory_Save_Path_Sub_Dir)){
  dir.create(Plot_Directory_Save_Path_Sub_Dir)
}

#### Indivual Replicate of Max Normalised Intensity
Cell_Summary_by_Track <- Table %>% 
  filter(
    PROTEIN != "MyD88",
    MAX_NORMALIZED_INTENSITY >= 1.0,
    NORMALIZED_INTENSITY >= 0.75
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    LIFETIME = max(TIME_ADJUSTED) - min (TIME_ADJUSTED)
  ) %>% 
  filter(
    LIFETIME >= 50
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID,
    IMAGE,
    COHORT,
    ORDER_NUMBER,
    SHORT_LABEL,
  ) %>% 
  summarise(
    MAX_NORMALIZED_INTENSITY = max(NORMALIZED_INTENSITY),
  ) %>% 
  mutate(
    MAX_NORMALIZED_INTENSITY_BIN = round(MAX_NORMALIZED_INTENSITY)
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
    levels = c("IRAK4 WT", "IRAK4 KD", "IRAK4 DelKDo")
  )

#### Image Summary
Cell_Summary_by_Image <- Cell_Summary_by_Track %>% 
  group_by(
    IMAGE,
    COHORT,
    SHORT_LABEL,
    ORDER_NUMBER
  ) %>% 
  summarise(
    MAX_NORMALIZED_INTENSITY_MEDIAN = median(MAX_NORMALIZED_INTENSITY),
    MAX_NORMALIZED_INTENSITY_MEAN = mean(MAX_NORMALIZED_INTENSITY)
  ) %>% 
  arrange(
    ORDER_NUMBER,
    SHORT_LABEL
  )

#### Protien-Cohort Summary
Cell_Summary_by_Cohort <- Cell_Summary_by_Track %>% 
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
  group_by(
    COHORT,
    SHORT_LABEL,
    MAX_NORMALIZED_INTENSITY_BIN
  ) %>%
  mutate(
    MAX_NORMALIZED_INTENSITY_BIN_MAX = n()
  ) %>%
  ungroup() %>% 
  group_by(
    SHORT_LABEL,
    COHORT,
    ORDER_NUMBER
  ) %>%
  summarise(
    TOTAL_NUMBER_OF_IMAGES = mean(TOTAL_NUMBER_OF_IMAGES),
    MAX_NORMALIZED_INTENSITY_MEAN = mean(MAX_NORMALIZED_INTENSITY),
    SEM = sd(MAX_NORMALIZED_INTENSITY)/TOTAL_NUMBER_OF_IMAGES,
    
    MAX_NORMALIZED_INTENSITY_BIN_MAX = max(MAX_NORMALIZED_INTENSITY_BIN_MAX),
    Y_POSITION = MAX_NORMALIZED_INTENSITY_BIN_MAX/2,
  ) %>% 
  mutate(
    XMIN_MEAN = MAX_NORMALIZED_INTENSITY_MEAN - SEM,
    XMIN_MEAN = case_when(
      XMIN_MEAN < 0 ~ 0.75,
      TRUE ~ XMIN_MEAN
    ),
    XMAX_MEAN = MAX_NORMALIZED_INTENSITY_MEAN + SEM,
    COLOR  = case_when(
      SHORT_LABEL == "IRAK4 WT" ~ "orange",
      SHORT_LABEL == "IRAK4 KD" ~ "lightblue",
      TRUE ~ "darkblue"
    )
  ) %>% 
  arrange(
    ORDER_NUMBER,
    SHORT_LABEL
  )


# Histogram Max Normalised Intensity --------------------------------------
BINWIDTH = 1               #Histogram Bindwidth
LOWER_LIMIT = 0           #Axis Lower Limit
UPPER_LIMIT = 11       #Axis Upper Limit
AXIS_BREAK_SEQ = 2        #X axis tick marks
FACET_ROW_NUMBERS = 3     #Number of facet rows
X_LABEL = "Max normalized intensity (a.u.)"
Y_LABEL = "Count (a.u.)"

Plot <- ggplot(
  data = Cell_Summary_by_Cohort,
  aes(
    fill = SHORT_LABEL
  )
) +
  geom_histogram(
    data = Cell_Summary_by_Track,
    aes(
      x = MAX_NORMALIZED_INTENSITY,
    ),
    color = "black",
    binwidth = BINWIDTH,
    alpha = 0.75,
    size = 0.1
  ) +
  scale_x_continuous(
    limits =  c(LOWER_LIMIT,UPPER_LIMIT), ## Sets x axis limts
    breaks = seq(LOWER_LIMIT,UPPER_LIMIT, by = AXIS_BREAK_SEQ),
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    expand = c(0, 0.1)
  ) +
  facet_rep_wrap(
    ~SHORT_LABEL,
    repeat.tick.labels = TRUE,
    scales = "free_y",
    nrow = FACET_ROW_NUMBERS
  ) +
  scale_fill_manual(
    values = c("orange", "lightblue", "darkblue")
  ) +
  labs(
    x = X_LABEL,
    y = Y_LABEL
  ) +
  theme_classic()

#### Plot with axis and median
Plot +
  geom_vline(
    data = Cell_Summary_by_Cohort,
    aes(
      xintercept = MAX_NORMALIZED_INTENSITY_MEAN,
    ),
    size = 1
  ) +
  geom_errorbar(
    data = Cell_Summary_by_Cohort,
    aes(
      x = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN,
      y = Cell_Summary_by_Cohort$Y_POSITION,
      xmin = Cell_Summary_by_Cohort$XMIN_MEAN, 
      xmax = Cell_Summary_by_Cohort$XMAX_MEAN
    ), 
    width = 0.2,
    position = position_dodge(0.9),
    linetype = "dashed"
  ) +
  annotate(
    "text", 
    x = 5, 
    y = 200,
    label = paste0(
      "IRAK4 WT Mean +/- SEM = ", signif(Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[1], digits = 4), " +/- ", signif(Cell_Summary_by_Cohort$SEM[1], digits = 4), "\n",
      "IRAK4 KD Mean +/- SEM = ", signif(Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[2], digits = 4), " +/- ", signif(Cell_Summary_by_Cohort$SEM[2], digits = 4), "\n",
      "IRAK4 DelKDo Mean +/- SEM = ", signif(Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[3], digits = 4), " +/- ", signif(Cell_Summary_by_Cohort$SEM[3], digits = 4)
    ),
    color="black",
    size=3
  ) +
  theme(
    legend.position="none",
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text = element_text(size = 20),
    strip.text.x = element_text(size = 15)
  )

Plot_Save_Path_1 <- "01_Max Normalized Intensity IRAK4 Histogram with mean and SEM.pdf"
Plot_Save_Path <- file.path(Plot_Directory_Save_Path_Sub_Dir, Plot_Save_Path_1)
ggsave(
  Plot_Save_Path,
  plot = last_plot(),
  height = 3*3,
  width = 5*4
)

#### Plot without axis
Plot +
  theme(
    legend.position="none",
    axis.title = element_blank(),
    axis.text = element_text(size = 5),
    strip.text.x = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  )

Plot_Save_Path_1 <- "02_Max Normalized Intensity IRAK4 Histogram without mean, SEM and no facet labels.pdf"
Plot_Save_Path <- file.path(Plot_Directory_Save_Path_Sub_Dir, Plot_Save_Path_1)
ggsave(
  Plot_Save_Path,
  plot = last_plot(),
  height = 40,
  width = 30,
  units = "mm"
)


# Cleanup -----------------------------------------------------------------
rm(
  Cell_Summary_by_Cohort,
  Cell_Summary_by_Track,
  Cell_Summary_by_Image,
  Plot,
  BINWIDTH,
  LOWER_LIMIT,
  UPPER_LIMIT,
  AXIS_BREAK_SEQ,
  FACET_ROW_NUMBERS,
  X_LABEL,
  Y_LABEL,
  Plot_Save_Path_1,
  Plot_Save_Path,
  Plot_Directory_Save_Path_Sub_Dir
)


