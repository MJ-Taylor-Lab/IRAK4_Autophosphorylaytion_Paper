# Creating Folder to save Max Normalised Intenisty Plots
Plot_Directory_Save_Path_Sub_Dir <- file.path(Plot_Directory_Save_Path, "02_Max Normalised Intensity IRAK1 Lifetime Greater Than 50s Average")
if(!file.exists(Plot_Directory_Save_Path_Sub_Dir)){
  dir.create(Plot_Directory_Save_Path_Sub_Dir)
}

#### Track Summary
Cell_Summary_by_Track <- Table %>% 
  filter(
    PROTEIN != "MyD88",
    MAX_NORMALIZED_INTENSITY >= 1.0,
    NORMALIZED_INTENSITY >= 0.75,
    IMAGE != "20230608 plate002_well7C_5nM_cl082_MyD88_IRAK1KI_inhibitor500nM_00"
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
    SHORT_LABEL,
    ORDER_NUMBER
  ) %>% 
  summarise(
    MAX_NORMALIZED_INTENSITY = max(NORMALIZED_INTENSITY),
  ) %>% 
  mutate(
    MAX_NORMALIZED_INTENSITY_BIN = round(MAX_NORMALIZED_INTENSITY)
  ) %>% 
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
    UNIVERSAL_TRACK_ID
  ) %>% 
  as.data.table()

Cell_Summary_by_Track$SHORT_LABEL <- 
  factor(
    Cell_Summary_by_Track$SHORT_LABEL,
   levels = c("DMSO", "Kinase Inhibitor 500 nM", "Kinase Inhibitor 20 uM")
  )

#### Image Summary
Cell_Summary_by_Image <- Cell_Summary_by_Track %>% 
  group_by(
    IMAGE,
    COHORT,
    SHORT_LABEL,
    ORDER_NUMBER,
    VIOLIN_SHAPE,
    VIOLIN_COLOR
  ) %>% 
  summarise(
    MAX_NORMALIZED_INTENSITY_MEDIAN = median(MAX_NORMALIZED_INTENSITY),
    MAX_NORMALIZED_INTENSITY_MEAN = mean(MAX_NORMALIZED_INTENSITY)
  ) %>% 
  arrange(
    ORDER_NUMBER,
    SHORT_LABEL
  )

#### Cohort Summary
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
  ungroup() %>% 
  group_by(
    SHORT_LABEL,
    COHORT,
    ORDER_NUMBER
  ) %>%
  summarise(
    TOTAL_NUMBER_OF_IMAGES = mean(TOTAL_NUMBER_OF_IMAGES),
    MAX_NORMALIZED_INTENSITY_MEAN = mean(MAX_NORMALIZED_INTENSITY),
    SEM = sd(MAX_NORMALIZED_INTENSITY)/TOTAL_NUMBER_OF_IMAGES
  ) %>% 
  mutate(
    XMIN_MEAN = MAX_NORMALIZED_INTENSITY_MEAN - SEM,
    XMIN_MEAN = case_when(
      XMIN_MEAN < 0 ~ 0.75,
      TRUE ~ XMIN_MEAN
    ),
    XMAX_MEAN = MAX_NORMALIZED_INTENSITY_MEAN + SEM,
    COLOR  = case_when(
      SHORT_LABEL == "DMSO" ~ "orange",
      TRUE ~ "lightblue"
    )
  ) %>% 
  arrange(
    ORDER_NUMBER,
    SHORT_LABEL
  )

# t-test ------------------------------------------------------------------
### DMSO vs Kinase Inhibitor 20 uM
p_value_Result <- Cell_Summary_by_Image %>% 
  filter(
    SHORT_LABEL != "Kinase Inhibitor 500 nM"
  )

p_value_Result_MyD88 <- wilcox.test(
  data = p_value_Result,
  MAX_NORMALIZED_INTENSITY_MEAN ~ SHORT_LABEL
)$p.value

p_value_Result <- signif(p_value_Result_MyD88, digits = 3)

df_p_val_DMSOvsKI20uM <- data.frame(
  group1 = "DMSO",
  group2 = "Kinase Inhibitor 20 uM",
  label = p_value_Result,
  y.position = 30.5
)

### Kinase Inhibitor 20 uM vs Kinase Inhibitor 500 nM
p_value_Result <- Cell_Summary_by_Image %>% 
  filter(
    SHORT_LABEL != "DMSO"
  )

p_value_Result_MyD88 <- wilcox.test(
  data = p_value_Result,
  MAX_NORMALIZED_INTENSITY_MEAN ~ SHORT_LABEL
)$p.value

p_value_Result <- signif(p_value_Result_MyD88, digits = 3)

df_p_val_KI20uMvsKI500nM <- data.frame(
  group1 = "Kinase Inhibitor 20 uM",
  group2 = "Kinase Inhibitor 500 nM",
  label = p_value_Result,
  y.position = 27.5
)

### DMSO vs Kinase Inhibitor 500 nM
p_value_Result <- Cell_Summary_by_Image %>% 
  filter(
    SHORT_LABEL != "Kinase Inhibitor 20 uM"
  )

p_value_Result_MyD88 <- wilcox.test(
  data = p_value_Result,
  MAX_NORMALIZED_INTENSITY_MEAN ~ SHORT_LABEL
)$p.value

p_value_Result <- signif(p_value_Result_MyD88, digits = 3)

df_p_val_DMSOvsKI500nM <- data.frame(
  group1 = "DMSO",
  group2 = "Kinase Inhibitor 500 nM",
  label = p_value_Result,
  y.position = 27.5
)

rm(
  p_value_Result,
  p_value_Result_MyD88
) 


# ggplot ------------------------------------------------------------------
Plot <- ggplot(
  data = Cell_Summary_by_Cohort
) +
  geom_violin(
    data = Cell_Summary_by_Track,
    aes(
      y = MAX_NORMALIZED_INTENSITY,
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
      x = 0.95, 
      xend = 1.05,
      y = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[1],
      yend = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[1]
    ),
    size = 0.5
  ) +
  geom_segment(
    aes(
      x = 1, 
      xend = 1,
      y = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[1] - Cell_Summary_by_Cohort$SEM[1],
      yend = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[1] + Cell_Summary_by_Cohort$SEM[1]
    ),
    size = 0.5
  ) +
  geom_segment(
    aes(
      x = 0.95, 
      xend = 1.05,
      y = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[1] - Cell_Summary_by_Cohort$SEM[1],
      yend = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[1] - Cell_Summary_by_Cohort$SEM[1]
    ),
    size = 0.5
  ) +
  geom_segment(
    aes(
      x = 0.95, 
      xend = 1.05,
      y = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[1] + Cell_Summary_by_Cohort$SEM[1],
      yend = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[1] + Cell_Summary_by_Cohort$SEM[1]
    ),
    size = 0.5
  ) +
  geom_segment(
    aes(
      x = 1.95, 
      xend = 2.05,
      y = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[2],
      yend = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[2]
    ),
    size = 0.5
  ) +
  geom_segment(
    aes(
      x = 2, 
      xend = 2,
      y = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[2] - Cell_Summary_by_Cohort$SEM[2],
      yend = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[2] + Cell_Summary_by_Cohort$SEM[2]
    ),
    size = 0.5
  ) +
  geom_segment(
    aes(
      x = 1.95, 
      xend = 2.05,
      y = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[2] - Cell_Summary_by_Cohort$SEM[2],
      yend = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[2] - Cell_Summary_by_Cohort$SEM[2]
    ),
    size = 0.5
  ) +
  geom_segment(
    aes(
      x = 1.95, 
      xend = 2.05,
      y = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[2] + Cell_Summary_by_Cohort$SEM[2],
      yend = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[2] + Cell_Summary_by_Cohort$SEM[2]
    ),
    size = 0.5
  ) +
  geom_segment(
    aes(
      x = 2.95, 
      xend = 3.05,
      y = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[3],
      yend = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[3]
    ),
    size = 0.5
  ) +
  geom_segment(
    aes(
      x = 3, 
      xend = 3,
      y = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[3] - Cell_Summary_by_Cohort$SEM[3],
      yend = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[3] + Cell_Summary_by_Cohort$SEM[3]
    ),
    size = 0.5
  ) +
  geom_segment(
    aes(
      x = 2.95, 
      xend = 3.05,
      y = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[3] + Cell_Summary_by_Cohort$SEM[3],
      yend = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[3] + Cell_Summary_by_Cohort$SEM[3]
    ),
    size = 0.5
  ) +
  geom_segment(
    aes(
      x = 2.95, 
      xend = 3.05,
      y = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[3] - Cell_Summary_by_Cohort$SEM[3],
      yend = Cell_Summary_by_Cohort$MAX_NORMALIZED_INTENSITY_MEAN[3] - Cell_Summary_by_Cohort$SEM[3]
    ),
    size = 0.5
  ) +
  scale_y_continuous(
    limits = c(0,36),#c(0,100),
    breaks = seq(0,36, by = 5)#seq(0,100, by = 10)
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
      y = MAX_NORMALIZED_INTENSITY_MEAN
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


Plot_Save_Path_1 <- "01_Max Normalized Intensity IRAK1 Lifetime Greater Than 50s Violin with axis labels and statistics.pdf"
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
      y = MAX_NORMALIZED_INTENSITY_MEAN
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

Plot_Save_Path_1 <- "02_Max Normalized Intensity IRAK1 Lifetime Greater Than 50s Violin without axis labels and statistics.pdf"
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
  Cell_Summary_by_Cohort,
  Cell_Summary_by_Image,
  Cell_Summary_by_Track,
  df_p_val_DMSOvsKI20uM,
  df_p_val_KI20uMvsKI500nM,
  df_p_val_DMSOvsKI500nM,
  Plot,
  Plot_Save_Path_1,
  Plot_Save_Path,
  Plot_Directory_Save_Path_Sub_Dir
)


