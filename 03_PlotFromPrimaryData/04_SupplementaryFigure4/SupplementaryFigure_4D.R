library(pacman)
pacman::p_load(readxl, dplyr, data.table, ggplot2)

# Setting directory paths
Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots_short"

Primary_Data <- "/Users/u_niranjan/Desktop/Cell_Figures/00_New_Figures/Primary_Data.xlsx"

Function_folder <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/03_PlotFromPrimaryData/Functions"

# Opening Data for relevant plot ----------------------------------------
Summary_by_BLot <- read_excel(Primary_Data, sheet = "Supplementary Figure 4D")

# Creating folder to save data --------------------------------------------
# Creating Folder to save Images
Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "04_SupplementaryFigure 4")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

# Creating Folder to save Images
Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "Supplementary Figure 4D")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

Summary_by_Cohort <- Summary_by_BLot %>% 
  group_by(
    stimulation_time
  ) %>% 
  summarise(
    pIRAK4toIRAK4_Normalized_mean = mean(pIRAK4toIRAK4_Normalized),
    pIRAK4toIRAK4_Normalized_sd = sd(pIRAK4toIRAK4_Normalized),
    pIRAK4toIRAK4_Normalized_max = pIRAK4toIRAK4_Normalized_mean + pIRAK4toIRAK4_Normalized_sd,
    pIRAK4toIRAK4_Normalized_min = case_when(
      pIRAK4toIRAK4_Normalized_mean - pIRAK4toIRAK4_Normalized_sd < 0 ~ 0,
      TRUE ~  pIRAK4toIRAK4_Normalized_mean - pIRAK4toIRAK4_Normalized_sd
    )
  ) 


# Plot IRAK4 and pIRAK4 log scale---------------------------------------------------
Plot <- ggplot() +
  geom_errorbar(
    data = Summary_by_Cohort,
    aes(
      x = stimulation_time,
      y = pIRAK4toIRAK4_Normalized_mean,
      ymin = pIRAK4toIRAK4_Normalized_min,
      ymax = pIRAK4toIRAK4_Normalized_max
    ),
    linewidth = 0.1
  ) +
  geom_line(
    data = Summary_by_Cohort,
    aes(
      x = stimulation_time,
      y = pIRAK4toIRAK4_Normalized_mean
    ),
    # linetype = "dashed",
    alpha = 0.5
  ) +
  geom_jitter(
    data = Summary_by_BLot,
    aes(
      x = stimulation_time,
      y = pIRAK4toIRAK4_Normalized
    ),
    shape = 21,
    stroke = 0.1,
    fill = "lightgrey",
    size = 0.5,
    position = position_jitter(0)
  ) +
  geom_jitter(
    data = Summary_by_Cohort,
    aes(
      x = stimulation_time,
      y = pIRAK4toIRAK4_Normalized_mean
    ),
    shape = 21,
    size = 0.3,
    position = position_jitter(0),
    fill = "black"
  ) +
  labs(
    x = "Stimulation Time (log2 scale, min)",
    y = "Signal Intensity"
  ) +
  scale_x_continuous(
    trans = "log2",
    breaks = Summary_by_BLot$stimulation_time,
    expand = expansion(add = 0.3)
  ) +
  scale_y_continuous(
    limits = c(0, 1.1),
    expand = c(0, 0),
    breaks = c(0, 0.5, 1),
  ) +
  theme_classic()

# Plot without axis labels
Plot +
  theme(
    axis.title = element_text(size = 7),
    axis.text.x = element_text(size = 6, hjust = 0.5),
    axis.text.y = element_text(size = 6),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.background = element_blank()
  )

Plot_Save_Path_1 <- "01_pIRAK4byIRAK4 Normalised without axis labels with log2 scale.pdf"
Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Plot_Save_Path_1)
ggsave(
  Plot_Save_Path,
  plot = last_plot(),
  width = 80,
  height = 35,
  units = "mm"
)

# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()