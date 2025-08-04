library(pacman)
pacman::p_load(readxl, dplyr, data.table, ggplot2)

# Setting directory paths
Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots_short"

Primary_Data <- "/Users/u_niranjan/Desktop/Cell_Figures/00_New_Figures/Primary_Data.xlsx"

Function_folder <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/03_PlotFromPrimaryData/Functions"

# Opening Data for relevant plot ----------------------------------------
Cell_Summary_by_Cell <- read_excel(Primary_Data, sheet = "Supplementary Figure 4B")

# Creating folder to save data --------------------------------------------
# Creating Folder to save Images
Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "04_SupplementaryFigure 4")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

# Creating Folder to save Images
Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "Supplementary Figure 4B")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

Plot <- ggplot(
  data = Cell_Summary_by_Cell,
  aes(
    x = TIME_STEP,
    y = NUMBER_OF_COLOCLIZED_TRACK
  )
) +
  geom_line(
    aes(
      colour = Cell_Summary_by_Cell$COLOR,
    ),
    linewidth = 0.2,
  ) +
  labs(
    x = "Time (s)",
    y = "# Colocalized puncta"
  ) +
  scale_x_continuous(
    limits = c(0,max(Cell_Summary_by_Cell$TIME_STEP)+1),
    breaks = seq(0,max(Cell_Summary_by_Cell$TIME_STEP)+1, by = 300),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    n.breaks = 3
  ) +
  theme_classic()

#### Plot with y axis and axis labels
Plot +
  theme(
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.background = element_blank(),
    legend.background = element_blank()
  )

Plot_Save_Path <-   file.path(Plot_Directory_Save_Path, "01_Number of Colocalisation Events per cell over time.pdf")
ggsave(
  Plot_Save_Path,
  plot = last_plot(),
  height = 38,
  width = 36,
  units = "mm"
)

# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()