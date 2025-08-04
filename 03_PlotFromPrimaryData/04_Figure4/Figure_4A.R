library(pacman)
pacman::p_load(readxl, dplyr, data.table, ggplot2)

# Setting directory paths
Plot_Directory_Save_Path <- "/Users/u_niranjan/Desktop/Cell_Figures/Plots_short"

Primary_Data <- "/Users/u_niranjan/Desktop/Cell_Figures/00_New_Figures/Primary_Data.xlsx"

Function_folder <- "/Users/u_niranjan/Desktop/Git Scripts/01_IRAK4 Phosphorylation Paper/IRAK4_Autophosphorylaytion_Paper/03_PlotFromPrimaryData/Functions"

# Opening Data for relevant plot ----------------------------------------
Cell_Summary_by_Cohort <- read_excel(Primary_Data, sheet = "Figure 4A")

# Creating folder to save data --------------------------------------------
# Creating Folder to save Images
Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "04_Figure 4")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}

# Creating Folder to save Images
Plot_Directory_Save_Path <- file.path(Plot_Directory_Save_Path, "Figure_4A")
if(!file.exists(Plot_Directory_Save_Path)){
  dir.create(Plot_Directory_Save_Path)
}


# Convert SHORT_LABEL to a factor with specified levels
Cell_Summary_by_Cohort$Combination <- factor(
  Cell_Summary_by_Cohort$Combination,
  levels = c("MyD88+IRAK4+pIRAK4", "MyD88+IRAK4", "MyD88_only")
  )

# Diet Bar Plot ----------------------------------------------------------------
BarPlot <- ggplot(
  data = Cell_Summary_by_Cohort,
  aes(
    y = Combination,
    x = Proportion,
    fill = Combination
  )
) +
  geom_bar(
    stat = "identity",
    color = "black",
    width = 0.53,
    size = 0.1
  ) +
  scale_fill_manual(
    values = c("#E176F7", "#F0B14A", "#00FFFF")
  ) +
  theme_classic()

BarPlot +
  theme(
    legend.position = "none",
    axis.title = element_blank()
  ) +
  annotate(
    "text",
    x = 15,
    y = 3,
    label = Cell_Summary_by_Cohort$Combination_Name_Full[3]
  ) +
  annotate(
    "text",
    x = 15,
    y = 2,
    label = Cell_Summary_by_Cohort$Combination_Name_Full[2]
  ) +
  annotate(
    "text",
    x = 15,
    y = 1,
    label = Cell_Summary_by_Cohort$Combination_Name_Full[1]
  )

Plot_Save_Path <- file.path(Plot_Directory_Save_Path, "01_Barplot skinny with labels.pdf")
ggsave(
  Plot_Save_Path,
  plot = last_plot(),
  height = 4*4,
  width = 4*4
)


BarPlot +
  scale_x_continuous(
    limits = c(0,50),
    breaks = seq(0,50, by = 20),
    # position="top",
    expand = c(0, 0)
  ) +
  labs(
    x = "% of events observed"
  ) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 7),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(size = 0.1),
    axis.line.y = element_blank(),
    axis.ticks = element_line(size = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(), # Remove the panel background
    plot.background = element_blank(),  # Remove the plot background
    legend.background = element_blank(), # Remove the legend background
  )

Plot_Save_Path <- file.path(Plot_Directory_Save_Path, "02_Barplot skinny without labels.pdf")
ggsave(
  Plot_Save_Path,
  plot = last_plot(),
  height = 35,
  width = 52,
  units = "mm"
)

# Cleanup -----------------------------------------------------------------
rm(list = ls())
gc()