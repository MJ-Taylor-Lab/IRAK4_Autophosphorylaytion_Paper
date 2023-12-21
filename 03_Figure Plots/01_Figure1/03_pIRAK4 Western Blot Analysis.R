# Table Creation ----------------------------------------------------------
Table_full <- Table %>% 
  mutate(
    IRAK4toActin = IRAK4/actin_IRAK4,
    pIRAK4toActin = pIRAK4/actin_pIRAK4,
    pIRAK4toIRAK4 = pIRAK4toActin/IRAK4toActin
  ) %>% 
  group_by(
    replicate
  ) %>% 
  mutate(
    pIRAK4toIRAK4_Normalized = pIRAK4toIRAK4 - min(pIRAK4toIRAK4),
    pIRAK4toIRAK4_Normalized = pIRAK4toIRAK4_Normalized/ max(pIRAK4toIRAK4_Normalized)
  )

Table_summary <- Table_full %>% 
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
    data = Table_summary,
    aes(
      x = stimulation_time,
      y = pIRAK4toIRAK4_Normalized_mean,
      ymin = pIRAK4toIRAK4_Normalized_min,
      ymax = pIRAK4toIRAK4_Normalized_max
    ),
    linewidth = 0.1
  ) +
  geom_line(
    data = Table_summary,
    aes(
      x = stimulation_time,
      y = pIRAK4toIRAK4_Normalized_mean
    ),
    # linetype = "dashed",
    alpha = 0.5
  ) +
  geom_jitter(
    data = Table_full,
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
    data = Table_summary,
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
    x = "Stimulation Time (log scale, min)",
    y = "pIRAK4/ IRAK4 (Normalized)"
  ) +
  scale_x_continuous(
    trans = "log2",
    breaks = Table_summary$stimulation_time,
    expand = expansion(add = 0.3)
  ) +
  scale_y_continuous(
    breaks = seq(0,1.1, by = 0.2)
  ) +
  theme_classic()

# Plot without axis labels
Plot +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5),
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
  width = 40,
  height = 25,
  units = "mm"
)

# Cleanup -----------------------------------------------------------------
rm(
  Plot,
  Table_full,
  Table_summary,
  Plot_Save_Path_1,
  Plot_Save_Path
)

