setwd(RETURN_TO_DIRECTORY)

TracksPlot <-
  ggplot(
    plot_spot,
    aes(
      x = t,
      y = NORMALIZED_INTENSITY
    )
  ) +
  geom_vline(
    xintercept = FOCUS_TIME,
    color = "magenta"
  ) +
  geom_line(
    color = "white"
  ) +
  labs(
    title = "Complete Track Intensities",
    x = "Time (s)",
    y = "Norm. Intensity (a.u.)"
  ) +
  scale_x_continuous(
    limits = c(min(RANGE), max(RANGE)),
  ) +
  scale_y_continuous(
    breaks = pretty_breaks()
  ) +
  dark_theme_classic()
