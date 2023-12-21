setwd(RETURN_TO_DIRECTORY)

LimitedTracksPlot <-
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
    title = "Limited Track Intensities",
    x = "Time (s)",
    y = "Norm. Intensity (a.u.)"
  ) +
  scale_y_continuous(
    limits = c(0, 8),
    breaks = pretty_breaks()
  ) +
  scale_x_continuous(
    limits = c(min(RANGE), max(RANGE))
  ) +
  dark_theme_classic()
