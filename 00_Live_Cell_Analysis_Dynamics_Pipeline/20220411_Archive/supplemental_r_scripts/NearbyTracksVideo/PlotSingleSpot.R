setwd(RETURN_TO_DIRECTORY)

SingleSpot <-
  ggplot() +
  geom_tile(
    data  = plot_img_spot %>% filter(t == FOCUS_TIME),
    aes(
      x = x,
      y = y,
      fill = z
    )
  ) +
  geom_point(
    data = plot_other_spots,
    aes(
      x = x,
      y = y
    ),
    color = "green",
    shape = 22,
    size = 2.5
  ) +
  geom_point(
    aes(
      x = 0,
      y = 0
    ),
    color = "white",
    shape = 22,
    size = 2.5
  ) +
  geom_text(
    aes(
      x = -13.5,
      y = -18,
      label = paste(FOCUS_TIME, "s")
    ),
    size = 4,
    color = "white"
  ) +
  scale_y_reverse(
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    expand = c(0, 0)
  ) +
  scale_fill_viridis(
    option = "inferno",
    na.value = "gray",
    limits = c(0, max(plot_img_spot$z))
  ) +
  facet_wrap(
    ~ t
  ) +
  labs(
    title = "Single Spot",
    x = "Pixels",
    y = "Pixels",
    fill = "Intensity (x GFP)"
  ) +
  coord_fixed() +
  dark_theme_classic() +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
