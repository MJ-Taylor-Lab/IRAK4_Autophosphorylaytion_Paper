setwd(RETURN_TO_DIRECTORY)

AllSpots <-
  ggplot() +
  geom_tile(
    data  = plot_all_img_spot,
    aes(
      x = x,
      y = y,
      fill = z
    )
  ) +
  geom_text(
    data  = plot_all_img_spot %>%select(t)%>%distinct(),
    aes(
      x = -7,
      y = -8,
      label = t
    ),
    size = 2.75,
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
    limits = c(0, 12)
  ) +
  facet_wrap(
    ~t,
    nrow = ROWS
  ) +
  labs(
    title = "All Spots, Limited Intensities",
    x = "",
    y = "",
    fill = "Intensity (x GFP)"
  ) +
  coord_fixed() +
  dark_theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
