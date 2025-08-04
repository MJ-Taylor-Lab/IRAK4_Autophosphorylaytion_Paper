setwd(RETURN_TO_DIRECTORY)

CellView <-
  ggplot() +
  geom_tile(
    data  = plot_img %>%filter(t == FOCUS_TIME),
    aes(
      x = x,
      y = y,
      fill = z
    )
  ) +
  geom_path(
    data = plot_tracks,
    aes(
      x = POSITION_X,
      y = POSITION_Y,
      group = TRACK_ID
    ),
    color = "green",
    size = 0.25
  ) +
  geom_point(
    data = plot_tracks %>% filter(FRAME == FOCUS_TIME),
    aes(
      x = POSITION_X,
      y = POSITION_Y
    ),
    color = "green",
    shape = 22,
    size = 4
  ) +
  geom_path(
    data = plot_tracks %>% filter(UNIVERSAL_TRACK_ID == FOCUS_TRACK),
    aes(
      x = POSITION_X,
      y = POSITION_Y,
      group = TRACK_ID
    ),
    color = "white",
    size = 0.25
  ) +
  geom_point(
    data = plot_spot %>% filter(t == FOCUS_TIME),
    aes(
      x = POSITION_X,
      y = POSITION_Y
    ),
    color = "white",
    shape = 22,
    size = 4
  ) +
  geom_point(
    data = plot_spot %>% filter(t == FOCUS_TIME),
    aes(
      x = POSITION_X,
      y = POSITION_Y
    ),
    color = "magenta",
    shape = 22,
    size = 4*7
  ) +
  geom_text(
    aes(
      x = max(plot_spot$POSITION_X, plot_spot$POSITION_Y)/6,
      y = max(plot_spot$POSITION_X, plot_spot$POSITION_Y)/9,
      label = paste(FOCUS_TIME, "s")
    ),
    size = max(plot_spot$POSITION_X, plot_spot$POSITION_Y)/5,
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
    title = "Cell View, Limited Intensities",
    x = "µm",
    y = "µm",
    fill = "Intensity (x GFP)"
  ) +
  coord_fixed() +
  dark_theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
