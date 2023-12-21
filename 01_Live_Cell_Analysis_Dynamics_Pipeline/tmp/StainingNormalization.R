ChannelFrequency <-
  IntensitiesRaw %>% 
  group_by(
    Stain,
    Condition,
    Cell_Line,
    Replicate,
    Frame,
    GFP
  ) %>% 
  summarize(
    N = sum(N)
  ) %>% 
  group_by(
    Stain,
    Condition,
    Cell_Line,
    Replicate,
    Frame
  ) %>% 
  mutate(
    Mode = getmode(GFP),
    GFP = GFP - Mode
  ) %>% 
  filter(
    GFP >= 0
  ) %>% 
  mutate(
    N = N/max(N),
    P25 = diagis::weighted_quantile(GFP, N, .25),
    GFP = GFP - P25,
    # Median = diagis::weighted_quantile(GFP, N, .5),
    # GFP = GFP/Median*1000,
    P85 = diagis::weighted_quantile(GFP, N, .85),
    GFP = GFP/P85*1000,
    P95 = diagis::weighted_quantile(GFP, N, .95)
  ) %>% 
  filter(
    GFP >= 0,
    GFP <= P95
  )

ChannelFrequency %>% 
  group_by(
    Frame
  ) %>% 
  summarize(
    Min = min(GFP),
    Max = max(GFP),
    Median = median(GFP)
  )

ToNormalize <-
  IntensitiesRaw %>% 
  group_by(
    Frame
  ) %>% 
  mutate(
    ID = 1:n()
  ) %>% 
  select(
   GFP,
   ID,
   Frame
  ) %>% 
  pivot_wider(
    names_from = Frame,
        id_cols = c(
          "ID"
        ),
    values_from = GFP
  ) %>% 
  select(-c(
    ID
  ))


quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

Test <- quantile_normalisation(ToNormalize)
Test <- as_tibble(Test)


frame_quantile_normalize <- function(TableX){
  TableX$GFP <- quantile_normalize(TableX$GFP)
  return(TableX)
}
Results <- mclapply(ToNormalize, frame_quantile_normalize)
Results <- rbindlist(Results)

Test2 <-
  Test %>% 
  pivot_longer(
    cols = everything()
  )

Test <- NULL
Test$Frame <- Test2$name
Test$GFP <- Test2$value
Test <- as_tibble(Test)

ggplot(
  Test,
  aes(
    x = GFP,
    y = ..ndensity..,
    group = Frame,
    color = Frame
  )
) +
  geom_density()+
  scale_x_log10() +
  dark_theme_classic()


ggplot(
  ChannelFrequency,
  aes(
    x = GFP,
    y = N,
    group = as.character(Frame),
    color = as.character(Frame)
  )
) +
  geom_path() +
  # facet_wrap(~Frame)+
  # scale_x_log10()+
  # scale_y_log10()+
  # facet_grid(
  #   Cell_Line+Condition~Replicate+Frame
  # ) +
  dark_theme_classic()
