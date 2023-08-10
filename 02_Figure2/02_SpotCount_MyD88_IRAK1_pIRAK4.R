# Reading List ------------------------------------------------------------
Image_List <- Table %>% 
  mutate(
    Final_Path = file.path(file_path, csv_files)
  ) %>% 
  as.data.table()


Spot_Counter <- data.frame(Counter = NA, Count = 0)
Spot_Counter <- Spot_Counter[-1,]

for(index in 1:nrow(Image_List)){
  csv_path = Image_List$Final_Path[index]
  csv_file = fread(csv_path)
  csv_file <- csv_file %>%  as.data.table()
  csv_file <- csv_file %>% 
    distinct(
      Counter,
      Count
    )
  
  Spot_Counter <- rbind(Spot_Counter, csv_file)
}

rm(
  csv_path,
  csv_file
)

Spot_Counter_compiled <- Spot_Counter %>% 
  group_by(
    Counter
  ) %>% 
  summarise(
    Count_Total = sum(Count)
  ) %>% 
  ungroup() %>% 
  mutate(
    Combination = case_when(
      Counter == 0 ~ "MyD88_only",
      Counter == 1 ~ "MyD88+IRAK1",
      Counter == 2 ~ "MyD88+pIRAK4",
      TRUE ~ "MyD88+pIRAK4+IRAK1"
    ),
    Combination_Colour = case_when(
      Counter == 0 ~ "cyan",
      Counter == 1 ~ "grey",
      Counter == 2 ~"magenta",
      TRUE ~ "red"
    )
  ) %>% 
  mutate(
    Total_Spots = sum(Count_Total)
  ) %>% 
  group_by(
    Counter
  ) %>% 
  mutate(
    Proportion = Count_Total*100/Total_Spots
  ) %>% 
  arrange(
    desc(Proportion)
  ) %>%
  ungroup() %>% 
  group_by(
    Combination
  ) %>% 
  mutate(
    Proportion = signif(Proportion, 4),
    Combination_Name_Full = paste0(Combination, " = ", Proportion)
  )

Spot_Counter_compiled$Combination <- 
  factor(
    Spot_Counter_compiled$Combination,
    levels = c("MyD88+IRAK1", "MyD88+pIRAK4+IRAK1", "MyD88+pIRAK4", "MyD88_only")
  )

# Diet Bar Plot ----------------------------------------------------------------
BarPlot <- ggplot(
  data = Spot_Counter_compiled,
  aes(
    y = Combination,
    x = Proportion,
    fill = Combination
  )
) +
  geom_bar(
    stat = "identity",
    color = "black",
    width = 0.5
  ) +
  scale_fill_manual(
    values = c("grey", "darkgreen", "magenta", "cyan")
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
    y = 2,
    label = Spot_Counter_compiled$Combination_Name_Full[1]
  ) +
  annotate(
    "text",
    x = 15,
    y = 4,
    label = Spot_Counter_compiled$Combination_Name_Full[2]
  ) +
  annotate(
    "text",
    x = 15,
    y = 3,
    label = Spot_Counter_compiled$Combination_Name_Full[3]
  ) +
  annotate(
    "text",
    x = 15,
    y = 1,
    label = Spot_Counter_compiled$Combination_Name_Full[4]
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
    limits = c(0,60),
    breaks = seq(0,60, by = 10),
    position="top",
    expand = c(0, 0)
  ) +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 5),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  )

Plot_Save_Path <- file.path(Plot_Directory_Save_Path, "02_Barplot skinny without labels.pdf")
ggsave(
  Plot_Save_Path,
  plot = last_plot(),
  height = 68,
  width = 60,
  units = "mm"
)


# Cleanup -----------------------------------------------------------------
rm(
  BarPlot,
  Image_List,
  Spot_Counter,
  Spot_Counter_compiled,
  Table2,
  index,
  Plot_Save_Path
)
