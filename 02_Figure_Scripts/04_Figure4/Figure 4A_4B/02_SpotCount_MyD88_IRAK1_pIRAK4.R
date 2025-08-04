# Reading List ------------------------------------------------------------
Image_List <- Table %>% 
  mutate(
    Final_Path = file.path(file_path, csv_files)
  ) %>% 
  as.data.table()


Spot_Counter <- data.frame(Counter = NA, Count = 0)
Spot_Counter <- Spot_Counter[-1,]

# for(index in 1:nrow(Image_List))
for(index in 1:13){
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
      Counter == 0 ~ "#00FFFF",
      Counter == 1 ~ "#BEBEBE",
      Counter == 2 ~"#CC66FF",
      TRUE ~ "#39B44A"
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


# Write csv ---------------------------------------------------------------
Source_Data_Path <- paste0(basename(Plot_Directory_Save_Path), ".csv")
Plot_Save_Path <- file.path(Plot_Directory_Save_Path, Source_Data_Path)
write.csv(Spot_Counter_compiled, Plot_Save_Path)


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
    width = 0.55,
    size = 0.1
  ) +
  scale_fill_manual(
    values = c("#BEBEBE", "#39B44A", "#CC66FF", "#00FFFF")
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
    limits = c(0,59),
    breaks = seq(0,59, by = 20),
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
  height = 40,
  width = 50,
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
