
# LinearPhasePortrait <- fread(paste0("/Users/u_deliz/Desktop/PhasePortraitAnalysis/PhasePortrait - LeadLag 5 - StepSize ", STEP_SIZE, ".csv.gz"))
# LinearPhasePortrait <- LinearPhasePortrait %>% filter(QUERY_PROTEIN != "TRAF6") %>% as.data.table()

FRAMES_SINCE_LANDING_CAT_BINS <- c(50)

# Variables list
Facets <- unique(LinearPhasePortrait$PLOT_FACETS)
x_variables <- c("ROUNDED_REFERENCE_TOTAL_INTENSITY", "ROUNDED_QUERY_TOTAL_INTENSITY")
y_variables <-c(
  "DELTA_REFERENCE_TOTAL_INTENSITY", "DELTA_QUERY_TOTAL_INTENSITY",
  # "ADJUSTED_DELTA_REFERENCE_TOTAL_INTENSITY", "ADJUSTED_DELTA_QUERY_TOTAL_INTENSITY",
  "FRAMES_ADJUSTED", "N",
  "DELTA_RATIO"
  # "ADJUSTED_DELTA_RATIO"
)
z_variables <- rev(x_variables)

PlotList <- NULL
# x-axis
rep_x_variables = rep(x_variables, each = NROW(y_variables))
rep_x_variables = rep(rep_x_variables, NROW(Facets))
PlotList$x_axis = rep_x_variables
# y-axis
rep_y_variables = rep(y_variables, NROW(x_variables))
PlotList$y_axis = rep(rep_y_variables, NROW(Facets))
# z-axis
rep_z_variables = rep(z_variables, each = NROW(y_variables))
rep_z_variables = rep(rep_z_variables, NROW(Facets))
PlotList$z_axis = rep_z_variables
# category
Facets <- sort(Facets)
PlotList$Facets = rep(Facets, each = NROW(x_variables)*NROW(y_variables))
# Combine
PlotList <- as_tibble(PlotList)

PlotFx <- function(FacetX){
  tryCatch({
    
    if(PlotList$y_axis[FacetX] == "N" || PlotList$y_axis[FacetX] == "FRAMES_ADJUSTED"||
       PlotList$y_axis[FacetX] == "THETA" ||PlotList$y_axis[FacetX] == "ADJUSTED_THETA"||
       PlotList$y_axis[FacetX] == "DELTA_RATIO" ||PlotList$y_axis[FacetX] == "ADJUSTED_DELTA_RATIO"
       ){
      mult_format <- function(l) {l}
    } else{
      mult_format <- function(l) {
        l <- format(l*1000, digits = 3, scientific = F)
        l <- paste0("'", l, "'%*%10^-3")
        parse(text=l)
      }
    }
    
    TempTable <-
      LinearPhasePortrait %>%
      as_tibble() %>% 
      filter(
        PLOT_FACETS == PlotList$Facets[FacetX],
        N >= 50
      ) %>% 
      as.data.table()
    
    if(NROW(grep("DELTA_", PlotList$y_axis[FacetX]))==1){
      # y_axis_min = gsub("_TOTAL_INTENSITY", "_Q1", PlotList$y_axis[FacetX])
      # y_axis_max = gsub("_TOTAL_INTENSITY", "_Q3", PlotList$y_axis[FacetX])
      y_axis_mad = paste0("MAD_", PlotList$y_axis[FacetX])
    } else{
      # y_axis_min = PlotList$y_axis[FacetX]
      # y_axis_max = PlotList$y_axis[FacetX]
      y_axis_mad = PlotList$y_axis[FacetX]
    }

    TempTable <-
      TempTable %>%
      filter(
        FRAMES_SINCE_LANDING_CAT == 100,
        # FRAMES_SINCE_LANDING_CAT %in% FRAMES_SINCE_LANDING_CAT_BINS,
        N_TEST == TRUE
      ) %>%
      group_by(
        PLOT_FACETS, IMAGENUMBER
      ) %>%
      mutate(
        # Define plot axis replacing variables with generic names
        x_axis = (!!as.name(PlotList$x_axis[FacetX])),
        y_axis = (!!as.name(PlotList$y_axis[FacetX])),
        z_axis = (!!as.name(PlotList$z_axis[FacetX])),
        # y_axis_min = (!!as.name(y_axis_min)),
        # y_axis_max = (!!as.name(y_axis_max))
        y_axis_sem =(!!as.name(y_axis_mad))/sqrt(N)
      ) %>%
      filter(
        x_axis <= 20,
        z_axis <= 8
      ) %>%
      as.data.table()
    
    # Labels
    x_title <- PlotList$x_axis[FacetX]
    x_title <- gsub("_TOTAL_INTENSITY", " Norm. Int.", x_title)
    x_title <- gsub("ROUNDED_", "", x_title)
    x_title <- gsub("REFERENCE", TempTable$REFERENCE_PROTEIN[1], x_title)
    x_title <- gsub("QUERY", TempTable$QUERY_PROTEIN[1], x_title)
    
    y_title <- PlotList$y_axis[FacetX]
    y_title <- gsub("_TOTAL_INTENSITY", " Norm. Int.", y_title)
    y_title <- gsub("ROUNDED_", "", y_title)
    y_title <- gsub("REFERENCE", TempTable$REFERENCE_PROTEIN[1], y_title)
    y_title <- gsub("QUERY", TempTable$QUERY_PROTEIN[1], y_title)
    y_title <- gsub("DELTA_", "Change ", y_title)
    y_title <- gsub("ADJUSTED_", "Scaled ", y_title)
    y_title <- gsub("FRAMES_ADJUSTED", "Age", y_title)
    y_title <- gsub("RATIO", "Ratio", y_title)
    
    z_title <- PlotList$z_axis[FacetX]
    z_title <- gsub("_TOTAL_INTENSITY", " Norm. Int.", z_title)
    z_title <- gsub("ROUNDED_", "", z_title)
    z_title <- gsub("REFERENCE", TempTable$REFERENCE_PROTEIN[1], z_title)
    z_title <- gsub("QUERY", TempTable$QUERY_PROTEIN[1], z_title)
    z_title <- gsub("DELTA_", "Change ", z_title)
    z_title <- gsub("ADJUSTED_", "Scaled ", z_title)
    z_title <- gsub("FRAMES_ADJUSTED", "Age ", z_title)
    
    PlotTitle <-
      bquote(
        "Cell Line " *
          .(TempTable$COHORT[1]) *
          " ~ " *
          .( TempTable$LIGAND_DENSITY_CAT[1]) *
          " mol. " *
          µm^2 *
          " ~ " *
          .(TempTable$FPS[1]) *
          " ~ Rep. "*
          .(TempTable$IMAGENUMBER[1])
      )
    
    InterpolationTable <-
      TempTable %>% 
      select(
        x_axis,
        y_axis,
        y_axis_sem,
        z_axis
      ) %>% 
      as_tibble() %>% 
      group_by(
        z_axis
      ) %>% 
      complete(
        x_axis = full_seq(c(min(x_axis), max(x_axis)+STEP_SIZE), period = STEP_SIZE/21),
      ) %>% 
      fill(
        one_of(
          "y_axis",
          "y_axis_sem"
        ),
        .direction = "downup"
      ) %>% 
      group_by(
        z_axis
      ) %>% 
      mutate(
        y_axis_old = y_axis,
        y_axis_sem_old = y_axis_sem,
        y_axis = roll_mean(y_axis, 11, fill = NA, align = "left"),
        y_axis_sem = roll_mean(y_axis_sem, 11, fill = NA, align = "left"),
        y_axis = roll_mean(y_axis, 11, fill = NA, align = "left"),
        y_axis_sem = roll_mean(y_axis_sem, 11, fill = NA, align = "left")
      ) %>% 
      tidyr::drop_na()
    
    
    ggplot() +
      geom_hline(
        yintercept = 0,
        size = 2,
        color = "black"
      ) +
      geom_point(
        data = TempTable,
        aes(
          x = x_axis,
          y = y_axis,
          color = z_axis,
          fill = z_axis,
          group = z_axis,
          size = N
        )
      ) +
      # geom_ribbon(
      #   aes(
      #     x = x_axis,
      #     y = y_axis,
      #     color = z_axis,
      #     fill = z_axis,
      #     group = z_axis,
      #     ymin = y_axis - y_axis_sem,
      #     ymax =  y_axis + y_axis_sem
      #   ),
      #   alpha = .25
    # ) +
    geom_ribbon(
      data = InterpolationTable,
      aes(
        x = x_axis,
        y = y_axis,
        color = z_axis,
        fill = z_axis,
        group = z_axis,
        ymin = y_axis - y_axis_sem,
        ymax =  y_axis + y_axis_sem
      ),
      alpha = .25
    ) +
      scale_size_continuous(
        trans = log2_trans(),
        breaks = trans_breaks("log2", function(x) 2^x),
        # labels = trans_format("log2", math_format(2^.x))
        guide = "none"
      ) +
      scale_color_viridis(
        option = "plasma"
      ) +
      scale_fill_viridis(
        option = "plasma"
      ) +
      scale_y_continuous(
        labels = mult_format
      ) +
      labs(
        title = PlotTitle,
        x = x_title,
        y = y_title,
        color = gsub(" Norm. Int.", "\nNorm. Int.", z_title),
        fill = gsub(" Norm. Int.", "\nNorm. Int.", z_title)
      ) +
      # facet_grid(
      #   ~z_axis
      # ) +
      theme_classic(
        base_size = 20
      ) +
      theme(
        legend.key.height= unit(1.25, 'cm'),
        legend.position = "right",
        rect = element_rect(fill = "transparent")
      )
    
    # Save folder
    COHORT_FOLDER = file.path(OUTPUT_DIRECTORY,
                              paste0(
                                "Cell Line ", TempTable$COHORT[1], " - ",
                                TempTable$LIGAND_DENSITY_CAT[1],  " mol. µm^2 - ",
                                TempTable$FPS[1], " Hz"
                              ))
    
    if(!file.exists(COHORT_FOLDER)){
      dir.create(COHORT_FOLDER)
    }
    
    # File name
    SaveName <-
      paste0(
        "1D - ",
        "y ", y_title, " - ",
        "x ", x_title, " - ",
        "z ", z_title, " - ",
        # "LeadLag ", LEAD_LAG, " - ",
        # "StepSize ", round(STEP_SIZE, 2), " - ",
        "Replicate ", TempTable$IMAGENUMBER[1],
        ".pdf"
      )
    
    SaveName <- gsub(" - \\.", "\\.", SaveName)
    
    ggsave(
      # Save vector image
      file.path(COHORT_FOLDER, SaveName),
      # height = 4.76*1.25,
      # width = 11.5*1.25
      height = 4.76,
      width = 11.5
    )
    
  }, error = function(e) {print(paste("Error with PlotFx. FacetX =",FacetX))})
}
PlotPaths <- mclapply(1:NROW(PlotList), PlotFx, mc.cores = detectCores(logical = F))
PlotPaths <- unlist(PlotPaths)

