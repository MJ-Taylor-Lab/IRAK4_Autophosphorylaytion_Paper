setwd(OUTPUT_DIRECTORY)

Facets <- unique(LinearPhasePortrait$PLOT_FACETS)

FRAMES_SINCE_LANDING_CAT_BINS <- c(50, 100)
FRAMES_SINCE_LANDING_CAT_BINS <- 100

Variables <-
  c(
    # "MAD_ADJUSTED_",
    "ADJUSTED_",
    ""
  )

PlotList <- NULL
PlotList$Facets = rep(Facets, each = NROW(Variables))
PlotList$Variables = rep(Variables, NROW(Facets))
PlotList <- as.data.table(PlotList)
PlotList <- PlotList %>% arrange(Facets) %>% as_tibble()
remove(Facets, Variables)

PlotFx <- function(FacetX){
  tryCatch({

    TempVariables <-
      paste0(
        PlotList$Variables[FacetX],
        c(
          "DELTA_REFERENCE_TOTAL_INTENSITY",
          "DELTA_QUERY_TOTAL_INTENSITY"
        )
      )

    TempTable <-
      LinearPhasePortrait %>%
      select(-c(
        FPS
      )) %>% 
      filter(
        FRAMES_SINCE_LANDING_CAT %in% FRAMES_SINCE_LANDING_CAT_BINS,
        PLOT_FACETS == PlotList$Facets[FacetX]
      ) %>%
      group_by(
        PLOT_FACETS, IMAGENUMBER
      ) %>%
      mutate(
        # Define plot axis replacing variables with generic names
        ref_int = ROUNDED_REFERENCE_TOTAL_INTENSITY,
        qry_int = ROUNDED_QUERY_TOTAL_INTENSITY,

        d_ref_int = (!!as.name(TempVariables[1])),
        d_qry_int = (!!as.name(TempVariables[2])),
        
        FPS = paste(as.character(LinearPhasePortrait$FPS[[1]]), "Hz")
      ) %>% 
      as.data.table()

    GreyTempTable <-
      TempTable %>%
      filter(
        N_TEST == FALSE
      ) %>%
      as.data.table()
    
    TempTable <-
      TempTable %>%
      filter(
        N_TEST == TRUE
      ) %>%
      as.data.table()

    PlotTitle <-
      paste0(
        "Cell Line ", TempTable$COHORT[[1]], " - ",
        TempTable$LIGAND_DENSITY_CAT[[1]],  " mol. µm^2 - ",
        TempTable$FPS[[1]],
        " - Replicate ",TempTable$IMAGENUMBER[[1]]
      )

    PlotTitle <-
      bquote(
        "Cell Line " *
          .(TempTable$COHORT[[1]]) *
          " ~ " *
          .( TempTable$LIGAND_DENSITY_CAT[[1]]) *
          " mol. " *
          µm^2 *
          " ~ " *
          .(TempTable$FPS[[1]]) *
          " ~ Rep. "*
          .(TempTable$IMAGENUMBER[[1]])
      )

    StreamTempTable <-
      TempTable %>%
      # bind_rows(TempTable, GreyTempTable) %>%
      as_tibble() %>%
      group_by(
        FRAMES_SINCE_LANDING_CAT
      ) %>%
      complete(
        ref_int = full_seq(c(min(ref_int), max(ref_int)), period = STEP_SIZE),
        qry_int = full_seq(c(min(qry_int), max(qry_int)), period = STEP_SIZE)
      ) %>%
      as.data.table()

    # Add derivatives
    StreamTempTable$d_ref_int[is.na(StreamTempTable$d_ref_int)] = 0
    StreamTempTable$d_qry_int[is.na(StreamTempTable$d_qry_int)] = 0

    if(NROW(GreyTempTable)==0){
      ggplot(
      ) +
        geom_streamline(
          data = StreamTempTable,
          aes(
            x = ref_int,
            y = qry_int,
            dx = d_ref_int,
            dy = d_qry_int,
            color = sqrt(d_ref_int^2 + d_qry_int^2),
            size = ..step..
            # alpha = ..step..
          ),
          arrow = NULL,
          # n = 5,
          # arrow.length = 0.3,
          # jitter = 4,
          L = 50, res = 10, lineend = "round"
        ) +
        geom_arrow(
          data = TempTable,
          aes(
            x = ref_int,
            y = qry_int,
            mag = 0.5,
            angle = atan2(d_qry_int, d_ref_int)*180/pi,
            # color = log(N+1)
            color = sqrt(d_ref_int^2 + d_qry_int^2)
            # color = FRAMES_ADJUSTED
          ),
          # color = "white",
          # size = .01, # Small arrow head
          # arrow.length = .5, # Small arrow head

          size = 1, # Big arrow head
          arrow.length = 1,# Big arrow head
          lineend = "square"
        ) +
        scale_size(range = c(.01, 1), guide = "none") +
        scale_alpha(guide = "none") +
        scale_mag(
          # max = 3, # Small arrow head
          max = 1, # Big arrow head
          guide = 'none'
        ) +
        scale_color_viridis(
          option = "plasma",
          # direction = -1,
          # trans = "log10",
          # breaks = trans_breaks("log10", function(x) 10^x),
          # labels = trans_format("log10", math_format(10^.x))
        ) +
        labs(
          title = PlotTitle,
          x = paste(TempTable$REFERENCE_PROTEIN[1], "Norm. Int."),
          y = paste(TempTable$QUERY_PROTEIN[1], "Norm. Int."),
          color = "N"
          # color = "Scaled\nMagnitude",
          # color = "Time\n(Frames)"
        ) +
        # Make dark to make the colors pop in a monitor
        # Use theme classic for printing
        theme_classic(
          base_size = 20
        ) +
        theme(
          legend.key.height= unit(1.25, 'cm'),
          legend.position = "right",
          rect = element_rect(fill = "transparent")
        ) +
        coord_fixed(
          xlim = c(0, max(TempTable$ref_int, GreyTempTable$ref_int)),
          ylim = c(0, max(TempTable$qry_int, GreyTempTable$qry_int))
        )

    } else{
      
      ggplot(
      ) +
        geom_streamline(
          data = StreamTempTable,
          aes(
            x = ref_int,
            y = qry_int,
            dx = d_ref_int,
            dy = d_qry_int,
            color = sqrt(..dx..^2 + ..dy..^2),
            size = ..step..
            # alpha = ..step..
          ),
          arrow = NULL,
          # n = 15,
          # arrow.length = 0.3,
          # jitter = 4,
          L = 10, res = 10, lineend = "round"
        ) +
        geom_arrow(
          data = GreyTempTable,
          aes(
            x = ref_int,
            y = qry_int,
            mag = 0.5,
            angle = atan2(d_qry_int, d_ref_int)*180/pi
          ),
          color = "black", # "#3A3B3C",
          # size = .01, # Small arrow head
          # arrow.length = .5, # Small arrow head

          size = 1, # Big arrow head
          arrow.length = 1,# Big arrow head
          lineend = "square"
        ) +
        geom_arrow(
          data = TempTable,
          aes(
            x = ref_int,
            y = qry_int,
            mag = 0.5,
            angle = atan2(d_qry_int, d_ref_int)*180/pi,
            color = sqrt(d_qry_int^2 + d_ref_int^2),
            # color = FRAMES_ADJUSTED
          ),
          # size = .01, # Small arrow head
          # arrow.length = .5, # Small arrow head

          size = 1, # Big arrow head
          arrow.length = 1,# Big arrow head
          lineend = "square"
        ) +
        scale_size(range = c(.01, 1), guide = "none") +
        scale_alpha(guide = "none") +
        scale_mag(
          # max = 3, # Small arrow head
          max = 1, # Big arrow head
          guide = 'none'
        ) +
        scale_color_viridis(
          option = "plasma",
          # direction = -1,
          # trans = "log10",
          # breaks = trans_breaks("log10", function(x) 10^x),
          # labels = trans_format("log10", math_format(10^.x))
        ) +
        labs(
          title = PlotTitle,
          x = paste(TempTable$REFERENCE_PROTEIN[1], "Norm. Int."),
          y = paste(TempTable$QUERY_PROTEIN[1], "Norm. Int."),
          color = "Magnitude"
          # color = "Scaled\nMagnitude",
          # color = "Time\n(Frames)"
        ) +
        facet_grid(
          ~FRAMES_SINCE_LANDING_CAT
        ) +
        # Make dark to make the colors pop in a monitor
        # Use theme classic for printing
        theme_classic(
          base_size = 20
        ) +
        theme(
          legend.key.height= unit(1.25, 'cm'),
          legend.position = "right",
          rect = element_rect(fill = "transparent")
        ) +
        coord_fixed(
          # xlim = c(0, max(TempTable$ref_int, GreyTempTable$ref_int)),
          # ylim = c(0, max(TempTable$qry_int, GreyTempTable$qry_int))
        )
      
    }

    # Save folder
    COHORT_FOLDER = file.path(OUTPUT_DIRECTORY,
                              paste0(
                                "Cell Line ", TempTable$COHORT[[1]], " - ",
                                TempTable$LIGAND_DENSITY_CAT[[1]],  " mol. µm^2 - ",
                                TempTable$FPS[[1]]
                              ))

    if(!file.exists(COHORT_FOLDER)){
      dir.create(COHORT_FOLDER)
    }

    # File name
    SaveName <-
      paste0(
        "Normalized ",
        PlotList$Variables[FacetX],
        "PhasePortrait - ",
        "LeadLag ", LEAD_LAG, " - ",
        "StepSize ", round(STEP_SIZE, 2), " - ",
        "Replicate ",TempTable$IMAGENUMBER[[1]],
        ".pdf"
      )

    SaveName <- gsub(" - \\.", "\\.", SaveName)

    ggsave(
      # Save vector image
      file.path(COHORT_FOLDER, SaveName),
      # height = 4.76*1.25,
      # width = 11.5*1.25
      height = 4.76*1.5,
      width = 11.5*1.5
    )
    
    print(file.path(COHORT_FOLDER, SaveName))
  }, error = function(e) {print(paste("Error with PlotFx", PlotList$Variables[FacetX], PlotList$Facets[FacetX]))})}
mclapply(1:NROW(PlotList), PlotFx, mc.cores = detectCores(logical = F))
