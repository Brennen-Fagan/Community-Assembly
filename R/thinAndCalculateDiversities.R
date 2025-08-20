thinAndCalculateDiversities <- function(loaded, nspecies,
                                        preferred_rows_per_event,
                                        divide_time_by) {
  # We can't handle all of the data that we are going to be looking at;
  # a small sample had ~180k rows for ~5.3k events = ~34 rows per event.
  # To reduce it, we will divide time up so that there are about the
  # preferred number of rows per event.
  bythin <- floor((nrow(loaded$Abundance)
                   / nrow(loaded$Events))
                  / preferred_rows_per_event)
  
  loaded$Abundance <- loaded$Abundance[seq(from = 1,
                                           to = nrow(loaded$Abundance),
                                           by = bythin), ]
  
  # Remove illegal values (that the numerical engine uses as inbetweens).
  toEliminate <- loaded$Abundance[, -1] <
    loaded$Parameters$EliminationThreshold & loaded$Abundance[, -1] > 0
  loaded$Abundance[, -1][toEliminate] <- 0
  loaded$Abundance[, 1] <- loaded$Abundance[, 1] / divide_time_by
  
  # Convert to Binary, since we are only interested in richness.
  loaded$Abundance[, -1] <- loaded$Abundance[, -1] > 0
  
  # loaded$Abundance is now prepared.
  # Now, we calculate richnesses.
  print("prep")
  
  ### Alpha Diversity: ##################################################
  diversity_alpha <- lapply(
    1:loaded$NumEnvironments,
    function(i, abund, numSpecies) {
      time <- abund[, 1]
      env <- abund[, 1 + 1:numSpecies + numSpecies * (i - 1)]
      env_basal <- env[, 1:nspecies[1]]
      env_consumer <- env[, nspecies[1] + 1:nspecies[2]]
      richness <- rowSums(env)
      richness_basal <- rowSums(env_basal)
      richness_consumer <- rowSums(env_consumer)
      species <- apply(
        env, MARGIN = 1,
        FUN = function(x) {
          toString(which(x > 0))
        }
      )
      species_basal <- apply(
        env_basal, MARGIN = 1,
        FUN = function(x) {
          toString(which(x > 0))
        }
      )
      species_consumer <- apply(
        env_consumer, MARGIN = 1,
        FUN = function(x) {
          toString(which(x > 0) + nspecies[1]) # + basals doesn't change it, but
          # results in a more intuitive numbering system of species.
        }
      )
      data.frame(Time = time,
                 Richness = richness,
                 Richness_Basal = richness_basal,
                 Richness_Consumer = richness_consumer,
                 Species = species,
                 Species_Basal = species_basal,
                 Species_Consumer = species_consumer,
                 Environment = i,
                 stringsAsFactors = FALSE)
    },
    abund = loaded$Abundance,
    numSpecies = sum(nspecies)
  )
  
  diversity_alpha <- dplyr::bind_rows(diversity_alpha)
  
  print("alpha")
  ### Gamma Diversity: ##################################################
  diversity_gamma <- diversity_alpha %>% dplyr::group_by(
    Time
  ) %>% dplyr::summarise(
    Mean = mean(Richness),
    Mean_Basal = mean(Richness_Basal),
    Mean_Consumer = mean(Richness_Consumer),
    Var = var(Richness),
    Var_Basal = var(Richness_Basal),
    Var_Consumer = var(Richness_Consumer),
    
    SpeciesTotal = toString(sort(unique(unlist(strsplit(paste(
      Species, collapse = ", "), split = ", ", fixed = TRUE))))),
    SpeciesTotal_Basal = toString(sort(unique(unlist(strsplit(paste(
      Species_Basal, collapse = ", "), split = ", ", fixed = TRUE))))),
    SpeciesTotal_Consumer = toString(sort(unique(unlist(strsplit(paste(
      Species_Consumer, collapse = ", "), split = ", ", fixed = TRUE))))),
    
    Gamma = unlist(lapply(
      strsplit(
        SpeciesTotal, split = ", ", fixed = TRUE
      ),
      function(x) length(x[x!=""])
    )),
    Gamma_Basal = unlist(lapply(
      strsplit(
        SpeciesTotal_Basal, split = ", ", fixed = TRUE
      ),
      function(x) length(x[x!=""])
    )),
    Gamma_Consumer = unlist(lapply(
      strsplit(
        SpeciesTotal_Consumer, split = ", ", fixed = TRUE
      ),
      function(x) length(x[x!=""])
    ))
  ) %>% tidyr::pivot_longer(
    cols = c(Mean, Var, Gamma),
    names_to = "Aggregation",
    values_to = "Richness"
  ) %>% dplyr::mutate(
    Basals = case_when(
      Aggregation == "Mean" ~ Mean_Basal,
      Aggregation == "Var" ~ Var_Basal,
      Aggregation == "Gamma" ~ as.numeric(Gamma_Basal),
      TRUE ~ -1
    ),
    Consumers = case_when(
      Aggregation == "Mean" ~ Mean_Consumer,
      Aggregation == "Var" ~ Var_Consumer,
      Aggregation == "Gamma" ~ as.numeric(Gamma_Consumer),
      TRUE ~ -1
    )
  ) %>% dplyr::select(
    -Mean_Basal, -Var_Basal, -Gamma_Basal,
    -Mean_Consumer, -Var_Consumer, -Gamma_Consumer,
    -SpeciesTotal, -SpeciesTotal_Basal, -SpeciesTotal_Consumer
  )
  
  print("gamma")
  ### Beta Diversity (Jaccard, Space): ##################################
  diversity_beta <- apply(
    loaded$Abundance,
    MARGIN = 1, # Rows
    function(row, envs) {
      time <- row[1]
      # Vegan complains about rows with all 0's.
      # The warning is generic, so we cannot silence it specifically.
      dists <- suppressWarnings(vegan::vegdist(
        method = "jaccard",
        x = matrix(row[-1] > 0, nrow = envs, byrow = TRUE)
      ))
      
      dataf <- expand.grid(
        Env1 = 1:envs,
        Env2 = 1:envs
      ) %>% dplyr::filter(
        Env1 < Env2
      ) %>% dplyr::arrange(
        Env1, Env2
      ) %>% dplyr::mutate(
        Time = time,
        Jaccard = as.vector(dists)
      )
      
      return(dataf)
    },
    envs = loaded$NumEnvironments
  )
  print("beta")
  
  ### Return Diversities: ###############################################
  return(list(
    alpha = diversity_alpha %>% dplyr::select(
      -Species_Basal, -Species_Consumer
    ),
    beta = diversity_beta,
    gamma = diversity_gamma
  ))
}