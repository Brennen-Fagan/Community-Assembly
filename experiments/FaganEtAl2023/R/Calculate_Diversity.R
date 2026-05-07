Calculate_Diversity <- function(loaded, nspecies) {
  loaded$Abundance[, -1] <-
    loaded$Abundance[, -1] > loaded$Parameters$EliminationThreshold
  
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
  diversity_gamma <- diversity_alpha |> dplyr::group_by(
    Time
  ) |> dplyr::summarise(
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
    
    Richness = unlist(lapply(
      strsplit(
        SpeciesTotal, split = ", ", fixed = TRUE
      ),
      function(x) length(x[x!=""])
    )),
    Richness_Basal = unlist(lapply(
      strsplit(
        SpeciesTotal_Basal, split = ", ", fixed = TRUE
      ),
      function(x) length(x[x!=""])
    )),
    Richness_Consumer = unlist(lapply(
      strsplit(
        SpeciesTotal_Consumer, split = ", ", fixed = TRUE
      ),
      function(x) length(x[x!=""])
    ))
  ) |> dplyr::select(
    -dplyr::starts_with("Species")
  ) |> tidyr::pivot_longer(
    cols = !Time,
    names_to = "Measurement",
    values_to = "Value"
  ) |> dplyr::mutate(
    Environment = "Gamma"
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
      ) |> dplyr::filter(
        Env1 < Env2
      ) |> dplyr::arrange(
        Env1, Env2
      ) |> dplyr::mutate(
        Time = time,
        Jaccard = as.vector(dists)
      )
      
      return(dataf)
    },
    envs = loaded$NumEnvironments
  )
  diversity_beta <- dplyr::bind_rows(
    diversity_beta
  ) |> tidyr::unite(
    "Environment", Env1, Env2, sep = " "
  ) |> dplyr::rename(
    Value = Jaccard
  ) |> dplyr::mutate(
    Measurement = "Jaccard"
  )
  
  diversity_beta_avg <- diversity_beta |> dplyr::group_by(
    Time, Measurement
  ) |> dplyr::summarise(
    Value = mean(Value)
  ) |> dplyr::ungroup(
  ) |> dplyr::mutate(
    Environment = "Mean"
  )
  print("beta")
  
  Diversities <- dplyr::bind_rows(
    diversity_alpha |> dplyr::select(
      -Species_Basal, -Species_Consumer, -Species
    ) |> tidyr::pivot_longer(
      c(-Time, -Environment),
      names_to = "Measurement",
      values_to = "Value"
    ) |> dplyr::mutate(
      Environment = as.character(Environment)
    ),
    diversity_beta,
    diversity_beta_avg,
    diversity_gamma
  )
  
  ### Return Diversities: ###############################################
  return(Diversities)
}