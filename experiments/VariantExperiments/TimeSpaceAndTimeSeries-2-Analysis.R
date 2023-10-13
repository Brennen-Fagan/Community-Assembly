calculationsBII <- FALSE
# Skip glmer.nb for BII from Scholes and Biggs?
calculationsLong <- FALSE
calculationsPlotLong <- FALSE

# Data Analysis: ##############################################################
# We'll want equivalents of the plots in the original TimeSpaceAndTimeSeries.R.
# Ines's conception of the summary statistic is the difference between Control
# and Experiment.
# This can be at the patch level (\Delta \alpha) or at the set of patches level
# (equivalent to \Delta \gamma).

### Alpha:#####################################################################
# Note the order of operations is very relevant here.
# We imagine that the researchers consider all samples across times to be
# samples from the same population (even if it isn't in reality), so long as
# the control and experiment are respected.
# We then believe they would look at the difference between the control and the
# experiment.
# Edit: This reflects In\'es's idea that it is the difference between number of
# species between control and experiment.
bootstrapSamplesDeltaAlpha <- bootstrapSamples %>% dplyr::group_by(
  # Average across patches first.
  Type, Control, Bootstrap
) %>% dplyr::summarise(
  AverageSamplingAlpha = mean(SamplingAlpha),
  AverageSamplingAlphaNative = mean(SamplingAlphaNative),
  AverageSamplingAlphaInvasive = mean(SamplingAlphaInvasive),
  .groups = "drop"
) %>% dplyr::group_by(
  # Then perform difference by converting control to negatives and adding.
  Type, Bootstrap
) %>% dplyr::mutate(
  dplyr::across(
    .cols = AverageSamplingAlpha : AverageSamplingAlphaInvasive,
    .fns = ~ ifelse(Control == "Control", -.x, .x)
  )
) %>% dplyr::summarise(
  DeltaAverageSamplingAlpha = sum(AverageSamplingAlpha),
  DeltaAverageSamplingAlphaNative = sum(AverageSamplingAlphaNative),
  DeltaAverageSamplingAlphaInvasive = sum(AverageSamplingAlphaInvasive),
  .groups = "drop"
) %>% dplyr::rename(
  `Overall` = DeltaAverageSamplingAlpha,
  `Detected in Control` = DeltaAverageSamplingAlphaNative,
  `Not Detected in Control` = DeltaAverageSamplingAlphaInvasive
) %>% tidyr::pivot_longer(
  cols = c(`Overall`, `Detected in Control`, `Not Detected in Control`),
  names_to = "Species Subset",
  values_to = "Difference of Average Number of Species in Patch"
)

### Gamma: ####################################################################
# (Note: no guarantee of agreement due to differing sampling.)
bootstrapSamplesDeltaGamma <- bootstrapSamples %>% dplyr::group_by(
  # "Average" (Collect) across patches first.
  Type, Control, Bootstrap
) %>% dplyr::group_modify(
  .f = function(x, y) {
    IDsNative <- unique(unlist(strsplit(
      x$SamplingIDsNative, ", ", fixed = TRUE
    )))
    IDsInvasive <- unique(unlist(strsplit(
      x$SamplingIDsInvasive, ", ", fixed = TRUE
    )))
    
    data.frame(
      SamplingGammaNative = length(IDsNative),
      IDsNative = toString(IDsNative),
      SamplingGammaInvasive = length(IDsInvasive),
      IDsInvasive = toString(IDsInvasive)
    ) %>% dplyr::mutate(
      SamplingGamma = SamplingGammaNative + SamplingGammaInvasive
    )
  }
) %>% dplyr::group_by(
  # Then perform difference by converting control to negatives and adding.
  Type, Bootstrap
) %>% dplyr::mutate(
  dplyr::across(
    .cols = dplyr::starts_with("SamplingGamma"),
    .fns = ~ ifelse(Control == "Control", -.x, .x)
  )
) %>% dplyr::summarise(
  DeltaSamplingGamma = sum(SamplingGamma),
  DeltaSamplingGammaNative = sum(SamplingGammaNative),
  DeltaSamplingGammaInvasive = sum(SamplingGammaInvasive),
  .groups = "drop"
)  %>% dplyr::rename(
  `Overall` = DeltaSamplingGamma,
  `Detected in Control` = DeltaSamplingGammaNative,
  `Not Detected in Control` = DeltaSamplingGammaInvasive
) %>% tidyr::pivot_longer(
  cols = c(`Overall`, `Detected in Control`, `Not Detected in Control`),
  names_to = "Species Subset",
  values_to = "Difference of Number of Species across Patches"
)

### Alpha Slope: ##############################################################
# Edit: In\'es also suggested something different from what I was expecting.
# In contrast to looking at the raw diversities, she said quite a few
# researchers are instead looking at the change in number of species through
# time. This probably entails fitting a linear model and reporting the slope.
# We can also decompose again into "native" slope and "invasive" slope.
# We'll start with the most basic linear model, although a better one would
# probably look like a Bayesian estimate of the number of categories in a
# multinomial that has Poisson sampling

bootstrapSamplesLMAlpha <- bootstrapSamples %>% dplyr::group_by(
  Type, Bootstrap
) %>% dplyr::group_modify(
  .f = function(.x, .y) {
    fit1 <- lme4::lmer(SamplingAlpha ~ TimeSinceStart : Control + (1 | Patch),
                       data = .x)
    fit2 <- lme4::lmer(SamplingAlphaNative ~ TimeSinceStart : Control + (1 | Patch),
                       data = .x)
    fit3 <- lme4::lmer(SamplingAlphaInvasive ~ TimeSinceStart : Control + (1 | Patch),
                       data = .x)
    # Pool Intercepts Across Patches
    cbind(data.frame(
      Intercept = rep(c(fit1@beta[1],
                        fit2@beta[1],
                        fit3@beta[1]), each = 2),
      Slope = c(fit1@beta[2], fit1@beta[3],
                fit2@beta[2], fit2@beta[3],
                fit3@beta[2], fit3@beta[3]),
      Control = rep(c("Control", "Experiment"), 3),
      Subset = rep(c("Overall",
                     "Detected in Control",
                     "Not Detected in Control"), each = 2)
    ))
  }
)

### Intactness and beta diversity: ############################################
if (calculationsBII) {
  ##### Source Data: ############################################################
  # https://adrianadepalma.github.io/BII_tutorial/bii_example.html#About_BII
  # Looking at the data they suggest, the equivalences are roughly:
  #  Species IDs =  c(Study_common_taxon, Rank, Taxon_number : Higher_taxon)
  #  Abundance = c(Diversity_metric..., Measurement:Effort_corrected_measurement)
  #  Time = c(Sample_start_earliest:Sample_date_resolution, Years_since_frag...)
  #  Patch = c(Study_number:SSBS, Transect_details:Wilderness_area)
  #  Control = Habitat_as_described:Use_intensity
  # Note: none of the time variables they use imply repeated sampling treated
  # as separate instances. Instead, they imply that they are combined.
  
  bootstrapSamplesBIISource <- bootstrapSamples %>% dplyr::group_by(
    # First format similarly to BII input data.
    Time, Patch, Type, Control, Bootstrap # together == Rowwise
  ) %>% dplyr::group_modify(
    .f = function(x, y) {
      IDs <- as.numeric(strsplit(x$SamplingIDs, ", ", fixed = TRUE)[[1]])
      if (length(IDs) == 0) return(data.frame())
      IDsInvasive <- as.numeric(strsplit(
        x$SamplingIDsInvasive, ", ", fixed = TRUE)[[1]])
      
      data.frame(table(IDs)) %>% dplyr::rename(
        ID = IDs, Abundance = Freq
      ) %>% dplyr::mutate(
        ID = as.numeric(as.character(ID)),
        Invasive = ID %in% IDsInvasive
      )
    }
  ) %>% dplyr::group_by(
    ID, Patch, Type, Control, Bootstrap
  ) %>% dplyr::summarise(
    # Across Time, estimate abundance, see note above.
    Abundance = sum(Abundance),
    Invasive = toString(unique(Invasive)),
    .groups = "drop"
  )
  
  ##### Abundance Data: #########################################################
  # Quite a weird thing here: they effectively mean to be summarising, but don't
  # in the guidance online (ID and Abundance are not used anymore and are
  # discarded for all species except one within each site). For safety and sanity
  # we discard them explicitly here, so if an error crops up, we will know.
  bootstrapSamplesBIIAbundance <- bootstrapSamplesBIISource %>% dplyr::group_by(
    # Rescale abundance by the maximum [abundance of all species] over patches.
    Patch, Type, Control, Bootstrap # "unique sites"
  ) %>% dplyr::mutate(
    TotalAbundance = sum(Abundance)
  ) %>% dplyr::ungroup(
  ) %>% dplyr::distinct(
    # This is effectively where the summary collapses the groups.
    Patch, Type, Control, Bootstrap, # "unique sites"
    .keep_all = TRUE
  ) %>% dplyr::group_by(
    Type, Bootstrap # "Study ID", not sure if Control should be in here.
  ) %>% dplyr::mutate(
    MaxAbundance = max(TotalAbundance)
  ) %>% dplyr::ungroup(
  ) %>% dplyr::mutate(
    RescaledAbundance = TotalAbundance/MaxAbundance
  ) %>% dplyr::select(-ID, -Abundance)
  
  # Interesting plots: it generally doesn't look like abundance is going up,
  # although there are some occasional spikes.
  # ggplot2::ggplot(
  #   bootstrapSamples,
  #   ggplot2::aes(x = Time, y = SamplingAbundance,
  #                group = interaction(Patch, Type, Control, Bootstrap),
  #                color = Bootstrap)
  #   ) + ggplot2::geom_line() + ggplot2::facet_wrap(Patch ~ Type)
  # Despite that, their measure of abundances
  # bootstrapSamplesBIIAbundance %>% ggplot2::ggplot(
  #   ggplot2::aes(x = interaction(Control, Type), y = RescaledAbundance)
  #   ) + ggplot2::geom_violin() + ggplot2::geom_boxplot(notch = TRUE)
  # And how I would do so if asked
  # bootstrapSamplesBIISource %>% dplyr::group_by(
  #   Type, Control, Bootstrap
  #   ) %>% dplyr::summarise(sum(Abundance)) %>% ggplot2::ggplot(
  #     ggplot2::aes(x = interaction(Control, Type), y = `sum(Abundance)`)
  #     ) + ggplot2::geom_violin() + ggplot2::geom_boxplot(notch = TRUE)
  # both show substantial differences, and theirs are exaggerated.
  # Nonetheless, I don't think *mine* are significant.
  
  ##### Site Data: ##############################################################
  bootstrapSamplesBIISites <- bootstrapSamplesBIISource %>% dplyr::group_by(
    Type, Bootstrap # "Study ID", but across controls
  ) %>% dplyr::group_modify(
    .f = function(x, y) {
      # "Pristine"/Control Data
      controlData <- x %>% dplyr::filter(Control == "Control")
      # Disturbed/Experimental Data.
      # disturbData <- x %>% dplyr::filter(Control != "Control")
      
      expand.grid(
        unique(paste(controlData$Patch, controlData$Control, sep = "-")),
        unique(paste(x$Patch, x$Control, sep = "-"))
      ) %>% dplyr::rename(
        s1 = Var1, s2 = Var2
      ) %>% dplyr::mutate(
        s1 = as.character(s1), s2 = as.character(s2),
        contrast = paste(s1, "vs", s2, sep = "_")
      ) %>% dplyr::filter(
        s1 != s2
      )
    }
  ) %>% tidyr::separate(
    col = s1,
    into = c("s1", "s1Control"),
    sep = "[-]", convert = TRUE
  ) %>% tidyr::separate(
    col = s2,
    into = c("s2", "s2Control"),
    sep = "[-]", convert = TRUE
  )
  
  ##### Jaccard Data: ###########################################################
  bootstrapSamplesBIIJaccard <- foreach::foreach(
    site = iterators::iter(bootstrapSamplesBIISites, by = "row"),
    .combine = dplyr::bind_rows
  ) %dopar% {
    data <- bootstrapSamplesBIISource %>% dplyr::filter(
      Type == site$Type, Bootstrap == site$Bootstrap
    )
    
    sor <- getJacAbSym(s1 = site$s1, s2 = site$s2,
                       s1Control = site$s1Control,
                       s2Control = site$s2Control,
                       data = data)
    
    site %>% dplyr::mutate(
      JacAbSym = sor
    )
  } %>% dplyr::ungroup()
  
  # As with abundance, something similar is happening here. Notably, time series
  # does NOT share the same behaviour as Space for time.
  # bootstrapSamplesBIIJaccard %>% ggplot2::ggplot(
  #   ggplot2::aes(x = interaction(s1Control, s2Control, Type), y = JacAbSym)
  # ) + ggplot2::geom_violin() + ggplot2::geom_boxplot(notch = TRUE, width = 0.1)
  # My guess is that the increase in "Time series" is driven by species that
  # appeared in the tail of the control period, which then stabilise through the
  # experiment period, giving the illusion of increasing intactness.
  # One way to test would be to decrease the rate of colonisations.
  
  ##### Bray-Curtis Data: #######################################################
  bootstrapSamplesBIIBray <- foreach::foreach(
    site = iterators::iter(bootstrapSamplesBIISites, by = "row"),
    .combine = dplyr::bind_rows
  ) %dopar% {
    data <- bootstrapSamplesBIISource %>% dplyr::filter(
      Type == site$Type, Bootstrap == site$Bootstrap
    )
    
    sor <- get_bray(s1 = site$s1, s2 = site$s2,
                    s1Control = site$s1Control,
                    s2Control = site$s2Control, data = data)
    
    site %>% dplyr::mutate(
      Bray = sor
    )
  } %>% dplyr::ungroup()
  
  # Different values (generally higher, but for highest?), but similar patterns.
  # bootstrapSamplesBIIBray %>% ggplot2::ggplot(
  #   ggplot2::aes(x = interaction(s1Control, s2Control, Type), y = Bray)
  # ) + ggplot2::geom_violin() + ggplot2::geom_boxplot(notch = TRUE, width = 0.1)
  
  ##### Convert to Predicted: ###################################################
  # This involves modelling each statistic (Abundance, Jaccard, Bray-Curtis),
  # followed by predicting for our patches based on the assumption that they are
  # all the pristine type instead of disturbed type (Control == "Control").
  # Finally, we'll be dividing the actual divided by the prediction.
  
  ####### Abundance: ############################################################
  # Note: Following plot is roughly linear for
  #       "MNA-ExampleExtProp-Result-Env10-Ring-Inf-1-1-ExtProp1.RData".
  #       This leads to a singular fit for Space for Time.
  # ggplot2::ggplot(
  #   bootstrapSamplesBIIAbundance %>% dplyr::arrange(
  #     RescaledAbundance
  #   ) %>% dplyr::mutate(
  #     number = 1:length(RescaledAbundance)
  #   ),
  #   ggplot2::aes(x = number, y = (RescaledAbundance)^2,
  #                color = interaction(Type, Control))
  # ) + ggplot2::geom_point(
  # ) + ggplot2::facet_wrap(
  #   Type ~ Control
  # )
  #
  # The oddity with the data is that we have fractions that include 1 (maybe 0).
  # Some looking around suggests Fractional Response Regression, as opposed to
  # logistic. (BII avoids this issue by ignoring the boundedness of abundance
  # and by adjusting the highest values to peak at 0.999 instead of 1.)
  # For now, at Jack and Tadhg's suggestion, we'll just do what they do.
  
  modelAbundanceSpaceTime <- lme4::lmer(
    sqrt(RescaledAbundance) ~ Control + (1 | Bootstrap),
    data = bootstrapSamplesBIIAbundance %>% dplyr::filter(
      Type == "Space for time"
    )
  )
  
  modelAbundanceTimeSeries <- lme4::lmer(
    sqrt(RescaledAbundance) ~ Control + (1 | Bootstrap),
    data = bootstrapSamplesBIIAbundance %>% dplyr::filter(
      Type == "Time series"
    )
  )
  
  predictedAbundanceSpaceTime <- data.frame(
    Control = unique(bootstrapSamplesBIIAbundance$Control),
    Type = "Space for time"
  ) %>% dplyr::mutate(
    RescaledAbundance = predict(modelAbundanceSpaceTime, ., re.form = NA)^2
  )
  
  predictedAbundanceTimeSeries <- data.frame(
    Control = unique(bootstrapSamplesBIIAbundance$Control),
    Type = "Time series"
  ) %>% dplyr::mutate(
    RescaledAbundance = predict(modelAbundanceTimeSeries, ., re.form = NA)^2
  )
  
  ####### Jaccard: ##############################################################
  bootstrapSamplesBIIJaccard <- bootstrapSamplesBIIJaccard %>% dplyr::mutate(
    logitCS = car::logit(JacAbSym, adjust = 0.001, percents = FALSE),
    contrastControl = factor(paste(s1Control, s2Control, sep = "-")),
    contrastControl = relevel(contrastControl, ref = "Control-Control")
  )
  
  modelJaccardSpaceTime <- lme4::lmer(
    logitCS ~ contrastControl + (1 | Bootstrap) + (1 | s2),
    data = bootstrapSamplesBIIJaccard %>% dplyr::filter(
      Type == "Space for time"
    )
  )
  
  modelJaccardTimeSeries <- lme4::lmer(
    logitCS ~ contrastControl + (1 | Bootstrap) + (1 | s2),
    data = bootstrapSamplesBIIJaccard %>% dplyr::filter(
      Type == "Time series"
    )
  )
  
  predictedJaccardSpaceTime <- data.frame(
    contrastControl = unique(bootstrapSamplesBIIJaccard$contrastControl),
    Type = "Space for time"
  ) %>% dplyr::mutate(
    JacAbSym = predict(
      modelJaccardSpaceTime, ., re.form = NA
    ) %>% inv_logit(a = 0.001)
  ) %>% tidyr::separate(
    col = contrastControl, remove = FALSE,
    into = c("s1Control", "Control"), sep = "[-]"
  )
  
  predictedJaccardTimeSeries <- data.frame(
    contrastControl = unique(bootstrapSamplesBIIJaccard$contrastControl),
    Type = "Time series"
  ) %>% dplyr::mutate(
    JacAbSym = predict(
      modelJaccardTimeSeries, ., re.form = NA
    ) %>% inv_logit(a = 0.001)
  ) %>% tidyr::separate(
    col = contrastControl, remove = FALSE,
    into = c("s1Control", "Control"), sep = "[-]"
  )
  
  ####### Bray-Curtis: ##########################################################
  bootstrapSamplesBIIBray <- bootstrapSamplesBIIBray %>% dplyr::mutate(
    logitCS = car::logit(Bray, adjust = 0.001, percents = FALSE),
    contrastControl = factor(paste(s1Control, s2Control, sep = "-")),
    contrastControl = relevel(contrastControl, ref = "Control-Control")
  )
  
  modelBraySpaceTime <- lme4::lmer(
    logitCS ~ contrastControl + (1 | Bootstrap) + (1 | s2),
    data = bootstrapSamplesBIIBray %>% dplyr::filter(
      Type == "Space for time"
    )
  )
  
  modelBrayTimeSeries <- lme4::lmer(
    logitCS ~ contrastControl + (1 | Bootstrap) + (1 | s2),
    data = bootstrapSamplesBIIBray %>% dplyr::filter(
      Type == "Time series"
    )
  )
  
  predictedBraySpaceTime <- data.frame(
    contrastControl = unique(bootstrapSamplesBIIBray$contrastControl),
    Type = "Space for time"
  ) %>% dplyr::mutate(
    Bray = predict(
      modelBraySpaceTime, ., re.form = NA
    ) %>% inv_logit(a = 0.001)
  ) %>% tidyr::separate(
    col = contrastControl, remove = FALSE,
    into = c("s1Control", "Control"), sep = "[-]"
  )
  
  predictedBrayTimeSeries <- data.frame(
    contrastControl = unique(bootstrapSamplesBIIBray$contrastControl),
    Type = "Time series"
  ) %>% dplyr::mutate(
    Bray = predict(
      modelBrayTimeSeries, ., re.form = NA
    ) %>% inv_logit(a = 0.001)
  ) %>% tidyr::separate(
    col = contrastControl, remove = FALSE,
    into = c("s1Control", "Control"), sep = "[-]"
  )
  
  ##### Final BII Predictions: ##################################################
  # Usage: Mutiply BII by proportion that is each land use type ("Control").
  predictedBII <- dplyr::bind_rows(
    predictedAbundanceSpaceTime,
    predictedAbundanceTimeSeries
  ) %>% dplyr::full_join(
    y = dplyr::bind_rows(
      predictedBraySpaceTime,
      predictedBrayTimeSeries
    ), by = c("Control", "Type")
  ) %>% dplyr::mutate(
    BII = RescaledAbundance * Bray
  )
  
  ### Scholes and Biggs (2005) BII: #############################################
  # Need to interpret:
  #   Taxon : Suggests group by ID. Could be explored more thoroughly with a
  #           phylogenetic model or functional niche.
  #   Ecosystem : Either a patch or a set of patches. Guess not important here.
  #   Land Use : Suggests type of disturbance, which we'll take to be Control.
  #   Richness of Taxon : We'll take this to be 1 species per ID.
  #   Area of Land Use in Ecosystem : If patch, 1, if set, number of patches.
  #   Abundance divided by Reference Abundance : Complicated:
  #      Note that they exclude "alien" species from the calculation, so
  #      we don't have to worry about reference abundance of 0.
  #      De Palma et al. would argue we should regress for reference abundance
  #      to help control errors in sampling. This seems reasonable.
  #      Unlike De Palma though, I think this modelling should be done at
  #      the taxon level.
  #      One major problem though is the autocorrelation; another is that
  #      multiple samples from the same place can vary greatly.
  
  # Note that variance != mean here upon initial check
  # %>% dplyr::group_by(ID) %>% dplyr::summarise(
  #        Mean = mean(Abundance),
  #        Var = var(Abundance)
  #    )
  if(calculationsLong) {
    modelSBAbundanceSpaceTime <- lme4::glmer.nb(
      Abundance ~ ID + (1 | Bootstrap) + (1 | Patch),
      data = bootstrapSamplesBIISource %>% dplyr::filter(
        Invasive == FALSE, Control == "Control", Type == "Space for time"
      ) %>% dplyr::mutate(
        ID = factor(ID)
      ) # samples in Abundance notation
    )
    
    modelSBAbundanceTimeSeries <-  lme4::glmer.nb(
      Abundance ~ ID + (1 | Bootstrap) + (1 | Patch),
      data = bootstrapSamplesBIISource %>% dplyr::filter(
        Invasive == FALSE, Control == "Control", Type == "Time series"
      ) %>% dplyr::mutate(
        ID = factor(ID)
      ) # samples in Abundance notation
    )
    
    predictSBAbundanceSpaceTimeDat <- bootstrapSamplesBIISource %>% dplyr::filter(
      Invasive == FALSE, Control == "Experiment",
      Type == "Space for time"
    ) %>% dplyr::mutate(
      ID = factor(ID)
    ) %>% dplyr::distinct(ID)
    
    predictSBAbundanceTimeSeriesDat <- bootstrapSamplesBIISource %>% dplyr::filter(
      Invasive == FALSE, Control == "Experiment",
      Type == "Time series"
    ) %>% dplyr::mutate(
      ID = factor(ID)
    ) %>% dplyr::distinct(ID)
    
    predictSBAbundanceSpaceTime <-
      predict(modelSBAbundanceSpaceTime,
              predictSBAbundanceSpaceTimeDat, re.form = NA, type = "response")
    
    predictSBAbundanceTimeSeries <-
      predict(modelSBAbundanceTimeSeries,
              predictSBAbundanceTimeSeriesDat, re.form = NA, type = "response")
    
    bootstrapSamplesBIISource <- bootstrapSamplesBIISource %>% dplyr::mutate(
      Invasive = dplyr::case_when(
        Invasive == "TRUE" ~ TRUE,
        Invasive == "FALSE" ~ FALSE
      )
    ) %>% dplyr::left_join(
      rbind(
        data.frame(
          ID = as.numeric(as.character(predictSBAbundanceSpaceTimeDat$ID)),
          Invasive = FALSE,
          Control = "Experiment",
          Type = "Space for time",
          AbundanceRef = predictSBAbundanceSpaceTime
        ),
        data.frame(
          ID = as.numeric(as.character(predictSBAbundanceTimeSeriesDat$ID)),
          Invasive = FALSE,
          Control = "Experiment",
          Type = "Time series",
          AbundanceRef = predictSBAbundanceTimeSeries
        )
      ),
      by = c("ID", "Invasive", "Control", "Type")
    )
    
    # Since each Taxon has one Species and we're comparing Control to Experiment,
    # We consider R(ij) = 1 and A(jk) = 5 regardles of i, j, k, so they cancel.
    # Then it is simply sum(native Abundace / native Reference Abundance)/num. nat.
    
    # I'm not perfectly convinced by how I've done it though. Going with gut here
    # and then we can refine it later.
    
    predictedSB <- bootstrapSamplesBIISource %>% dplyr::filter(
      !is.na(AbundanceRef)
    ) %>% dplyr::group_by(Type, Bootstrap) %>% dplyr::summarise(
      BII = mean(Abundance/AbundanceRef)
    )
  }
  
  ### Bray-Curtis's: ############################################################
  
  bootstrapSamplesBrayAll <- foreach::foreach(
    site = iterators::iter(bootstrapSamplesBIISites, by = "row"),
    .combine = dplyr::bind_rows
  ) %dopar% {
    data <- bootstrapSamplesBIISource %>% dplyr::filter(
      Type == site$Type, Bootstrap == site$Bootstrap
    )
    
    sor <- get_bray_all(s1 = site$s1, s2 = site$s2,
                        s1Control = site$s1Control,
                        s2Control = site$s2Control, data = data)
    
    cbind(site, sor)
  } %>% dplyr::ungroup()
}

# Cleanup: ####################################################################
parallel::stopCluster(clust)

# Plotting: ###################################################################
### Plot 0: Sense Checking: ###################################################
# Richness Comparison
plot_0_Richness <- ggplot2::ggplot(
  bootstrapSamples %>% dplyr::mutate(
    TrueRichness = unlist(lapply(
      strsplit(SamplingNonZeroSpecies, split = ", ", fixed = TRUE), length
    ))
  ),
  ggplot2::aes(x = Time, y = TrueRichness, color = Patch, group = Patch)
) + ggplot2::geom_line(
) + ggplot2::geom_line(
  ggplot2::aes(y = SamplingAlpha * 3),
  alpha = 0.25,
  linetype = 2
) + ggplot2::scale_y_continuous(
  sec.axis = ggplot2::sec_axis(~ . / 3, name = "Richness Observed")
) + ggplot2::facet_wrap(. ~ Patch)

# But correlation remains low? (but significant) (For file 5-1-1)
# with(bootstrapSamples %>% dplyr::mutate(
#   TrueRichness = unlist(lapply(
#     strsplit(SamplingNonZeroSpecies, split = ", ", fixed = TRUE), length
#   ))
# ), cor.test(TrueRichness, SamplingAlpha))

plot_0_Abundance <- ggplot2::ggplot(
  bootstrapSamples,
  ggplot2::aes(x = Time, y = SamplingAbundance, color = Patch, group = Patch)
) + ggplot2::geom_line(
) + ggplot2::geom_line(
  ggplot2::aes(y = SamplingObserved * 100),
  alpha = 0.25,
  linetype = 2
) + ggplot2::scale_y_continuous(
  name = "True Total Abundance",
  sec.axis = ggplot2::sec_axis(~ . / 100, name = "Abundance Observed")
) + ggplot2::facet_wrap(. ~ Patch)

# Abundance is (unsurprisingly) much more correlated. (For file 5-1-1)
# with(bootstrapSamples, cor.test(SamplingAbundance, SamplingObserved))

if (calculationsPlotLong) {
  plot_0_AbundanceSpecies <- lapply(
    1:result$NumEnvironments,
    function(i) RMTRCode2::LawMorton1996_PlotAbundance(
      result$Abundance[
        result$Abundance[, 1] > min(bootstrapSamples$Time) &
          result$Abundance[, 1] < max(bootstrapSamples$Time)
        , c(1, (i-1) * 200 + 1:200 + 1)], guides = FALSE
    ) + ggplot2::scale_y_log10()
  ) %>% patchwork::wrap_plots()
}

# bootstrapSamplesOccupancy <- bootstrapSamples %>% dplyr::group_by(
#   Bootstrap, Time, Type, Control, TimeActualRow, TimeActual, TimeSinceStart
# ) %>% dplyr::group_map(
#   function(x, y) {
#     # Per Bootstrap, Time, Patch Set (Control/Experiment)
#     # We are reproducing our Figure 3a in principle.
#     # Note that that would mean No. of Bootstrap plots!
#     # It doesn't make sense to average these on a per species basis.
#     # It might make sense, however, to average on the dominant niche
#     # (according to the pool) though.
#     # How do we minimise the number of plots? We can reduce to 100 + 2
#     # by having Fig 3a equivalent have nowhere present as NA, and present
#     # as -5 : 5 (all Control : all Experiment) and the 2 as richness plots.
#
#     # Second thoughts: maybe we only need the true occupancy.
#
#     TruePresent <- unique(unlist(
#       strsplit(x$SamplingNonZeroSpecies, split = ", ")
#       ))
#
#   }
# )

plot_0_OccupancyTrue <- (
  speciesAbundances <- RMTRCode2::Calculate_Species(
    result
  )
) %>% dplyr::group_by(
  Time, Species
) %>% dplyr::summarise(
  Count = dplyr::n()
) %>% dplyr::left_join(
  Pool %>% dplyr::arrange(Size) %>% dplyr::mutate(SizeID = 1:nrow(Pool)),
  by = c("Species" = "ID")
) %>% ggplot2::ggplot(
  ggplot2::aes(x = Time, y = SizeID, color = Count)
) + ggplot2::geom_point(
  shape = '.'
) + ggplot2::scale_color_viridis_c(
  direction = -1, limits = c(1, 10)
) + ggplot2::geom_hline(
  yintercept = nrow(Pool %>% dplyr::filter(Type == Type[1])) + 0.5,
  color = "red"
) + ggplot2::labs(
  y = "Species by Size"
) + ggplot2::theme_bw()


(speciesAbundanceStats <- speciesAbundances %>% dplyr::group_by(
  Time, Species
) %>% dplyr::summarise(
  Count = dplyr::n(),
  AvgAbundance = mean(Abundance),
  AvgLogAbundance = mean(log(Abundance))
))

plot_0_AbundanceAvg <- speciesAbundanceStats %>% ggplot2::ggplot(
  ggplot2::aes(x = Time, y = AvgAbundance,
               color = as.character(Species), group = Species)
) + ggplot2::geom_line(
) + ggplot2::labs(
  y = "Average Abundance"
) + ggplot2::facet_wrap(
  . ~ Count
) + ggplot2::theme_bw() + ggplot2::scale_y_log10(
) + ggplot2::guides(color = "none")

plot_0_AbundanceGeoAvg <- speciesAbundanceStats %>% ggplot2::ggplot(
  ggplot2::aes(x = Time, y = exp(AvgLogAbundance),
               color = as.character(Species), group = Species)
) + ggplot2::geom_line(
) + ggplot2::labs(
  y = "Average Abundance (Geometric)"
) + ggplot2::facet_wrap(
  . ~ Count
) + ggplot2::theme_bw() + ggplot2::scale_y_log10(
) + ggplot2::guides(color = "none")


# Still need to do Observed Species Abundances, consumer richness,
# basal richness, and to add correlations to the plots.
#
# plot_0_Consumers <-

### Plot 1: Change in Local Richness: #########################################
plot_1_DeltaAlpha <- ggplot2::ggplot(
  bootstrapSamplesDeltaAlpha,
  ggplot2::aes(
    x = Type,
    y = `Difference of Average Number of Species in Patch`
  )
) + ggplot2::geom_violin(
) + ggplot2::geom_boxplot(
  width = 0.1,
  notch = TRUE
) + ggplot2::geom_line(
  ggplot2::aes(group = Bootstrap),
  alpha = 0.1
) + ggplot2::facet_wrap(
  `Species Subset` ~ .
) + ggplot2::labs(
  title = "Delta Alpha",
  subtitle = "Difference is Experiment - Control",
  caption = paste0("file: ", file_result)
)

# Note: Can try to establish if there is correlation here.
# with(
#   bootstrapSamplesDeltaAlpha %>% dplyr::filter(
#     `Species Subset` == "Overall"
#   ) %>% tidyr::pivot_wider(
#     names_from = "Type",
#     values_from = "Difference of Average Number of Species in Patch",
#     id_cols = "Bootstrap"
#   ), cor.test(`Space for time`, `Time series`)
# )

### Plot 2: ###################################################################
plot_2_DeltaGamma <- ggplot2::ggplot(
  bootstrapSamplesDeltaGamma,
  ggplot2::aes(
    x = Type,
    y = `Difference of Number of Species across Patches`
  )
) + ggplot2::geom_violin(
) + ggplot2::geom_boxplot(
  width = .1,
  notch = TRUE
) + ggplot2::geom_line(
  ggplot2::aes(group = Bootstrap),
  alpha = 0.1
) + ggplot2::facet_wrap(
  `Species Subset` ~ .
) + ggplot2::labs(
  title = paste0("Delta Gamma"),
  subtitle = "Difference is Experiment - Control",
  caption = paste0("file: ", file_result)
)

### Plot 3: ###################################################################
if (exists("predictedBII")) {
  plot_3_dodge <- ggplot2::position_dodge(0.8)
  plot_3_BIITotal <- ggplot2::ggplot(
    predictedBII,
    ggplot2::aes(x = Type, y = BII, color = Control)
  ) + ggplot2::geom_point(
    position = plot_3_dodge,
    shape = 4, size = 4
  ) + ggplot2::labs(
    title = "BII Weights"#,
    # subtitle = "Crosses are model predictions.",
    # caption = paste0("file: ", file_result)
  )
  
  plot_3_BIIAbundance <- ggplot2::ggplot(
    bootstrapSamplesBIIAbundance,
    ggplot2::aes(x = Type, y = RescaledAbundance, fill = Control)
  ) + ggplot2::geom_violin(
    position = plot_3_dodge
  ) + ggplot2::geom_boxplot(
    position = plot_3_dodge,
    width = .1,
    notch = TRUE
  ) + ggplot2::geom_point(
    data = dplyr::bind_rows(predictedAbundanceSpaceTime,
                            predictedAbundanceTimeSeries),
    position = plot_3_dodge,
    shape = 4
  ) + ggplot2::labs(
    title = "Rescaled Abundance"#,
    # subtitle = "Crosses are model predictions.",
    # caption = paste0("file: ", file_result)
  )
  
  plot_3_BIIJaccard <- ggplot2::ggplot(
    bootstrapSamplesBIIJaccard,
    ggplot2::aes(x = Type, y = JacAbSym, fill = contrastControl)
  ) + ggplot2::geom_violin(
    position = plot_3_dodge
  ) + ggplot2::geom_boxplot(
    position = plot_3_dodge,
    width = .1,
    notch = TRUE
  ) + ggplot2::geom_point(
    data = dplyr::bind_rows(predictedJaccardSpaceTime,
                            predictedJaccardTimeSeries),
    position = plot_3_dodge,
    shape = 4
  ) + ggplot2::labs(
    title = "Asymmetric Jaccard Index"#,
    # subtitle = "Crosses are model predictions.",
    # caption = paste0("file: ", file_result)
  )
  
  plot_3_BIIBray <- ggplot2::ggplot(
    bootstrapSamplesBIIBray,
    ggplot2::aes(x = Type, y = Bray, fill = contrastControl)
  ) + ggplot2::geom_violin(
    position = plot_3_dodge
  ) + ggplot2::geom_boxplot(
    position = plot_3_dodge,
    width = .1,
    notch = TRUE
  ) + ggplot2::geom_point(
    data = dplyr::bind_rows(predictedBraySpaceTime,
                            predictedBrayTimeSeries),
    position = plot_3_dodge,
    shape = 4
  ) + ggplot2::labs(
    title = "Balanced Bray-Curtis"#,
    # subtitle = "Crosses are model predictions.",
    # caption = paste0("file: ", file_result)
  )
  
  plot_3_BII <- (
    (plot_3_BIITotal + plot_3_BIIAbundance)/(plot_3_BIIBray + plot_3_BIIJaccard)
  ) + patchwork::plot_layout(
    guides = "collect"
  ) + patchwork::plot_annotation(
    title = "Biodiversity Intactness Index Breakdown",
    subtitle = "Crosses are model predictions.",
    caption = paste0("file: ", file_result)
  )
}

### Plot 4: ###################################################################
if (exists("predictedSB")) {
  plot_4_SB <- ggplot2::ggplot(
    predictedSB,
    ggplot2::aes(x = Type, y = BII)
  ) + ggplot2::geom_violin(
    fill = "#00BFC4"
  ) + ggplot2::geom_boxplot(
    width = .1,
    notch = TRUE,
    fill = "#00BFC4"
  ) + ggplot2::geom_hline(
    yintercept = 1, color = "#F8766D"
  ) + ggplot2::labs(
    title = "Scholes and Biggs: BII",
    subtitle = "mean(Abundance/Modelled Abundance)",
    caption = paste0("file: ", file_result)
  )
}

### Plot 5: ###################################################################
if (exists("bootstrapSamplesBrayAll")) {
  plot_5_Bray <- ggplot2::ggplot(
    bootstrapSamplesBrayAll %>% tidyr::pivot_longer(
      cols = bray.bal : bray,
      names_to = "Component",
      values_to = "Dissimilarity"
    ) %>% dplyr::mutate(
      contrastControl = s2Control
    ),
    ggplot2::aes(x = Type, y = Dissimilarity, fill = contrastControl)
  ) + ggplot2::geom_violin(
    position = plot_3_dodge
  ) + ggplot2::geom_boxplot(
    position = plot_3_dodge,
    width = .1,
    notch = TRUE
    # ) + ggplot2::geom_point( # Not implemented predictions yet.
    #   data = dplyr::bind_rows(
    #     predictedBraySpaceTime,
    #     predictedBrayTimeSeries
    #   ) %>% dplyr::rename(
    #     Dissimilarity = Bray
    #   ) %>% dplyr::mutate(
    #     contrastControl = gsub("Control-", "", contrastControl, fixed = TRUE),
    #     Component = "bray"
    #   ),
    #   position = plot_3_dodge,
    #   shape = 4
  ) + ggplot2::labs(
    title = "Balanced Bray-Curtis"#,
    # subtitle = "Crosses are model predictions.",
    # caption = paste0("file: ", file_result)
  ) + ggplot2::facet_wrap(
    . ~ Component
  )
}

### Plot 6 (Alpha Slope): #####################################################
# This plot is a little too on the nose for what we're trying to convey.
# We need to step back a bit.
# ggplot2::ggplot(
#   bootstrapSamples %>% tidyr::pivot_longer(
#     cols = c(SamplingAlpha, SamplingAlphaNative, SamplingAlphaInvasive),
#     names_to = "Subset", values_to = "Patch Species Richness"
#   ) %>% dplyr::mutate(
#     Subset = dplyr::case_when(
#       Subset == "SamplingAlpha" ~ "All",
#       Subset == "SamplingAlphaNative" ~ "Native",
#       Subset == "SamplingAlphaInvasive" ~ "Invasive"
#     )),
#   ggplot2::aes(
#     x = TimeSinceStart,
#     y = `Patch Species Richness`
#   )
# ) + ggplot2::geom_point(
#   alpha = 0.2
# ) + ggplot2::facet_grid(
#   Subset ~ Type
# ) + ggplot2::geom_abline(
#   data = bootstrapSamplesLMAlpha,
#   mapping = ggplot2::aes(
#     slope = Slope,
#     intercept = Intercept,
#     group = interaction(Bootstrap, Control, Subset)
#   ),
#   alpha = 0.2
# )
plot_6_LMSlope <- ggplot2::ggplot(
  bootstrapSamplesLMAlpha,
  ggplot2::aes(
    x = Type,
    y = Slope
  )
) + ggplot2::geom_violin(
) + ggplot2::geom_boxplot(
  width = 0.1,
  notch = TRUE
) + ggplot2::geom_line(
  ggplot2::aes(group = Bootstrap),
  alpha = 0.1
) + ggplot2::facet_wrap(
  Subset ~ .
) + ggplot2::labs(
  title = "Linear Fitted Slopes",
  subtitle = "Patch Alpha ~ Time Since Recording : Control + (1 | Patch)"#,
  #caption = paste0("file: ", file_result)
)

plot_6_LMIntercept <- ggplot2::ggplot(
  bootstrapSamplesLMAlpha,
  ggplot2::aes(
    x = Type,
    y = Intercept
  )
) + ggplot2::geom_violin(
) + ggplot2::geom_boxplot(
  width = 0.1,
  notch = TRUE
) + ggplot2::geom_line(
  ggplot2::aes(group = Bootstrap),
  alpha = 0.1
) + ggplot2::facet_wrap(
  Subset ~ .
) + ggplot2::labs(
  title = "Linear Fitted Intercepts",
  #subtitle = "Patch Alpha ~ Time Since Recording : Control + (1 | Patch)",
  caption = paste0("file: ", file_result)
)

plot_6_LM <- plot_6_LMSlope / plot_6_LMIntercept
