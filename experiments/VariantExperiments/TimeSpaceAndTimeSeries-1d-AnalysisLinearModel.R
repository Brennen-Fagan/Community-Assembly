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