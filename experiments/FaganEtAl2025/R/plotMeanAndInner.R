plotMeanAndInner <- function(
  data, CIs = c(0.5, 0.95), digits = -2,
  facets = as.formula(Intervention ~ SpeciesPreferences)
) {
  # Correct for a problem with how to handle NAs by converting to strings
  data$Subset <- ifelse(is.na(data$Subset), "NA", data$Subset)

  # Create base with particular attention to grouping structure.
  baseplot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = Time, y = Value,
      group = interaction(
        Subset,
        Intervention, InterventionInitial, InterventionFinal, SpeciesPreferences
      ),
      color = Intervention, fill = Intervention, linetype = SpeciesPreferences
    )
  )

  # Plot each CI overlaid. Overlaying => the innermost have the darkest alpha.
  for (CI in CIs) {
    baseplot <- baseplot + ggplot2::geom_ribbon(
      data = data |> tidytable::mutate(
        Time = round(Time, digits = digits)
      ) |> tidytable::group_by(
        Time, Subset,
        Intervention, InterventionInitial, InterventionFinal, SpeciesPreferences
      ) |> tidytable::summarise(
        top = quantile(Value, probs = CI+(1-CI)/2, na.rm = TRUE),
        bot = quantile(Value, probs = (1-CI)-(1-CI)/2, na.rm = TRUE)
      ), mapping = ggplot2::aes(
        x = Time, ymin = bot, ymax = top,
        group = interaction(
          Subset,
          Intervention, InterventionInitial, InterventionFinal, SpeciesPreferences
        ),
        fill = Intervention,
        color = Intervention
      ), inherit.aes = FALSE,
      alpha = 0.25, linewidth = 0.25 #, linetype = "dotted"
    )
  }

  # Add an average line and handle the meta-details.
  baseplot <- baseplot + ggplot2::geom_line(
    data = data |> tidytable::mutate(
      Time = round(Time, digits = digits)
    ) |> tidytable::group_by(
      Time, Subset,
      Intervention, InterventionInitial, InterventionFinal, SpeciesPreferences
    ) |> tidytable::summarise(
      Value = mean(Value, na.rm = TRUE)
    )
  ) + ggplot2::facet_grid(
    facets
  ) + ggplot2::scale_color_manual(
    values = colorPalette, aesthetics = c("color", "fill"),
    name = "Habitat's Land-use"
  ) + ggplot2::scale_linetype_manual(
    name = "Species Preferences",
    values = linetypePalette
  ) + ggplot2::theme_minimal(
  ) + ggplot2::guides(
    color = if (length(CIs)>0) {"none"} else {"legend"},
    fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
  )

  return(baseplot)
}
