hillWrapper <- function(env, i) {
  if (length(dim(env)) >= 1 && ncol(env) == 1) {
    metrics <- do.call(rbind, lapply(env, function(i)
      if(i > 0) {
        c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
      } else {
        c(NA, NaN, NaN, NA, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
      }
    ))
    metrics <- data.frame(metrics)
    colnames(metrics) <-  c("0", "0.25", "0.5", "1", "2", "4",
                            "8", "16", "32", "64", "Inf")
  } else {
    metrics <- vegan::renyi(env, hill = TRUE) # returns data.frame, but only if
    if (length(dim(env)) <= 1 || nrow(env) == 1) {# nrow > 1.
      tempnames <- names(metrics)
      metrics <- data.frame(matrix(metrics, nrow = 1))
      colnames(metrics) <- tempnames
    }
  }
  names(metrics) <- paste0("Hill:", names(metrics))
  cbind(Environment1 = i,
        Environment2 = NA,
        metrics,
        stringsAsFactors = FALSE) %>% tidyr::pivot_longer(
          cols = tidyr::contains("Hill"),
          names_to = "Metric", values_to = "Value"
        )
}
