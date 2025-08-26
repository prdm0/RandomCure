plot_kaplan <- function(
  dados,
  surv,
  title = "Kaplan–Meier Survival Curve",
  xlab = "Time",
  ylab = "Survival Probability",
  color = "#1F77B4",
  ci_color = "#1F77B4",
  cure_color = "#D62728",
  base_size = 11,
  ci_alpha = 0.2,
  line_size = 0.7,
  show_cure = TRUE,
  show_percent = TRUE,
  cure_label = "Theoretical cure (π): ",
  extra_label = "Extra censoring among susceptibles: ",
  total_label = "Overall no-event proportion (δ=0): ",
  ...
) {
  # Cure fraction
  cure_prob <- tryCatch(
    do.call(surv, c(list(.Machine$double.xmax), list(...))),
    error = function(e) surv(.Machine$double.xmax, ...)
  )

  # Proportions
  overall_zeros <- mean(dados$delta == 0)
  extra_censor <- overall_zeros - cure_prob
  extra_censor <- max(0, extra_censor)

  fmt <- function(x) {
    if (show_percent) sprintf("%.1f%%", 100 * x) else sprintf("%.3f", x)
  }

  # KM fit
  dados_surv <- Surv(time = dados$t, event = dados$delta)
  km_fit <- survfit(dados_surv ~ 1)
  km_data <- data.frame(
    time = km_fit$time,
    surv = km_fit$surv,
    upper = km_fit$upper,
    lower = km_fit$lower
  )

  if (nrow(km_data) > 0 && tail(km_data$surv, 1) < cure_prob) {
    km_data <- km_data[-nrow(km_data), ]
  }
  max_time <- if (nrow(km_data)) max(km_data$time) else 0

  # Subtitle: multiline
  plot_subtitle <- if (show_cure) {
    paste0(
      cure_label,
      fmt(cure_prob),
      "\n",
      extra_label,
      fmt(extra_censor),
      "\n",
      total_label,
      fmt(overall_zeros)
    )
  } else {
    NULL
  }

  # Plot
  p <- ggplot(km_data, aes(x = time)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = ci_color,
      alpha = ci_alpha
    ) +
    geom_step(aes(y = surv), color = color, linewidth = line_size) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0.02)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = xlab, y = ylab, title = title, subtitle = plot_subtitle) +
    theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2),
      plot.subtitle = element_text(
        face = "bold",
        hjust = 0,
        color = "black",
        size = base_size,
        margin = margin(b = 10)
      ),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    )

  if (show_cure && is.finite(max_time) && max_time > 0) {
    p <- p +
      annotate(
        "segment",
        x = 0,
        xend = max_time,
        y = cure_prob,
        yend = cure_prob,
        color = cure_color,
        linewidth = line_size
      )
  }
  return(p)
}
