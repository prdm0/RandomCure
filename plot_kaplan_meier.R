library(ggplot2)
library(survival)

# Função usada para construir o Kaplan-Meier mais fácil no ggplot2
plot_kaplan <- function(
  dados,
  surv,
  title = "Kaplan-Meier Survival Curve",
  xlab = "Time",
  ylab = "Survival Probability",
  color = "#1F77B4",
  ci_color = "#1F77B4",
  cure_color = "#D62728", # Vivid red
  base_size = 11,
  ci_alpha = 0.2,
  line_size = 0.7,
  show_cure = TRUE,
  cure_label = "Cure fraction: ",
  ...
) {
  # Calculate cure probability
  cure_prob <- surv(.Machine$double.xmax, ...)

  # Fit Kaplan-Meier model
  dados_surv <- Surv(time = dados$t, event = dados$delta)
  km_fit <- survfit(dados_surv ~ 1)

  # Prepare plot data
  km_data <- data.frame(
    time = km_fit$time,
    surv = km_fit$surv,
    upper = km_fit$upper,
    lower = km_fit$lower
  )

  max_time <- max(km_data$time)

  # Create legend data
  legend_data <- data.frame(
    label = c("Kaplan-Meier", "Cure fraction"),
    color = c(color, cure_color),
    linetype = c("solid", "solid")
  )

  # Create base plot
  p <- ggplot(km_data, aes(x = time)) +
    # Confidence interval
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = ci_color,
      alpha = ci_alpha
    ) +
    # Main survival curve
    geom_step(aes(y = surv), color = color, linewidth = line_size) +
    # Cure line
    {
      if (show_cure)
        annotate(
          "segment",
          x = 0,
          xend = max_time,
          y = cure_prob,
          yend = cure_prob,
          color = cure_color,
          linewidth = line_size,
          linetype = "solid"
        )
    } +
    # Scales
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      expand = c(0, 0.02)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    # Legend
    geom_segment(
      data = legend_data,
      aes(
        x = max_time * 0.7,
        xend = max_time * 0.75,
        y = 0.95 - (0:1) * 0.08,
        yend = 0.95 - (0:1) * 0.08,
        color = color,
        linetype = linetype
      ),
      linewidth = line_size
    ) +
    geom_text(
      data = legend_data,
      aes(x = max_time * 0.78, y = 0.95 - (0:1) * 0.08, label = label),
      hjust = 0,
      size = base_size / 3
    ) +
    scale_color_identity() +
    scale_linetype_identity() +
    # Labels
    labs(
      x = xlab,
      y = ylab,
      title = title,
      subtitle = if (show_cure)
        paste0(cure_label, format(round(cure_prob, 3), nsmall = 3)) else NULL
    ) +
    # Theme
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        linewidth = 0.3,
        color = "#D3D3D3",
        linetype = "dashed"
      ),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      plot.title = element_text(
        hjust = 0,
        face = "bold",
        size = base_size + 2,
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(
        hjust = 0,
        face = "bold",
        color = "black",
        size = base_size,
        margin = margin(b = 10)
      ),
      legend.position = "none",
      plot.margin = margin(10, 15, 10, 10),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )

  return(p)
}
