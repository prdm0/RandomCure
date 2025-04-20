library(ggplot2)
library(survival)

plot_kaplan <- function(
  dados,
  surv,
  title = "Kaplan-Meier Survival Curve",
  xlab = "Time",
  ylab = "Survival Probability",
  color = "#1F77B4",
  ci_color = "#1F77B4",
  cure_color = "#D62728",
  base_size = 11,
  ci_alpha = 0.2,
  line_size = 0.7,
  show_cure = TRUE,
  cure_label = "Cure: ",
  censoring_label = "Censoring: ",
  ...
) {
  # Teórica fração de cura
  cure_prob <- surv(.Machine$double.xmax, ...)

  # Proporções de censura
  total_censored <- mean(dados$delta == 0)
  additional_censoring <- total_censored - cure_prob

  # Ajuste Kaplan-Meier
  dados_surv <- Surv(time = dados$t, event = dados$delta)
  km_fit <- survfit(dados_surv ~ 1)

  # Dados da curva
  km_data <- data.frame(
    time = km_fit$time,
    surv = km_fit$surv,
    upper = km_fit$upper,
    lower = km_fit$lower
  )

  # Remove queda final artificial, se necessário
  if (tail(km_data$surv, 1) < cure_prob) {
    km_data <- km_data[-nrow(km_data), ]
  }

  max_time <- max(km_data$time)

  # Subtítulo informativo
  plot_subtitle <- if (show_cure) {
    bquote(
      .(cure_label) * .(sprintf("%.3f", abs(cure_prob))) ~ "+" ~
        .(censoring_label) * .(sprintf("%.3f", abs(additional_censoring))) ~
        "=" ~
        .(sprintf("%.3f", total_censored))
    )
  }

  # Construção do gráfico
  p <- ggplot(km_data, aes(x = time)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = ci_color,
      alpha = ci_alpha
    ) +
    geom_step(aes(y = surv), color = color, linewidth = line_size) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0.02)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      x = xlab,
      y = ylab,
      title = title,
      subtitle = plot_subtitle
    ) +
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

  # Adiciona linha da fração de cura com annotate()
  if (show_cure) {
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
