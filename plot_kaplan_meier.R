library(ggplot2)

# Função usada para construir o Kaplan-Meier mais fácil no ggplot2
plot_kaplan <- function(dados, cura_real = NULL) {
  dados_surv <- Surv(time = dados$t, event = dados$delta)
  km_fit <- survfit(dados_surv ~ 1)

  km_data <- data.frame(
    time = km_fit$time,
    surv = km_fit$surv,
    upper = km_fit$upper,
    lower = km_fit$lower
  )

  ggplot(km_data, aes(x = time, y = surv)) +
    geom_step() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    labs(x = "Tempo", y = "Probabilidade de Sobrevivência", title = "Curva de Kaplan-Meier") +
    stat_summary(geom = "line", fun = \(i) cura_real, aes(color = "Cura real"), linewidth = 1.5, linetype = 2) +
    scale_y_continuous(limits = c(0, 1), n.breaks = 10) +
    labs(x = "t") +
    guides(color = guide_legend(title = ""))
}