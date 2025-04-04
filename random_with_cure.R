library(purrr) # install.packages("purrr")
library(rlang) # install.packages("rlang")
library(tibble) # install.packages("tibble")
library(ggplot2) # install.packages("ggplot2")
library(cowplot) # install.packages("cowplot")
library(survival) # install.packages("survival")
library(glue) # install.pakcages("glue")

rm(list = ls(all = TRUE))

# Funcao que usaremos para gerar os dados
random_cure <-
  function(n = 1L,
           cure_fraction,
           quantile_function,
           cdf = NULL,
           t_max = NULL,
           uniform_censoring = TRUE,
           lower = 0,
           upper = 50,
           args_model) {
    find_root <- function(u, ...) {
      uniroot(\(x) cdf(x, ...) - u, lower = lower, upper = upper, tol = .Machine$double.eps^0.5)$root
    }
    find_root <- Vectorize(FUN = find_root, vectorize.args = c("u"))

    c <- rbinom(n = n, size = 1L, prob = 1 - cure_fraction)
    u <- runif(n = n, min = 0, max = 1 - cure_fraction)

    if (is.function(quantile_function)) {
      quantile_function <- purrr::partial(.f = quantile_function, !!!args_model)
    } else if (is.function(cdf)) {
      quantile_function <- purrr::partial(.f = find_root, !!!args_model)
    } else {
      stop("Invalid input: You must provide either a valid quantile function (quantile_function) or a cumulative distribution function (cdf).")
    }

    t <- quantile_function(u)
    t_c <- ifelse(c, yes = t, no = .Machine$double.xmax)

    if (!is.null(t_max)) {
      t_observed <- pmin(t_c, t_max)
      delta <- ifelse(t_observed < t_max, yes = 1L, no = 0L)
    } else {
      t_max <- max(t * c)
      if (uniform_censoring) {
        w <- runif(n = n, min = 0, max = t_max)
      } else {
        w <- rexp(n = n, rate = t_max)
      }
      t_observed <- pmin(t_c, w)
      delta <- ifelse(t_observed < w, yes = 1L, no = 0L)
    }
    return(tibble::tibble(t = t_observed, delta = delta))
  }

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