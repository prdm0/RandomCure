library(purrr)
library(rlang)
library(tibble)
library(survival)
library(ggplot2)

# Função usada para construir o Kaplan-Meier mais fácil no ggplot2
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
  cure_label = "Cure fraction: ",
  ...
) {
  require(ggplot2)
  require(survival)

  # Fração de cura teórica
  cure_prob <- surv(.Machine$double.xmax, ...)

  # Kaplan-Meier
  dados_surv <- Surv(time = dados$t, event = dados$delta)
  km_fit <- survfit(dados_surv ~ 1)

  # Dados da curva
  km_data <- data.frame(
    time = km_fit$time,
    surv = km_fit$surv,
    upper = km_fit$upper,
    lower = km_fit$lower
  )

  # Remover o último ponto se ele for menor que a fração de cura
  if (tail(km_data$surv, 1) < cure_prob) {
    km_data <- km_data[-nrow(km_data), ]
  }

  max_time <- max(km_data$time)

  # Legenda
  legend_data <- data.frame(
    label = c("Kaplan-Meier", "Cure fraction"),
    color = c(color, cure_color),
    linetype = c("solid", "solid")
  )

  # Gráfico
  p <- ggplot(km_data, aes(x = time)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = ci_color,
      alpha = ci_alpha
    ) +
    geom_step(aes(y = surv), color = color, linewidth = line_size) +
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
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0.02)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
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
    labs(
      x = xlab,
      y = ylab,
      title = title,
      subtitle = if (show_cure)
        paste0(cure_label, format(round(cure_prob, 4), nsmall = 1)) else NULL
    ) +
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

# Waring Fragility Model - Random Activation
# G is the baseline distribution
# rho > 2. Se rho < 2, a variância da distribuição é infinita.
# a > 0
surv_wfm_random <- function(G) {
  function(t, rho, a, ...) {
    cure <- rho / (a + rho)
    cure + (1 - cure) * (1 - G(t, ...))
  }
}
desity_wfm_random <- function(g) {
  function(t, rho, a, ...) {
    cure <- rho / (a + rho)
    (1 - cure) * g(t, ...)
  }
}

# Criando automaticamente a função de sobrevivência do modelo com
# ativação aleatória, com G ~ Weibull
# Aqui você cria qualquer sobrevivencia, distribuição acumulada e densidade
# do seu modelo.
surv_wfm_random_weibull <- surv_wfm_random(G = pweibull) # Basta passar a G

cdf_wfm_random_weibull <- function(t, rho, a, shape, scale) {
  1 - surv_wfm_random_weibull(t, rho, a, shape = shape, scale = scale)
}

pdf_wfm_random_weibull <- desity_wfm_random(g = dweibull) # Basta passar a g

random_cure <- function(
  n = 1L,
  quantile_function = NULL,
  surv,
  censoring_cdf = NULL,
  args_censoring_cdf = NULL,
  args_model,
  prop_censoring = NULL
) {
  # Função auxiliar para encontrar o limite de cura
  find_cure_limit <- function(surv, args_model, tol = 1e-6, t_max = 1e5) {
    s <- function(t) do.call(surv, c(list(t), args_model))
    pi_hat <- s(t_max)

    t_upper <- 1
    while (abs(s(t_upper) - pi_hat) > tol && t_upper < t_max) {
      t_upper <- t_upper * 2
    }

    if (t_upper >= t_max) return(t_max)

    uniroot(
      function(t) abs(s(t) - pi_hat) - tol,
      lower = 1,
      upper = t_upper,
      tol = .Machine$double.eps^0.5
    )$root
  }

  # Obtendo a fração de cura
  surv_func <- purrr::partial(.f = surv, !!!args_model)
  cure_fraction <- surv_func(.Machine$double.xmax)

  # Verificação do parâmetro prop_censoring
  if (!is.null(prop_censoring)) {
    if (prop_censoring < cure_fraction) {
      warning(
        "prop_censoring (",
        prop_censoring,
        ") cannot be smaller than cure fraction (",
        cure_fraction,
        ")"
      )
    }
    if (prop_censoring > 1) {
      stop("prop_censoring cannot be greater than 1")
    }
  }

  # Gerando os tempos de falha reais
  if (is.function(quantile_function)) {
    quantile_func <- purrr::partial(.f = quantile_function, !!!args_model)
    u <- runif(n = n, min = 0, max = 1 - cure_fraction)
    t <- quantile_func(u)
  } else {
    find_root <- function(u, func) {
      vapply(
        u,
        function(u) {
          uniroot(
            function(t) func(t) - u,
            lower = 0,
            upper = find_cure_limit(surv, args_model),
            tol = .Machine$double.eps^0.5
          )$root
        },
        double(1L)
      )
    }
    u <- cure_fraction + (1 - cure_fraction) * runif(n = n)
    t <- find_root(u, surv_func)
  }

  # Determinando status de cura
  c <- rbinom(n = n, size = 1L, prob = 1 - cure_fraction)
  real_t <- ifelse(c, yes = t, no = .Machine$double.xmax)

  # Calculando o número necessário de censuras
  if (!is.null(prop_censoring)) {
    n_cens <- round(n * prop_censoring)
    n_cure <- sum(!c)
    n_suscept_cens <- max(0, n_cens - n_cure)

    # Ordenando os suscetíveis pelo tempo de falha
    suscept_idx <- which(c == 1L)
    ordered_idx <- suscept_idx[order(t[suscept_idx])]

    # Determinando quem será censurado
    if (n_suscept_cens > 0L) {
      cens_idx <- ordered_idx[1L:n_suscept_cens]
    } else {
      cens_idx <- integer(0L)
    }

    # Todos os curados são automaticamente censurados
    cure_idx <- which(c == 0L)

    # Juntando todos os índices de censura
    all_cens_idx <- c(cure_idx, cens_idx)

    # Gerando tempos de censura
    if (is.function(censoring_cdf)) {
      censoring_func <- purrr::partial(
        .f = censoring_cdf,
        !!!args_censoring_cdf
      )
      inv_censoring <- function(p) {
        uniroot(
          function(t) censoring_func(t) - p,
          lower = 0,
          upper = max(t),
          tol = .Machine$double.eps^0.5
        )$root
      }
      time_c <- sapply(runif(n), inv_censoring)
    } else {
      t_max <- max(real_t[c == 1L])
      time_c <- runif(n, min = 0, max = t_max)
    }

    # Aplicando censura determinística
    delta <- rep(1L, n)
    delta[all_cens_idx] <- 0L
    t_observed <- real_t
    t_observed[all_cens_idx] <- pmin(real_t[all_cens_idx], time_c[all_cens_idx])
  } else {
    # Comportamento original quando prop_censoring não é especificado
    if (is.function(censoring_cdf)) {
      censoring_func <- purrr::partial(
        .f = censoring_cdf,
        !!!args_censoring_cdf
      )
      inv_censoring <- function(p) {
        uniroot(
          function(t) censoring_func(t) - p,
          lower = 0,
          upper = max(t),
          tol = .Machine$double.eps^0.5
        )$root
      }
      time_c <- sapply(runif(n), inv_censoring)
    } else {
      t_max <- max(t[c == 1L])
      time_c <- runif(n, min = 0, max = t_max)
    }

    t_observed <- pmin(real_t, time_c)
    delta <- ifelse(t_observed < real_t, yes = 0L, no = 1L)
  }

  message("Cure fraction: ", round(cure_fraction, 4))

  return(tibble::tibble(t = t_observed, delta = delta))
}

# Função auxiliar para quantile_censoring (definida fora para evitar repetição)
quantile_censoring <- function(u, cdf, t_max) {
  u_adjusted <- u * cdf(t_max)
  vapply(
    X = u_adjusted,
    FUN = \(p) {
      uniroot(
        \(t) cdf(t) - p,
        lower = 0,
        upper = t_max,
        tol = .Machine$double.eps^0.5,
        extendInt = "upX"
      )$root
    },
    FUN.VALUE = double(1L)
  )
}

q_dagum <- function(u, theta, alpha, beta) {
  ((u * theta) / (beta * (theta - u)))^(1 / alpha)
}

surv_dagum <- function(t, theta, beta, alpha) {
  1 - theta * beta / (beta + theta * t^(-alpha))
}

theta <- 0.4
beta = 1.2
alpha <- 2.5

set.seed(0)
dados <-
  random_cure(
    n = 1000L,
    surv = surv_dagum,
    quantile_function = q_dagum,
    args_model = c(theta = theta, beta = beta, alpha = alpha),
    prop_censoring = 0.74 # Proporção de zeros
  )

dados |>
  dplyr::group_by(delta) |>
  dplyr::summarise(n = dplyr::n(), prop = n / nrow(dados))

plot_kaplan(dados, surv = surv_dagum, theta = theta, alpha = alpha, beta = beta)
