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
  cure_color = "#D62728", # Vivid red
  base_size = 11,
  ci_alpha = 0.2,
  line_size = 0.7,
  show_cure = TRUE,
  cure_label = "Cure fraction: ",
  ...
) {
  require(ggplot2)
  require(survival)

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
        paste0(cure_label, format(round(cure_prob, 4), nsmall = 1)) else NULL
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

# Funcao que usaremos para gerar os dados
random_cure <-
  function(
    n = 1L,
    quantile_function = NULL,
    surv,
    censoring_cdf = NULL, # Se não for passada, é a distribuição uniforme (censura não-informativa)
    args_censoring_cdf = NULL, # Lista de valores dos argumentos da distribução da censura.
    args_model,
    t_max = NULL
  ) {
    find_cure_limit <- function(
      surv,
      args_model,
      tol = 1e-6,
      tol_time = 1e-6,
      t_init = 10,
      t_max = 1e5
    ) {
      s <- function(t) do.call(surv, c(list(t), args_model))

      # Aproximação da fração de cura
      pi_hat <- s(t_max)

      # Busca exponencial para achar intervalo
      t_lower <- t_init
      t_upper <- t_lower

      repeat {
        diff <- abs(s(t_upper) - pi_hat)
        if (diff < tol || t_upper >= t_max) break
        t_lower <- t_upper
        t_upper <- min(t_upper * 2, t_max)
      }

      if (t_upper >= t_max && abs(s(t_upper) - pi_hat) > tol) return(NA_real_)

      # Busca binária refinada
      while ((t_upper - t_lower) > tol_time) {
        t_mid <- 0.5 * (t_lower + t_upper)
        s_mid <- s(t_mid)
        if (abs(s_mid - pi_hat) < tol) {
          t_upper <- t_mid
        } else {
          t_lower <- t_mid
        }
      }
      # Retorna o menor valor de t tal que a sobrevivência converge
      # para cura do modelo.
      return(t_upper)
    }

    upper <- find_cure_limit(surv = surv, args_model)

    find_root_scalar <- function(u, func) {
      uniroot(
        \(t) func(t) - u,
        lower = 0,
        upper = upper,
        tol = .Machine$double.eps^0.5
      )$root
    }

    # Obtendo a fração de cura da função de sobrevivência
    surv <- purrr::partial(.f = surv, !!!args_model)
    cure_fraction <- surv(.Machine$double.xmax)

    find_root <- function(u, func) {
      vapply(
        X = u,
        FUN = function(u) find_root_scalar(u, func),
        FUN.VALUE = numeric(1)
      )
    }

    # find_root <- Vectorize(FUN = find_root_scalar, vectorize.args = c("u"))

    if (is.function(quantile_function)) {
      quantile_function <- purrr::partial(.f = quantile_function, !!!args_model)
      u <- runif(n = n, min = 0, max = 1 - cure_fraction)
    } else if (is.function(surv)) {
      u <- cure_fraction + (1 - cure_fraction) * runif(n = n, min = 0, max = 1) # Utilizado para gerar os tempos reais
      quantile_function <- purrr::partial(.f = find_root, func = surv)
    } else {
      stop(
        "Invalid input: You must provide either a valid quantile function (quantile_function) and/or a survival function (surv)."
      )
    }

    # Dividindo a população em idivíduos suscetíveis e curados
    c <- rbinom(n = n, size = 1L, prob = 1 - cure_fraction) # Divisão da população

    # Calculando o tempo de sobrevivência real
    t <- quantile_function(u) # Tempo real

    # Atribuindo tempo máximo (maior real da arquitetura) para os indivíduos
    # curados.
    real_t <- ifelse(c, yes = t, no = .Machine$double.xmax)

    if (is.null(t_max)) {
      t_max <- max(t * c)
    }

    if (is.function(censoring_cdf)) {
      censoring_cdf <- purrr::partial(.f = censoring_cdf, !!!args_censoring_cdf)

      # Função vetorizada que já considera t_max
      quantile_censoring <- function(u) {
        # Ajusta os valores de u para não ultrapassar censoring_cdf(t_max)
        u_adjusted <- u * censoring_cdf(t_max)

        # Encontra todos os roots de uma vez (mais eficiente)
        vapply(
          X = u_adjusted,
          FUN = \(p) {
            uniroot(
              \(t) censoring_cdf(t) - p,
              lower = 0,
              upper = t_max * 10, # Limite ampliado para evitar falhas
              tol = .Machine$double.eps^0.5,
              extendInt = "upX"
            )$root
          },
          FUN.VALUE = double(1L)
        )
      }

      time_c <- quantile_censoring(runif(n = n))
    } else {
      time_c <- runif(n = n, min = 0, max = t_max)
    }

    t_observed <- pmin(real_t, time_c)
    delta <- ifelse(t_observed < real_t, yes = 0L, no = 1L) # 1 para falha e 0 para censura
    print(cure_fraction)
    return(tibble::tibble(t = t_observed, delta = delta))
  }

# Exemplo 1 (PVF sem especificar a quantílica) ---------------------------

rho <- 2.5 # rho > 2
a <- 1.5
shape <- 1.5
scale <- 2.5

set.seed(0)
dados <-
  random_cure(
    n = 1000L,
    surv = surv_wfm_random_weibull,
    args_model = c(rho = rho, a = a, shape = shape, scale = scale),
    t_max = NULL
  )

plot_kaplan(
  dados,
  surv = surv_wfm_random_weibull,
  rho = rho,
  a = a,
  shape = shape,
  scale = scale
)

# Exemplo 2 (Dagum sem especificar a quantílica) -------------------------

surv_dagum <- function(t, theta, beta, alpha) {
  1 - theta * beta / (beta + theta * t^(-alpha))
}

theta <- 0.61
beta = 1.2
alpha <- 2.5

set.seed(0)
dados <-
  random_cure(
    n = 1000L,
    surv = surv_dagum,
    args_model = c(theta = theta, beta = beta, alpha = alpha),
    t_max = NULL
  )

dados |>
  dplyr::group_by(delta) |>
  dplyr::summarise(n = dplyr::n(), prop = n / nrow(dados))

plot_kaplan(
  dados,
  surv = surv_dagum,
  theta = theta,
  alpha = alpha,
  beta = beta
)

# Exemplo 3 (Dagum especificando a quantílica) ---------------------------

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
    t_max = NULL
  )

plot_kaplan(
  dados,
  surv = surv_dagum,
  theta = theta,
  alpha = alpha,
  beta = beta
)
dados |>
  dplyr::group_by(delta) |>
  dplyr::summarise(n = dplyr::n(), prop = n / nrow(dados))


# Exemplo 4 (Dagum especificando quantilica e outra distribuição para a censura ----

q_dagum <- function(u, theta, alpha, beta) {
  ((u * theta) / (beta * (theta - u)))^(1 / alpha)
}

surv_dagum <- function(t, theta, beta, alpha) {
  1 - theta * beta / (beta + theta * t^(-alpha))
}

theta <- 0.5
beta = 1.2
alpha <- 2.5

set.seed(0)
dados <-
  random_cure(
    n = 1000L,
    surv = surv_dagum,
    quantile_function = q_dagum,
    censoring_cdf = pweibull,
    args_censoring_cdf = list(shape = 1.5, scale = 1.3),
    args_model = c(theta = theta, beta = beta, alpha = alpha)
  )

plot_kaplan(
  dados,
  surv = surv_dagum,
  theta = theta,
  alpha = alpha,
  beta = beta
)

dados |>
  dplyr::group_by(delta) |>
  dplyr::summarise(n = dplyr::n(), prop = n / nrow(dados))
