library(purrr)
library(rlang)
library(tibble)
library(survival)
library(ggplot2)
library(glue)

plot_kaplan <- function(dados, surv, ...) {
  cura_real <- surv(.Machine$double.xmax, ...)
  dados_surv <- Surv(time = dados$t, event = dados$delta)
  km_fit <- survfit(dados_surv ~ 1)

  km_data <- data.frame(
    time = km_fit$time,
    surv = km_fit$surv,
    upper = km_fit$upper,
    lower = km_fit$lower
  )

  ggplot(km_data, aes(x = time, y = surv)) +
    geom_step(linewidth = 1.2, color = "#2C3E50") +  # Azul escuro elegante
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#3498DB", alpha = 0.15) +  # Azul leve para IC
    geom_hline(aes(yintercept = cura_real, color = "Fração de cura"), 
               linetype = "dashed", linewidth = 1.2) +
    scale_color_manual(values = c("Fração de cura" = "#E74C3C")) +  # Vermelho suave
    guides(color = guide_legend(title = "")) +
    labs(
      x = "Tempo",
      y = "Probabilidade de Sobrevivência",
      title = "Curva de Kaplan-Meier com Fração de Cura"
    ) +
    scale_y_continuous(limits = c(0, 1), n.breaks = 10) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      legend.position = "top",
      legend.text = element_text(size = 12),
      panel.grid.major = element_line(color = "gray85"),
      panel.grid.minor = element_blank()
    )
}

# Modelo Waring com ativação aleatória
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

surv_wfm_random_weibull <- surv_wfm_random(G = pweibull)
pdf_wfm_random_weibull <- desity_wfm_random(g = dweibull)

# Gerador de dados com fração de cura e censura controlada
random_cure <- function(
  n = 1L,
  quantile_function = NULL,
  surv,
  censoring_cdf = NULL,
  args_censoring_cdf = NULL,
  args_model,
  prop_cens = NULL
) {
  find_cure_limit <- function(surv, args_model, tol = 1e-6, tol_time = 1e-6, t_init = 10, t_max = 1e5) {
    s <- function(t) do.call(surv, c(list(t), args_model))
    pi_hat <- s(t_max)
    t_lower <- t_init
    t_upper <- t_lower
    repeat {
      diff <- abs(s(t_upper) - pi_hat)
      if (diff < tol || t_upper >= t_max) break
      t_lower <- t_upper
      t_upper <- min(t_upper * 2, t_max)
    }
    if (t_upper >= t_max && abs(s(t_upper) - pi_hat) > tol) return(NA_real_)
    while ((t_upper - t_lower) > tol_time) {
      t_mid <- 0.5 * (t_lower + t_upper)
      s_mid <- s(t_mid)
      if (abs(s_mid - pi_hat) < tol) {
        t_upper <- t_mid
      } else {
        t_lower <- t_mid
      }
    }
    return(t_upper)
  }

  upper <- find_cure_limit(surv = surv, args_model)

  find_root <- function(u, func) {
    uniroot(
      \(t) func(t) - u,
      lower = 0,
      upper = upper,
      tol = .Machine$double.eps^5
    )$root
  }

  surv <- purrr::partial(.f = surv, !!!args_model)
  cure_fraction <- surv(.Machine$double.xmax)

  if(prop_cens < cure_fraction){
    warning(glue::glue("O valor de prop_cens deve ser maior que a fração de cura de {cure_fraction}"))
  }

  find_root <- Vectorize(FUN = find_root, vectorize.args = c("u"))

  if (is.function(quantile_function)) {
    quantile_function <- purrr::partial(.f = quantile_function, !!!args_model)
  } else if (is.function(surv)) {
    quantile_function <- purrr::partial(.f = find_root, func = surv)
  } else {
    stop("Invalid input: You must provide either a valid quantile function or a cumulative distribution function.")
  }

  c <- rbinom(n = n, size = 1L, prob = 1 - cure_fraction)

  u <- cure_fraction + (1 - cure_fraction) * runif(n = n)
  t <- quantile_function(u)

  idx_susc <- which(c == 1)
  idx_cure <- which(c == 0)
  n_susc <- length(idx_susc)
  n_cure <- length(idx_cure)

  t_susc <- t[idx_susc]
  t_max <- max(t_susc)

  t_observed <- numeric(n)
  delta <- integer(n)

  # ✅ Curados censurados no tempo máximo
  t_observed[idx_cure] <- t_max
  delta[idx_cure] <- 0L

  t_observed[idx_susc] <- t_susc
  delta[idx_susc] <- 1L

  if (!is.null(prop_cens)) {
    n_cens_total <- floor(n * prop_cens)
    n_cens_susc <- max(0, n_cens_total - n_cure)

    if (n_cens_susc > n_susc) {
      warning("Proporção de censura excede o possível dado o número de curados. Ajustando.")
      n_cens_susc <- n_susc
    }

    if (n_cens_susc > 0) {
      idx_cens_susc <- sample(idx_susc, size = n_cens_susc)
      t_observed[idx_cens_susc] <- runif(n_cens_susc, min = 0, max = t[idx_cens_susc])
      delta[idx_cens_susc] <- 0L
    }
  }

  print(cure_fraction)
  return(tibble::tibble(t = t_observed, delta = delta))
}

# Função de sobrevivência exemplo: Dagum
surv_dagum <- function(t, theta, beta, alpha) {
  1 - theta * beta / (beta + theta * t^(-alpha))
}

# Parâmetros
theta <- 0.6
beta <- 1.2
alpha <- 2.5

set.seed(0)
# Geração dos dados
dados <- random_cure(
  n = 1000L,
  surv = surv_dagum,
  args_model = c(theta = theta, beta = beta, alpha = alpha),
  prop_cens = 0.4 # Proporção de censura
)

# Proporção de censura observada
dplyr::count(dados, delta) |> dplyr::mutate(prop = n / sum(n))

# Plotando Kaplan-Meier com curva de cura real
plot_kaplan(
  dados,
  surv = surv_dagum,
  theta = theta,
  alpha = alpha,
  beta = beta
)
