library(purrr)
library(rlang)
library(tibble)
library(survival)
library(ggplot2)

# Função usada para construir o Kaplan-Meier mais fácil no ggplot2
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
    geom_step() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    labs(
      x = "Tempo",
      y = "Probabilidade de Sobrevivência",
      title = "Curva de Kaplan-Meier"
    ) +
    stat_summary(
      geom = "line",
      fun = \(i) cura_real,
      aes(color = "Cura real"),
      linewidth = 1.5,
      linetype = 2
    ) +
    scale_y_continuous(limits = c(0, 1), n.breaks = 10) +
    labs(x = "t") +
    guides(color = guide_legend(title = ""))
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

    find_root <- function(u, func) {
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

    find_root <- Vectorize(FUN = find_root, vectorize.args = c("u"))

    if (is.function(quantile_function)) {
      quantile_function <- purrr::partial(.f = quantile_function, !!!args_model)
      u <- runif(n = n, min = 0, max = 1 - cure_fraction)
    } else if (is.function(surv)) {
      u <- cure_fraction + (1 - cure_fraction) * runif(n = n, min = 0, max = 1) # Utilizado para gerar os tempos reais
      quantile_function <- purrr::partial(.f = find_root, func = surv)
    } else {
      stop(
        "Invalid input: You must provide either a valid quantile function (quantile_function) or a cumulative distribution function (cdf)."
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

      quantile_censoring_single <- purrr::partial(
        .f = find_root,
        func = censoring_cdf
      )

      # Garante que o tempo de censura nunca ultrapasse t_max
      generate_censored_time <- function(u) {
        repeat {
          t <- quantile_censoring_single(u)
          if (t <= t_max) break
          u <- runif(1) # novo u aleatório
        }
        t
      }

      # Vetoriza corretamente
      quantile_censoring <- Vectorize(generate_censored_time)

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

dados <-
  random_cure(
    n = 1000L,
    surv = surv_wfm_random_weibull,
    args_model = c(rho = rho, a = a, shape = shape, scale = scale),
    t_max = 25
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

theta <- 0.6
beta = 1.2
alpha <- 2.5

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

dados <-
  random_cure(
    n = 2000L,
    surv = surv_dagum,
    quantile_function = q_dagum,
    censoring_cdf = pweibull,
    args_censoring_cdf = list(shape = 1.2, scale = 1.5),
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
