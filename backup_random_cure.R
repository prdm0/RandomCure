library(purrr)
library(rlang)
library(tibble)
library(dplyr)
library(survival)

source("plot_kaplan_meier.R")

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


# ------------------------------------------------------------------------------
# random_cure: Geração com fração de cura e proporção total de zeros,
# usando Uniform(0, T_true) quando censoring_cdf = NULL,
# e default prop_zeros = cure_fraction se NULL
# ------------------------------------------------------------------------------
random_cure <- function(
  n = 1L,
  quantile_function = NULL, # função quantil u→t (ou NULL para inversão via surv)
  surv, # função S(t)
  censoring_cdf = NULL, # função F_C(c) (ou NULL → Uniform)
  args_censoring_cdf = NULL, # parâmetros da CDF de censura (lista)
  args_model,
  prop_zeros = NULL # proporção total de zeros (cura + censura)
) {
  # 1) fração de cura
  surv_fn <- purrr::partial(surv, !!!args_model)
  cure_fraction <- surv_fn(.Machine$double.xmax)

  # 2) default de prop_zeros
  if (is.null(prop_zeros)) {
    prop_zeros <- cure_fraction
  }
  if (prop_zeros < cure_fraction || prop_zeros > 1) {
    stop(sprintf(
      "valor inválido para prop_zeros = %.3f (%.3f <= prop_zeros <= 1), pois cure_fraction = %.3f!",
      prop_zeros,
      cure_fraction,
      cure_fraction
    ))
  }

  # 3) gera tempos de falha T_true
  if (!is.null(quantile_function)) {
    qf <- partial(quantile_function, !!!args_model)
    U1 <- runif(n, 0, 1 - cure_fraction)
    T_true <- qf(U1)
  } else {
    U1 <- cure_fraction + (1 - cure_fraction) * runif(n)
    T_true <- vapply(
      U1,
      function(u) {
        uniroot(
          function(t) surv_fn(t) - u,
          lower = 0,
          upper = 1e6,
          extendInt = "downX"
        )$root
      },
      double(1)
    )
  }

  # 4) curados vs suscetíveis
  is_susc <- rbinom(n, 1L, 1 - cure_fraction) == 1L
  T_true[!is_susc] <- Inf

  # 5) define fração de censura entre suscetíveis
  p_c <- (prop_zeros - cure_fraction) / (1 - cure_fraction)

  # 6) gera tempos de censura apenas entre suscetíveis
  T_cens <- rep(Inf, n)
  if (p_c > 0) {
    idx_s <- which(is_susc)
    m <- length(idx_s)
    k <- round(m * p_c)
    cens_idx <- sample(idx_s, size = k)

    if (!is.null(censoring_cdf) && !is.null(args_censoring_cdf)) {
      # inverte pela raiz na CDF
      F_c <- purrr::partial(censoring_cdf, !!!args_censoring_cdf)
      U2 <- runif(k)
      Fc_T <- F_c(T_true[cens_idx])
      p0 <- U2 * Fc_T # garante C < T_true
      T_cens[cens_idx] <- vapply(
        seq_along(p0),
        function(ii) {
          uniroot(
            function(t) F_c(t) - p0[ii],
            lower = 0,
            upper = T_true[cens_idx[ii]],
            tol = .Machine$double.eps^0.5
          )$root
        },
        double(1)
      )
    } else {
      # caso uniforme: C ~ Uniform(0, T_true)
      U2 <- runif(k)
      T_cens[cens_idx] <- U2 * T_true[cens_idx]
    }
  }

  # 7) observado e variável indicadora delta
  T_obs <- pmin(T_true, T_cens)
  delta <- as.integer((T_true <= T_cens) & is_susc)

  tibble(t = T_obs, delta = delta)
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
    prop_zeros = NULL
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
    prop_zeros = 0.5
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
    prop_zeros = 0.7
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
    n = 5000L,
    surv = surv_dagum,
    quantile_function = q_dagum,
    censoring_cdf = pweibull,
    args_censoring_cdf = list(shape = 2.5, scale = 2.3),
    args_model = c(theta = theta, beta = beta, alpha = alpha),
    prop_zeros = 0.6
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
