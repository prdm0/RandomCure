library(rlang)
library(tibble)
library(dplyr)
library(survival)
library(purrr)

rm(list = ls())

# Exemplo 4 (Dagum especificando quantílica e outra distribuição para a censura)
q_dagum <- function(u, theta, alpha, beta) {
  ((u * theta) / (beta * (theta - u)))^(1 / alpha)
}

surv_dagum <- function(t, theta, beta, alpha) {
  1 - theta * beta / (beta + theta * t^(-alpha))
}

theta <- 0.5
beta <- 1.2
alpha <- 2.5

n <- 5000L
surv <- surv_dagum
quantile_function <- q_dagum
censoring_cdf <- pweibull
args_censoring_cdf <- list(shape = 2.5, scale = 2.3)
args_model <- c(theta = theta, beta = beta, alpha = alpha)
prop_zeros <- 0.6

source("plot_kaplan_meier.R")

# Waring Fragility Model - Random Activation
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
cdf_wfm_random_weibull <- function(t, rho, a, shape, scale) {
  1 - surv_wfm_random_weibull(t, rho, a, shape = shape, scale = scale)
}
pdf_wfm_random_weibull <- desity_wfm_random(g = dweibull)

#' Geração de dados censurados com fração de cura
#'
#' Esta função gera um conjunto de observações (t, δ) sob um modelo de fração
#' de cura, permitindo ao usuário:
#' * controlar a proporção total de zeros (`prop_zeros`, soma de curados e censurados),
#' * usar uma CDF de censura arbitrária (ou regra uniforme interna),
#' * fornecer ou não a função quantílica de T.
#'
#' @param n Integer. Número de observações a gerar.
#' @param quantile_function Função quantílica de T em `u -> QT(u, ...)`. Se `NULL`,
#'   o algoritmo resolve `S(t) = π + u` numericamente.
#' @param surv Função de sobrevivência de T em `t -> ST(t, ...)`.
#' @param censoring_cdf Função CDF de censura `c -> FC(c, ...)`. Se `NULL`, usa
#'   censura uniforme em [0, T_i].
#' @param args_censoring_cdf Lista de argumentos adicionais para `censoring_cdf`.
#' @param args_model Vetor nomeado com parâmetros a passar a `surv` e
#'   `quantile_function` (se fornecida).
#' @param prop_zeros Numeric em [π, 1], proporção total de zeros desejada
#'   (curados + censurados). Se `NULL`, adota-se `prop_zeros = π`.
#'
#' @return Um `tibble` com colunas:
#'   * `t`: tempo observado (mínimo entre T e C, ou `Inf` para curados puros),
#'   * `delta`: indicador de evento (1 se falha, 0 se censurado ou curado).
#'
#' @details
#' 1. Calcula a fração de cura teórica π = ST(∞) via busca binária em `find_cure_limit`.
#' 2. Define `p_c = (prop_zeros - π)/(1 - π)`, fração de censura entre suscetíveis.
#' 3. Gera tempos verdadeiros Tᵢ por quantilização ou raiz numérica.
#' 4. Divide em curados (Tᵢ = ∞) e suscetíveis (Mi ~ Bernoulli(1 - π)).
#' 5. Escolhe `k = round(|I| * p_c)` suscetíveis para censurar.
#' 6. Gera Ci:
#'    - Se `censoring_cdf` é função: Vi ~ U(0, FC(Tᵢ)), Ci = F⁻¹_C(Vi);
#'    - Senão: Ci ~ U(0, Tᵢ).
#' 7. Observação final: tᵢ = min(Tᵢ, Ci), δᵢ = 1{Tᵢ ≤ Ci}.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' dados <- random_cure(
#'   n = 1000,
#'   surv = surv_dagum,
#'   quantile_function = q_dagum,
#'   censoring_cdf = pweibull,
#'   args_censoring_cdf = list(shape = 1.5, scale = 2),
#'   args_model = c(theta = 0.5, beta = 1.2, alpha = 2.5),
#'   prop_zeros = 0.6
#' )
#' }
#' @export
random_cure <- function(
  n = 1L,
  quantile_function = NULL,
  surv,
  censoring_cdf = NULL,
  args_censoring_cdf = NULL,
  args_model,
  prop_zeros = NULL
) {
  #------------------------------------------------------------
  # 1) Encontrar um "upper" grande para o domínio de S(t)
  #------------------------------------------------------------
  find_cure_limit <- function(
    surv,
    tol = 1e-6,
    tol_time = 1e-6,
    t_init = 10,
    t_max = 1e5
  ) {
    # Aproximação da fração de cura
    pi_hat <- surv(t_max)

    # Busca exponencial para achar intervalo
    t_lower <- t_init
    t_upper <- t_lower
    repeat {
      diff <- abs(surv(t_upper) - pi_hat)
      if (diff < tol || t_upper >= t_max) {
        break
      }
      t_lower <- t_upper
      t_upper <- min(t_upper * 2, t_max)
    }

    if (t_upper >= t_max && abs(surv(t_upper) - pi_hat) > tol) {
      return(NA_real_)
    }

    # Busca binária refinada
    while ((t_upper - t_lower) > tol_time) {
      t_mid <- 0.5 * (t_lower + t_upper)
      s_mid <- surv(t_mid)
      if (abs(s_mid - pi_hat) < tol) {
        t_upper <- t_mid
      } else {
        t_lower <- t_mid
      }
    }
    # Retorna o menor valor de t tal que a sobrevivência converge
    # para a cura do modelo.
    t_upper
  }

  # "Congela" os parâmetros do modelo em surv
  surv <- purrr::partial(surv, !!!args_model)

  upper <- find_cure_limit(surv = surv)
  if (is.na(upper)) {
    upper <- 1e5
  } # fallback defensivo
  cure_fraction <- surv(upper * 10) # π ≈ S(∞)

  #------------------------------------------------------------
  # 2) Validar prop_zeros
  #    (por padrão, só curados; sem censura extra entre suscetíveis)
  #------------------------------------------------------------
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

  #------------------------------------------------------------
  # 3) Utilidades de inversão (root-finding) com upper
  #------------------------------------------------------------
  find_root_no_vec <- function(u, func) {
    uniroot(
      \(t) func(t) - u,
      lower = 0,
      upper = upper,
      tol = .Machine$double.eps^0.5
    )$root
  }
  find_root <- function(u, func) {
    vapply(X = u, FUN = \(ui) find_root_no_vec(ui, func), FUN.VALUE = double(1))
  }

  #------------------------------------------------------------
  # 4) Preparar CDF de censura (se fornecida)
  #------------------------------------------------------------
  if (is.function(censoring_cdf)) {
    censoring_cdf <- purrr::partial(censoring_cdf, !!!args_censoring_cdf)
  }

  #------------------------------------------------------------
  # 5) Quantil de T (tempo verdadeiro para suscetíveis)
  #------------------------------------------------------------
  if (is.function(quantile_function)) {
    quantile_function_real_time <- function(n) {
      q <- purrr::partial(quantile_function, !!!args_model)
      # amostra apenas no segmento dos suscetíveis (massa 1 - π)
      q(runif(n, 0, 1 - cure_fraction))
    }
  } else {
    # inversão via sobrevivência: S(T) = u, com u ~ (π, 1)
    quantile_function_real_time <- function(n) {
      u <- cure_fraction + (1 - cure_fraction) * runif(n)
      purrr::partial(find_root, func = surv)(u)
    }
  }

  #------------------------------------------------------------
  # 6) Quantil de C (para suscetíveis selecionados para censura)
  #    - Se houver CDF: C|T=t por inversão em [0, F_C(t)]
  #    - Senão: U(0, t)
  #------------------------------------------------------------
  if (is.function(censoring_cdf)) {
    quantile_function_censoring <- function(t_i) {
      Fc <- pmin(1, censoring_cdf(t_i))
      u <- runif(length(t_i), 0, Fc)
      purrr::partial(find_root, func = censoring_cdf)(u)
    }
  } else {
    quantile_function_censoring <- function(t_i) {
      runif(length(t_i), 0, t_i)
    }
  }

  #------------------------------------------------------------
  # 7) Fração de censura entre suscetíveis
  #    (proporção extra de zeros só entre suscetíveis)
  #------------------------------------------------------------
  p_c <- (prop_zeros - cure_fraction) / (1 - cure_fraction)

  #------------------------------------------------------------
  # 8) Geração dos tempos verdadeiros e suscetibilidade
  #------------------------------------------------------------
  t_true <- quantile_function_real_time(n)
  m <- rbinom(n = n, size = 1L, prob = 1 - cure_fraction) # 1=suscetível, 0=curado
  t_true <- ifelse(m == 0L, Inf, t_true) # curado tem T=Inf (conceitual)

  #------------------------------------------------------------
  # 9) Seleção de suscetíveis a censurar (excesso de zeros)
  #------------------------------------------------------------
  n_susc <- sum(m == 1L)
  k <- round(n_susc * p_c) # censura só entre suscetíveis
  id <- if (k > 0 && n_susc > 0) {
    sample(which(m == 1L), size = k)
  } else {
    integer(0)
  }

  #------------------------------------------------------------
  # 10) Vetor de censura t_c
  #     - suscetíveis sorteados: censura condicional (ou U(0,T))
  #     - curados: censura administrativa tardia => tempos finitos, nunca Inf
  #------------------------------------------------------------
  t_c <- rep(Inf, n)

  # 10a) Censura dos suscetíveis sorteados
  if (length(id)) {
    t_c[id] <- quantile_function_censoring(t_true[id])
  }

  # 10b) Censura administrativa tardia para curados (não altera contagem de zeros)
  idx_cured <- which(m == 0L)
  if (length(idx_cured)) {
    # horizonte administrativo: maior T verdadeiro entre suscetíveis (ou 'upper' se não houver)
    t_max_susc <- if (any(is.finite(t_true[m == 1L]))) {
      max(t_true[m == 1L])
    } else {
      upper
    }
    t_admin <- t_max_susc
    t_c[idx_cured] <- t_admin # todos os curados censurados no fim do seguimento
  }

  #------------------------------------------------------------
  # 11) Observação final
  #------------------------------------------------------------
  t <- pmin(t_true, t_c)
  delta <- ifelse(t_true < t_c, 1L, 0L)

  tibble::tibble(t = t, delta = delta)
}

# Teste rápido
rho <- 2.5
a <- 1.5
shape <- 1.5
scale <- 2.5
set.seed(0)
dados <- random_cure(
  n = 10000L,
  surv = surv_wfm_random_weibull,
  args_model = c(rho = rho, a = a, shape = shape, scale = scale),
  prop_zeros = 0.8
)

plot_kaplan(
  dados,
  surv = surv_wfm_random_weibull,
  rho = rho,
  a = a,
  shape = shape,
  scale = scale
)

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
    n = 4000L,
    surv = surv_dagum,
    args_model = c(theta = theta, beta = beta, alpha = alpha),
    prop_zeros = 0.7
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
    prop_zeros = 0.67
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
    prop_zeros = NULL
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
