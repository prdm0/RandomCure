library(purrr)
library(rlang)
library(tibble)

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