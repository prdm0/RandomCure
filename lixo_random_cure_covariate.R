library(purrr)
library(rlang)
library(tibble)
library(ggplot2)
library(survival)

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

# Funcao que usaremos para gerar os dados
random_cure_covariate <-
  function(n = 1L,
           quantile_function,
           cdf = NULL,
           t_max = NULL,
           uniform_censoring = TRUE,
           lower = 0,
           upper = 50,
           covariate_names = NULL,
           args_model) {
    # Geração da covariável x
    x <- rexp(n, rate = 1.5) #rbinom(n = n, size = 1L, prob = 0.5)

    # Ajuste dos parâmetros em args_model com base nas covariáveis
    if (!is.null(covariate_names)) {
      for (name in covariate_names) {
        const_1 <- args_model[[paste0(name, "_0")]]
        const_2 <- args_model[[paste0(name, "_1")]]
        args_model[[name]] <- mean(exp(const_1 + const_2 * x))
      }
    }

    find_root <- function(u, ...) {
      uniroot(\(x) cdf(x, ...) - u, lower = lower, upper = upper, tol = .Machine$double.eps^0.5)$root
    }
    find_root <- Vectorize(FUN = find_root, vectorize.args = c("u"))

    cure_fraction <- args_model$rho / (args_model$rho + args_model$a)

    c <- rbinom(n = n, size = 1L, prob = 1 - cure_fraction)
    u <- runif(n = n, min = 0, max = 1 - cure_fraction)

    args_model_filtered <- args_model[!grepl("_(0|1)$", names(args_model))]

    if (is.function(quantile_function)) {
      quantile_function <- purrr::partial(.f = quantile_function, !!!args_model_filtered)
    } else if (is.function(cdf)) {
      quantile_function <- purrr::partial(.f = find_root, !!!args_model_filtered)
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
    return(tibble::tibble(t = t_observed, delta = delta, x = x))
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

# Valores verdadeiros
rho_0 <- 1.1
rho_1 <- 1.1
a_0 <- 1.5
a_1 <- 1.5
shape_0 <- 0.1
shape_1 <- 0.1
scale_0 <- 0.5
scale_1 <- 0.5

set.seed(0)
dados <- random_cure_covariate(
  n = 5000L,
  quantile_function = NULL,
  cdf = cdf_wfm_random_weibull,
  t_max = NULL,
  lower = 0,
  upper = 30,
  uniform_censoring = TRUE,
  covariate_names = c("rho", "a", "shape", "scale"),
  args_model = list(
    a_0 = a_0,
    a_1 = a_1,
    rho_0 = rho_0,
    rho_1 = rho_1,
    shape_0 = shape_0,
    shape_1 = shape_1,
    scale_0 = scale_0,
    scale_1 = scale_1
  )
)

x <- dados$x
a <- mean(exp(a_0 + a_1 * x))
rho <- mean(exp(rho_0 + rho_1 * x))
shape <- mean(exp(shape_0 + shape_1 * x))
scale <- mean(exp(scale_0 + scale_1 * x))
cura <- rho / (a + rho)

p1 <- plot_kaplan(dados, cura)

# Estimação ---------------------------------------------------------------
# Função de log-verossimilhança
log_likelihood <- function(t, delta, rho, a, shape, scale) {
  term_1 <- delta * log(pdf_wfm_random_weibull(t = t, rho = rho, a = a, shape = shape, scale = scale))
  term_2 <- (1 - delta) * log(surv_wfm_random_weibull(t = t, rho = rho, a = a, shape = shape, scale = scale))
  -sum(term_1 + term_2)
}

try_log_likelihood <- function(...) {
  tryCatch(
    log_likelihood(...),
    error = function(e) Inf
  )
}

emv <-
  optim(
    par = c(3, 1.5, 1.5, 1.5),
    fn = \(par) try_log_likelihood(t = dados$t, delta = dados$delta, rho = par[1], a = par[2], shape = par[3], scale = par[4]),
    method = "BFGS"
  )

# Adiciona ao gráfico uma curva de sobrevivência estimada
p1 +
  stat_function(
    fun = \(t) surv_wfm_random_weibull(t, rho = emv$par[1], a = emv$par[2], shape = emv$par[3], scale = emv$par[4]),
    color = "blue",
    linewidth = 1,
    aes(color = "Sobrevivência estimada")
  )