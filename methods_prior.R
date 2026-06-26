## ============================================================
## Prior-induced median and PESS for baseline toxicity curve
## Models:
##   1) Regular/AIDE power CRM:      p_j = q_j^exp(theta)
##   2) alpha-CRM baseline:          p_j = S(d_j)^exp(theta) = q_j^exp(theta)
##   3) logistic cumulative/IPCRM:   logit(p_j) = beta0 + beta1*x_j + beta2*0
##   4) CFO/PRIDE:                   logit(p_ij) = beta_j + W_i
##
## Baseline only:
##   alpha-CRM: D_ij = current dose d_j
##   IPCRM/cumu CRM: cumu.d = 0
## ============================================================

set.seed(20260617)

## -----------------------------
## User-specified inputs
## -----------------------------
K <- 5

## skeleton for power CRM and alpha-CRM
# q_skeleton <- c()
q_skeleton <- c(0.15, 0.20, 0.30, 0.35, 0.45)
## alpha-CRM actual dose amounts, mg
# dose_alpha_mg <- c(10, 20, 30, 40, 50, 60, 70, 80)
dose_alpha_mg <- c(15, 20, 30, 35, 45)
## IPCRM / cumulative CRM current dose scores
# dose_ipcrm <- c(0.1, 0.3, 0.5, 0.7, 0.9)
dose_ipcrm <- c(15, 20, 30, 35, 45)
# dose_ipcrm <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5)
dose_ipcrm <- dose_ipcrm / (2 * sd(dose_ipcrm))

## prior draws
n_prior <- 200000

## -----------------------------
## Priors
## -----------------------------

## Power CRM / alpha-CRM prior:
## theta ~ N(0, 2)
theta_mean <- 0
theta_sd   <- sqrt(2)

## IPCRM / logistic cumulative CRM priors from your code:
## beta0 ~ t(fixed.intercept, precision = 4, df = 1)
## beta1 ~ Gamma(5.83, 1.21)
## beta2 ~ Exp(1), but beta2 drops out at baseline because cumu.d = 0
fixed_intercept <- -2.8
beta0_prec <- 2
beta0_df   <- 1

beta1_shape <- 2.5
beta1_rate  <- 1.6

## CFO / PRIDE priors:
## beta_j ~ N(mu_j, sigma2_beta), mu_j = logit(c_j)
## inv_sigma2_w ~ Gamma(eta, eta), W_i ~ N(0, sigma2_w)
cfo_skeleton <- c(0.005, 0.01, 0.05, 0.1, 0.3)
cfo_sigma2_beta <- 30
cfo_eta <- 1
## -----------------------------
## Helper functions
## -----------------------------

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

moment_pess <- function(p) {
  m <- mean(p)
  v <- stats::var(p)
  pess <- m * (1 - m) / v - 1
  pess
}

prior_summary <- function(pmat, model_name, dose_label) {
  ## pmat: n_prior by K matrix
  out <- data.frame(
    model = model_name,
    dose = dose_label,
    # mean = apply(pmat, 2, mean),
    median = apply(pmat, 2, median),
    # q025 = apply(pmat, 2, quantile, probs = 0.025),
    # q975 = apply(pmat, 2, quantile, probs = 0.975),
    # var = apply(pmat, 2, var),
    PESS = apply(pmat, 2, moment_pess),
    row.names = NULL
  )
  out
}

interp_S_vec <- function(D, d_grid, s_grid) {
  K <- length(d_grid)
  out <- numeric(length(D))
  
  idx_lo <- D <= d_grid[1]
  out[idx_lo] <- s_grid[1]
  
  idx_hi <- D >= d_grid[K]
  if (any(idx_hi)) {
    t <- (D[idx_hi] - d_grid[K - 1]) / (d_grid[K] - d_grid[K - 1])
    out[idx_hi] <- pmin(
      1,
      s_grid[K - 1] + t * (s_grid[K] - s_grid[K - 1])
    )
  }
  
  idx_mid <- !(idx_lo | idx_hi)
  if (any(idx_mid)) {
    Dm <- D[idx_mid]
    k <- findInterval(Dm, d_grid)
    k <- pmax(1, pmin(k, K - 1))
    t <- (Dm - d_grid[k]) / (d_grid[k + 1] - d_grid[k])
    out[idx_mid] <- s_grid[k] + t * (s_grid[k + 1] - s_grid[k])
  }
  
  out
}

## -----------------------------
## 1) Regular / AIDE power CRM
## -----------------------------
theta <- rnorm(n_prior, mean = theta_mean, sd = theta_sd)

p_power <- sapply(q_skeleton, function(qj) {
  qj ^ exp(theta)
})

summary_power <- prior_summary(
  pmat = p_power,
  model_name = "power_crm",
  dose_label = paste0("dose", seq_len(K))
)

## -----------------------------
## 2) alpha-CRM baseline
## -----------------------------
## Baseline alpha-CRM assumes D_ij = current dose d_j.
## Therefore S(D_ij) = q_j at the grid points.
S_baseline <- interp_S_vec(
  D = dose_alpha_mg,
  d_grid = dose_alpha_mg,
  s_grid = q_skeleton
)

p_alpha <- sapply(S_baseline, function(Sj) {
  Sj ^ exp(theta)
})

summary_alpha <- prior_summary(
  pmat = p_alpha,
  model_name = "alpha_crm_baseline",
  dose_label = paste0(dose_alpha_mg, " mg")
)

## -----------------------------
## 3) Logistic cumulative CRM / IPCRM baseline
## -----------------------------
## Baseline means cumu.d = 0, so beta2 does not enter.
beta0 <- fixed_intercept + rt(n_prior, df = beta0_df) / sqrt(beta0_prec)
beta1 <- rgamma(n_prior, shape = beta1_shape, rate = beta1_rate)

p_ipcrm <- sapply(dose_ipcrm, function(xj) {
  inv_logit(beta0 + beta1 * xj)
})

summary_ipcrm <- prior_summary(
  pmat = p_ipcrm,
  model_name = "ipcrm_logistic_baseline",
  dose_label = paste0("x=", dose_ipcrm)
)

## -----------------------------
## 4) CFO / PRIDE induced prior
## -----------------------------
## Match PRIDE.R's prior-only approximation:
## p_j = E_W[logit^{-1}(beta_j + W)] is approximated by
## logit^{-1}(beta_j / sqrt(1 + (pi^2 / 3) * sigma_w^2)).
mu_cfo <- stats::qlogis(cfo_skeleton)
beta_cfo <- sapply(mu_cfo, function(muj) {
  stats::rnorm(n_prior, mean = muj, sd = sqrt(cfo_sigma2_beta))
})

sigma2_w_cfo <- 1 / stats::rgamma(n_prior, shape = cfo_eta, rate = cfo_eta)
sigma_w_cfo <- sqrt(sigma2_w_cfo)
denom_cfo <- sqrt(1 + (pi^2 / 3) * sigma_w_cfo^2)

p_cfo <- sweep(beta_cfo, 1, denom_cfo, "/")
p_cfo <- inv_logit(p_cfo)

summary_cfo <- prior_summary(
  pmat = p_cfo,
  model_name = "cfo_pride",
  dose_label = paste0("dose", seq_len(K))
)

## -----------------------------
## Combine results
## -----------------------------
prior_induced_summary <- rbind(
  summary_power,
  summary_alpha,
  summary_ipcrm,
  summary_cfo
)

prior_induced_summary$median <- round(prior_induced_summary$median, 4)
# prior_induced_summary$var    <- round(prior_induced_summary$var, 6)
prior_induced_summary$PESS   <- round(prior_induced_summary$PESS, 2)

print(prior_induced_summary)

## Optional: save
write.csv(
  prior_induced_summary,
  "prior_induced_baseline_median_PESS_four_models.csv",
  row.names = FALSE
)
