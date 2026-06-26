## ------------------------------------------------------------
## Run AIDE-CRM OC simulation for 100 trials
## ------------------------------------------------------------

rm(list = ls())

source("AIDE_BOIN_helper.R")
source("AIDE_CRM_helper_modified.R")
source("AIDE_modified.R")

set.seed(20260612)

TARGET <- 0.30

p.true <- c(0.07, 0.12, 0.17, 0.22, 0.30)

r_carry <- 0

## Since r_carry = 0, IPDE truth is the same as regular-patient truth
p.true_ipde <- p.true
## equivalently:
## p.true_ipde <- r_carry + (1 - r_carry) * p.true

crm_skeleton <- c(0.15, 0.20, 0.30, 0.35, 0.45)

fit100 <- get_oc_sim_AIDE(
  target = TARGET,
  p.true = p.true,
  p.true_ipde = p.true_ipde,
  
  ntrial = 1,
  seed = 9,
  
  model = "BOIN",
  ipde_design = 2,
  
  N_pat = 30L,
  Nmax_eff = 30L,
  C = 3L,
  T_assess = 28,
  cycle_max = 2L,
  
  arrival_rate = 1 / 28,
  t0 = 0,
  
  cutoff = 0.95,
  
  d.cap = 100,
  dose_cap = 3L,
  day_obs = 0,
  
  r_carry = r_carry,
  r_estimator = "r_mle",
  crm_r_model = "level",
  decision_method = "approx2",
  mtd_method = "approx2",
  crm_model_file = "random_CRM_level.bug",
  crm_skeleton = crm_skeleton,
  crm_alpha_sd = 2,
  crm_a_r = 1,
  crm_b_r = 9,
  crm_n_chains = 2,
  crm_n_adapt = 500,
  crm_n_burnin = 500,
  crm_n_iter = 2000,
  crm_thin = 1,
  
  store_raw = TRUE,
  verbose = TRUE
)

