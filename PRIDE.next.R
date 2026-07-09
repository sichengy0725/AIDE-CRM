PRIDE.next <- function(target, DLT_res, currdose, ays, ans, iterations, 
                     sigma_prior = c("invgamma", "halfcauchy"), 
                     prior.l = list(mu_k=mu_k, sigma_beta=10, eta=1, tau=1),
                     prior.para=list(alp.prior=target, bet.prior=1-target),
                     cutoff.eli, early.stop){
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
  post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.1){
    alp <- alp.prior + y 
    bet <- bet.prior + n - y
    1 - pbeta(phi, alp, bet)
  }
  
  overdose.fn <- function(phi, threshold, prior.para=list()){
    y <- prior.para$y
    n <- prior.para$n
    alp.prior <- prior.para$alp.prior
    bet.prior <- prior.para$bet.prior
    pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
    # print(data.frame("prob of overdose" = pp))
    if ((pp >= threshold) & (prior.para$n>=3)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  
  gamma.cal <- function(n1, n2, phi, type, alp.prior, bet.prior){
    prob.int <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior){
      alp1 <- alp.prior + y1
      alp2 <- alp.prior + y2
      bet1 <- bet.prior + n1 - y1
      bet2 <- bet.prior + n2 - y2
      fn.min <- function(x){
        dbeta(x, alp1, bet1)*(1-pbeta(x, alp2, bet2)) 
      }
      fn.max <- function(x){
        pbeta(x, alp1, bet1)*dbeta(x, alp2, bet2)
      }
      const.min <- integrate(fn.min, lower=0, upper=0.99, subdivisions=1000, rel.tol = 1e-10)$value
      const.max <- integrate(fn.max, lower=0, upper=1, rel.tol = 1e-10)$value
      p1 <- integrate(fn.min, lower=0, upper=phi)$value/const.min
      p2 <- integrate(fn.max, lower=0, upper=phi)$value/const.max
      
      list(p1=p1, p2=p2)
    }
    
    OR.values <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior, type){
      ps <- prob.int(phi, y1, n1, y2, n2, alp.prior, bet.prior)
      if (type=="L"){
        pC <- 1 - ps$p2
        pL <- 1 - ps$p1
        oddsC <- pC/(1-pC)
        oddsL <- pL/(1-pL)
        OR <- oddsC*oddsL
        
      }else if (type=="R"){
        pC <- 1 - ps$p1
        pR <- 1 - ps$p2
        oddsC <- pC/(1-pC)
        oddsR <- pR/(1-pR)
        OR <- (1/oddsC)/oddsR
      }
      return(OR)
    }
    
    All.OR.table <- function(phi, n1, n2, type, alp.prior, bet.prior){
      ret.mat <- matrix(rep(0, (n1+1)*(n2+1)), nrow=n1+1)
      for (y1cur in 0:n1){
        for (y2cur in 0:n2){
          ret.mat[y1cur+1, y2cur+1] <- OR.values(phi, y1cur, n1, y2cur, n2, alp.prior, bet.prior, type)
        }
      }
      ret.mat
    }
    
    # compute the marginal prob when lower < phiL/phiC/phiR < upper
    # i.e., Pr(Y=y|lower<phi<upper)
    margin.phi <- function(y, n, lower, upper){
      C <- 1/(upper-lower)
      fn <- function(phi) {
        dbinom(y, n, phi)*C
      }
      integrate(fn, lower=lower, upper=upper)$value
    }
    
    # Obtain the table of marginal distribution of (y1, y2) 
    # after intergrate out (phi1, phi2)
    # under H0 and H1
    # H0: phi1=phi, phi < phi2 < 2phi
    # H1: phi2=phi, 0   < phi1 < phi
    margin.ys.table <- function(n1, n2, phi, hyperthesis){
      if (hyperthesis=="H0"){
        p.y1s <- dbinom(0:n1, n1, phi)
        p.y2s <- sapply(0:n2, margin.phi, n=n2, lower=phi, upper=min(1,2*phi))
      }else if (hyperthesis=="H1"){
        p.y1s <- sapply(0:n1, margin.phi, n=n1, lower=0, upper=phi)
        p.y2s <- dbinom(0:n2, n2, phi)
      }
      p.y1s.mat <- matrix(rep(p.y1s, n2+1), nrow=n1+1)
      p.y2s.mat <- matrix(rep(p.y2s, n1+1), nrow=n1+1, byrow=TRUE)
      margin.ys <- p.y1s.mat * p.y2s.mat
      margin.ys
    }
    
    optim.gamma.fn <- function(n1, n2, phi, type, alp.prior, bet.prior){
      OR.table <- All.OR.table(phi, n1, n2, type, alp.prior, bet.prior) 
      ys.table.H0 <- margin.ys.table(n1, n2, phi, "H0")
      ys.table.H1 <- margin.ys.table(n1, n2, phi, "H1")
      
      argidx <- order(OR.table)
      sort.OR.table <- OR.table[argidx]
      sort.ys.table.H0 <- ys.table.H0[argidx]
      sort.ys.table.H1 <- ys.table.H1[argidx]
      n.tol <- length(sort.OR.table)
      
      if (type=="L"){
        errs <- rep(0, n.tol-1)
        for (i in 1:(n.tol-1)){
          err1 <- sum(sort.ys.table.H0[1:i])
          err2 <- sum(sort.ys.table.H1[(i+1):n.tol])
          err <- err1 + err2
          errs[i] <- err
        }
        min.err <- min(errs)
        if (min.err > 1){
          gam <- 0
          min.err <- 1
        }else {
          minidx <- which.min(errs)
          gam <- sort.OR.table[minidx]
        }
      }else if (type=='R'){
        errs <- rep(0, n.tol-1)
        for (i in 1:(n.tol-1)){
          err1 <- sum(sort.ys.table.H1[1:i])
          err2 <- sum(sort.ys.table.H0[(i+1):n.tol])
          err <- err1 + err2
          errs[i] <- err
        }
        min.err <- min(errs)
        if (min.err > 1){
          gam <- 0
          min.err <- 1
        }else {
          minidx <- which.min(errs)
          gam <- sort.OR.table[minidx]
        }
      }
      list(gamma=gam, min.err=min.err)
    }
    
    res <- optim.gamma.fn(n1, n2, phi, type, alp.prior, bet.prior)
    return(res)
  }
  
  p.values <- function(iterations=1000, DLT_res, mu_k, sigma_beta=10, eta=0.1){
    trial_pts <- sum(apply(DLT_res, 1, function(x) any(!is.na(x)))) #total number of patients in trial
    
    # Initialize parameters
    beta <- mu_k  # Initial values for vector beta
    accept_rate_beta <- numeric(length(mu_k))
    
    # Initial values for sigma2
    upper_limit <- qinvgamma(0.95, shape = eta, scale = eta)
    repeat {
      sigma2 <- rinvgamma(1, shape = eta, scale = eta)
      if (sigma2 <= upper_limit) break
    }
    
    W <- rnorm(trial_pts, 0, sqrt(sigma2)) # Initial values for vector W
    accept_rate_W <- numeric(trial_pts)
    
    # Store samples
    beta_samples <- matrix(0, nrow=iterations, ncol=length(mu_k))
    W_samples <- matrix(0, nrow=iterations, ncol=trial_pts)
    sigma2_samples <- numeric(iterations)
    
    log_posterior_beta <- function(beta_k, W, y, mu_k, sigma_beta) {
      log_likelihood <- 0
      for (ii in 1:length(W)) {
        for (jj in 1:length(y[[ii]])) {
          log_likelihood <- log_likelihood + y[[ii]][jj] * (beta_k + W[ii]) - log(1 + exp(beta_k + W[ii]))
        }
      }
      log_prior <- - (beta_k - mu_k)^2 / (2 * sigma_beta^2)
      return(log_likelihood + log_prior)
    }
    
    log_posterior_W <- function(W_i, beta, y, sigma2) {
      log_likelihood <- 0
      for (ii in 1:length(beta)) {
        for (jj in 1:length(y[[ii]])) {
          log_likelihood <- log_likelihood + y[[ii]][jj] * (beta[ii] + W_i) - log(1 + exp(beta[ii] + W_i))
        }
      }
      log_prior <- - W_i^2 / (2 * sigma2)
      return(log_likelihood + log_prior)
    }
    
    proposal_sd_beta <- 1
    proposal_sd_W <- 1
    for (iter in 1:iterations) {
      # update each W_i
      for (i in 1:trial_pts) {
        current_W_i <- W[i]
        proposed_W_i <- rnorm(1, current_W_i, proposal_sd_W)
        index <- which(!is.na(DLT_res[i,]))
        
        log_accept_ratio <- log_posterior_W(proposed_W_i, beta[index], DLT_res[i,index], sigma2) -
          log_posterior_W(current_W_i, beta[index], DLT_res[i,index], sigma2)
        
        if (log(runif(1)) < log_accept_ratio) {
          accept_rate_W[i] <- accept_rate_W[i] + 1
          W[i] <- proposed_W_i
        }
      }
      
      # update sigma2
      shape <- eta + trial_pts / 2
      rate <- eta + sum(W^2) / 2
      sigma2 <- rinvgamma(1, shape = shape, scale = rate)
      
      # update each beta_k
      for (k in 1:length(ans)) {
        current_beta_k <- beta[k]
        proposed_beta_k <- rnorm(1, current_beta_k, proposal_sd_beta)
        index <- which(!is.na(DLT_res[,k]))
        
        if (length(index) != 0) {
          log_accept_ratio <- log_posterior_beta(proposed_beta_k, W[index], DLT_res[index,k], mu_k[k], sigma_beta) -
            log_posterior_beta(current_beta_k, W[index], DLT_res[index,k], mu_k[k], sigma_beta)
          
          if (log(runif(1)) < log_accept_ratio) {
            accept_rate_beta[k] <- accept_rate_beta[k] + 1
            beta[k] <- proposed_beta_k
          }
        } else {
          beta[k] <- rnorm(1, mu_k[k], sigma_beta)
        }
      }
      
      # store samples
      beta_samples[iter, ] <- beta
      #W_samples[iter, ] <- W
      #sigma2_samples[iter] <- sigma2
    }
    
    p_samples <- exp(beta_samples) / (1 + exp(beta_samples))
    #column_means <- colMeans(p_samples, na.rm = TRUE)
    #print(column_means)
    return(p_samples)
  }
  
  p.values.halfcauchy <- function(iterations = 1000, DLT_res, mu_k, sigma_beta = 10, tau = 1) {
    
    trial_pts <- sum(apply(DLT_res, 1, function(x) any(!is.na(x)))) # number of patients
    dose_levels <- ncol(DLT_res)
    
    # Initialize parameters
    beta <- mu_k
    accept_rate_beta <- numeric(length(mu_k))
    
    sigma <- 1  # Initialize sigma (standard deviation)
    W <- rnorm(trial_pts, 0, sigma)
    accept_rate_W <- numeric(trial_pts)
    
    # Store samples
    beta_samples <- matrix(0, nrow = iterations, ncol = length(mu_k))
    W_samples <- matrix(0, nrow = iterations, ncol = trial_pts)
    sigma2_samples <- numeric(iterations)
    
    # Log posterior for beta_k
    log_posterior_beta <- function(beta_k, W, y, mu_k, sigma_beta) {
      log_likelihood <- 0
      for (i in 1:length(W)) {
        for (j in 1:length(y[[i]])) {
          log_likelihood <- log_likelihood + y[[i]][j] * (beta_k + W[i]) - log(1 + exp(beta_k + W[i]))
        }
      }
      log_prior <- - (beta_k - mu_k)^2 / (2 * sigma_beta^2)
      return(log_likelihood + log_prior)
    }
    
    # Log posterior for W_i
    log_posterior_W <- function(W_i, beta, y, sigma) {
      log_likelihood <- 0
      for (j in 1:length(beta)) {
        for (k in 1:length(y[[j]])) {
          log_likelihood <- log_likelihood + y[[j]][k] * (beta[j] + W_i) - log(1 + exp(beta[j] + W_i))
        }
      }
      log_prior <- - W_i^2 / (2 * sigma^2)
      return(log_likelihood + log_prior)
    }
    
    # Log posterior for sigma (Half-Cauchy prior)
    log_posterior_sigma <- function(sigma, W, tau) {
      if (sigma <= 0) return(-Inf)
      n <- length(W)
      log_prior <- log(2 / (pi * tau)) - log(tau^2 + sigma^2)
      log_likelihood <- -n * log(sigma) - sum(W^2) / (2 * sigma^2)
      return(log_prior + log_likelihood)
    }
    
    # Proposal SDs
    proposal_sd_beta <- 1
    proposal_sd_W <- 1
    proposal_sd_sigma <- 0.1  # on log-scale
    
    for (iter in 1:iterations) {
      # Update W_i
      for (i in 1:trial_pts) {
        current_W_i <- W[i]
        proposed_W_i <- rnorm(1, current_W_i, proposal_sd_W)
        index <- which(!is.na(DLT_res[i, ]))
        
        log_ratio <- log_posterior_W(proposed_W_i, beta[index], DLT_res[i, index], sigma) -
          log_posterior_W(current_W_i, beta[index], DLT_res[i, index], sigma)
        
        if (log(runif(1)) < log_ratio) {
          W[i] <- proposed_W_i
          accept_rate_W[i] <- accept_rate_W[i] + 1
        }
      }
      
      # Update sigma using MH (log-normal proposal)
      proposed_sigma <- rlnorm(1, log(sigma), proposal_sd_sigma)
      
      log_ratio_sigma <- log_posterior_sigma(proposed_sigma, W, tau) -
        log_posterior_sigma(sigma, W, tau) +
        log(proposed_sigma) - log(sigma)  # Jacobian adjustment for log-normal proposal
      
      if (log(runif(1)) < log_ratio_sigma) {
        sigma <- proposed_sigma
      }
      
      # Update beta_k
      for (k in 1:dose_levels) {
        current_beta_k <- beta[k]
        proposed_beta_k <- rnorm(1, current_beta_k, proposal_sd_beta)
        index <- which(!is.na(DLT_res[, k]))
        
        if (length(index) > 0) {
          log_ratio <- log_posterior_beta(proposed_beta_k, W[index], DLT_res[index, k], mu_k[k], sigma_beta) -
            log_posterior_beta(current_beta_k, W[index], DLT_res[index, k], mu_k[k], sigma_beta)
          
          if (log(runif(1)) < log_ratio) {
            beta[k] <- proposed_beta_k
            accept_rate_beta[k] <- accept_rate_beta[k] + 1
          }
        } else {
          # If no data for this dose level, sample from prior
          beta[k] <- rnorm(1, mu_k[k], sigma_beta)
        }
      }
      
      # Store samples
      beta_samples[iter, ] <- beta
    }
    
    p_samples <- exp(beta_samples) / (1 + exp(beta_samples))
    #column_means <- colMeans(p_samples, na.rm = TRUE)
    #print(column_means)
    return(p_samples)
  }
  
  OR.values <- function(phi, currdose, type = c("L", "R"), iterations=1000, DLT_res, 
                        mu_k, sigma_beta=10, eta=1, tau=1, sigma_prior = c("invgamma", "halfcauchy")){
    type <- match.arg(type)
    sigma_prior <- match.arg(sigma_prior)
    
    if (sigma_prior == "invgamma") {
      p_samples <- p.values(iterations, DLT_res, mu_k, sigma_beta, eta)
    } else if (sigma_prior == "halfcauchy") {
      p_samples <- p.values.halfcauchy(iterations, DLT_res, mu_k, sigma_beta, tau)
    }
    
    if (type=="L"){
      pC <- mean(p_samples[,currdose] > phi)
      pL <- mean(p_samples[,(currdose-1)] > phi)
      oddsC <- pC/(1-pC)
      oddsL <- pL/(1-pL)
      OR <- oddsC*oddsL
      
    }else if (type=="R"){
      pC <- mean(p_samples[,currdose] > phi)
      pR <- mean(p_samples[,(currdose+1)] > phi)
      oddsC <- pC/(1-pC)
      oddsR <- pR/(1-pR)
      OR <- (1/oddsC)/oddsR
    }
    return(OR)
  }
  
  ###############################################################################
  ############################MAIN DUNCTION######################################
  ###############################################################################
  sigma_prior <- match.arg(sigma_prior)
  ndose <- length(ans)
  if (is.null(prior.para$alp.prior)){
    prior.para <- c(prior.para, list(alp.prior=target, bet.prior=1-target))
  }
  alp.prior <- prior.para$alp.prior
  bet.prior <- prior.para$bet.prior
  
  mu_k <- prior.l$mu_k; sigma_beta<-prior.l$sigma_beta; eta<-prior.l$eta; tau <- prior.l$tau
  
  tover.doses <- rep(0, ndose)
  for (i in 1:ndose){
    cy <- ays[i]
    cn <- ans[i]
    prior.para <- c(list(y=cy, n=cn), list(alp.prior=alp.prior, bet.prior=bet.prior))
    if (overdose.fn(target, cutoff.eli, prior.para)){
      tover.doses[i:ndose] <- 1
      break()
    }
  }
  
  tover.prob <- rep(0, ndose)
  for (i in 1:ndose){
    cy <- ays[i]
    cn <- ans[i]
    tover.prob[i] <- post.prob.fn(target, cy, cn, alp.prior, bet.prior)
  }
  
  if (cutoff.eli != early.stop) {
    cy <- ays[1]
    cn <- ans[1]
    prior.para <- c(list(y=cy, n=cn),list(alp.prior=alp.prior, bet.prior=bet.prior))
    if (overdose.fn(target, early.stop, prior.para)){
      tover.doses[1:ndose] <- 1
    }
  }
  
  if (currdose!=1){
    cys <- ays[(currdose-1):(currdose+1)]
    cns <- ans[(currdose-1):(currdose+1)]
    cover.doses <- tover.doses[(currdose-1):(currdose+1)]
    #cover.doses <- c(0, 0, 0) # No elimination rule
  }else{
    cys <- c(NA, ays[1:(currdose+1)])
    cns <- c(NA, ans[1:(currdose+1)])
    cover.doses <- c(NA, tover.doses[1:(currdose+1)])
    #cover.doses <- c(NA, 0, 0) # No elimination rule
  }
  
  position <- which(tover.doses == 1)[1]
  prior.para <- c(list(alp.prior=alp.prior, bet.prior=bet.prior))
  if ((tover.doses[1] == 1) & (position == 1)){
    index <- NA
    decision <- "stop"
  } else {
    if (cover.doses[2] == 1){
      index <- -1
      decision <- "de-escalation"
    }
    else{
      if (is.na(cys[1]) & (cover.doses[3]==1)){
        index <- 0
        decision <- "stay"
      }
      else  if (is.na(cys[1]) & (!(cover.doses[3]==1))){
        OR.v2 <- OR.values(target, currdose, 'R', iterations, DLT_res, mu_k, sigma_beta, eta,
                           tau, sigma_prior)
        gam2 <- gamma.cal(cns[2], cns[3], target, "R", alp.prior, bet.prior)$gamma
        if (OR.v2>gam2){
          index <- 1
          decision <- "escalation"
        }else{
          index <- 0
          decision <- "stay"
        }
      }
      else  if (is.na(cys[3]) | (cover.doses[3]==1)){
        gam1 <- gamma.cal(cns[1], cns[2], target, "L", alp.prior, bet.prior)$gamma
        OR.v1 <- OR.values(target, currdose, 'L', iterations, DLT_res, mu_k, sigma_beta, eta,
                           tau, sigma_prior)
        if (OR.v1>gam1){
          index <- -1
          decision <- "de-escalation"
        }else{
          index <- 0
          decision <- "stay"
        }
      }
      else  if (!(is.na(cys[1]) | is.na(cys[3]) | cover.doses[3]==1)){
        gam1 <- gamma.cal(cns[1], cns[2], target, "L", alp.prior, bet.prior)$gamma
        gam2 <- gamma.cal(cns[2], cns[3], target, "R", alp.prior, bet.prior)$gamma
        OR.v1 <- OR.values(target, currdose, 'L', iterations, DLT_res, mu_k, sigma_beta, eta,
                           tau, sigma_prior)
        OR.v2 <- OR.values(target, currdose, 'R', iterations, DLT_res, mu_k, sigma_beta, eta,
                           tau, sigma_prior)
        v1 <- OR.v1 > gam1
        v2 <- OR.v2 > gam2
        if (v1 & !v2){
          index <- -1
          decision <- "de-escalation"
        }else if (!v1 & v2){
          index <- 1
          decision <- "escalation"
        }else{
          index <- 0
          decision <- "stay"
        }
      }
    }
  }
  
  if (decision=='stop'){
    nextdose <- 99
  }else{
    nextdose <- currdose+index
  }
  
  out <- list(target=target, DLT_res=DLT_res, decision=decision, currdose = currdose, 
              nextdose=nextdose, overtox=position, toxprob=tover.prob)
  return(out)
}




