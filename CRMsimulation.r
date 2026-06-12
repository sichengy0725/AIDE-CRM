##########################################################################
#                         CRM  Design  simulations                       #
##########################################################################
# -----------------------------------------------------------------------#
#     Main calculation function                                          #
#     PI        --> Ture toxicity                                        #
#     ntrial    --> Number of simulations                                #
#     TARGET    --> Target toxicity                                      #
#     p         --> Initial guess of toxicity probability                #
# -----------------------------------------------------------------------#


# posterior = likelihood x prior
posterior <- function(alpha, p, y, d)
{
  sigma2 = 2;
  lik=1;
  for(i in 1:length(y))
  {
    pi = p[d[i]]^(exp(alpha));
    lik = lik*pi^y[i]*(1-pi)^(1-y[i]);
  }
  return(lik*exp(-0.5*alpha*alpha/sigma2));
}

# used to calculate the posterior mean of pi
posttoxf <- function(alpha, p, y, d, j) { p[j]^(exp(alpha))*posterior(alpha, p, y, d); }

CRM<-function(PI, TARGET=0.3, p, COHORTSIZE=3, ncohort=10, ntrial)
{
  TOXSTOP = 0.9; ## toxicity stopping boundary
  set.seed(6);
  ndose = length(p);
  pi.hat = numeric(ndose); # estimate of toxicity prob
  dose.select=numeric(ndose);  # dose selection
  ntox=rep(0, ndose);  # number of toxicity at each dose level
  ntrted=rep(0, ndose); # number of patient at each dose level
  nstop = 0; # number of trial stopped due to high toxicity 
  
  t.start=Sys.time();
  for(trial in 1:ntrial)
  {
    y=NULL;  #binary outcome
    d=NULL;  #dose level
    dose.curr = 1;  # current dose level
    stop=0; #indicate if trial stops early
    for(i in 1:ncohort)
    {
      # generate data for the new patient
      y = c(y, rbinom(COHORTSIZE, 1, PI[dose.curr]));
      d = c(d, rep(dose.curr, COHORTSIZE));
      
      # calculate posterior mean of toxicity probability at each dose leavel
      marginal=integrate(posterior,lower=-Inf,upper=Inf,p,y,d)$value
      for(j in 1:ndose) { pi.hat[j] = integrate(posttoxf,lower=-Inf,upper=Inf,p,y,d,j)$value/marginal;}
      
      # calculate pr(pi_1>target)
      p.overtox = integrate(posterior,lower=-Inf,upper=log(log(TARGET)/log(p[1])),p,y,d)$value/marginal;	
      if(p.overtox>TOXSTOP) { stop=1; break;}
      
      diff = abs(pi.hat-TARGET);
      dose.best = min(which(diff==min(diff)));
      if(dose.best>dose.curr && dose.curr != ndose) dose.curr = dose.curr+1;
      if(dose.best<dose.curr && dose.curr != 1) dose.curr = dose.curr-1;
    }
    if(stop==1) { nstop=nstop+1; }
    else { dose.select[dose.best] = dose.select[dose.best]+1; }
    # record summary statistics
    for(j in 1:ndose) { ntox[j] = ntox[j] + sum(y[d==j]); }
    for(j in 1:ndose) { ntrted[j] = ntrted[j] + sum(d==j); }
  }
  t.end=Sys.time();
  cat("selection probability: ", dose.select/ntrial*100, "\n");
  cat("number of toxicity:    ", c(ntox/ntrial, sum(ntox[PI>TARGET]/ntrial)), "\n");
  cat("number of treated:     ", c(ntrted/ntrial, sum(ntrted[PI>TARGET]/ntrial)), "\n");
  cat("percentage of stop:    ", nstop/ntrial*100, "\n");
  cat("duration:  ", ncohort*T, "\n");
  print(t.end-t.start);
}


## Your 5-dose scenarios
scenarios <- rbind(
  SC1 = c(0.07, 0.12, 0.17, 0.22, 0.30),
  SC2 = c(0.05, 0.10, 0.18, 0.30, 0.40),
  SC3 = c(0.15, 0.20, 0.30, 0.35, 0.45),
  SC4 = c(0.15, 0.30, 0.38, 0.45, 0.55),
  SC5 = c(0.30, 0.35, 0.40, 0.45, 0.50),
  SC6 = c(0.50, 0.55, 0.60, 0.65, 0.70)
)

## CRM skeleton used in AIDE paper for 5-dose setting
p_skeleton <- c(0.15, 0.20, 0.30, 0.35, 0.45)

for (s in seq_len(nrow(scenarios))) {
  cat("\n====================================\n")
  cat("Scenario", s, "\n")
  cat("True PI:", paste(scenarios[s, ], collapse = ", "), "\n")
  cat("Skeleton:", paste(p_skeleton, collapse = ", "), "\n")
  cat("====================================\n")
  
  CRM(
    PI = scenarios[s, ],
    TARGET = 0.3,
    p = p_skeleton,
    COHORTSIZE = 3,
    ncohort = 10,
    ntrial = 2000
  )
}

