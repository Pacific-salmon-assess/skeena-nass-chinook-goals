// Stan version of age-structured state-space spawner-recruitment model with AR-1 process variation (adapted from Fleischman et al. CJFAS. 2013)
data{
  int nyrs;           // number of calender years
  int a_min;          // minimum age class
  int a_max;          // maximum age class
  int A;              // number of age classes
  int nRyrs;          // number of recruitment years
  int A_obs[nyrs, A]; // observed age composition in counts by age class
  vector[nyrs] S_obs; // observed spawners
  vector[nyrs] H_obs; // observed harvest
  vector[nyrs] S_cv;  // spawner observation error CV
  vector[nyrs] H_cv;  // harvest observation error CV
}
parameters{
  real<lower=0,upper=1> F_rw; //fraction of variance attributed to random walk in alpha
  real<lower=0> sigma_tot;
  real ln_alpha0;                         // initial productivity (on log scale)
  vector[nyrs] alpha_dev;                // time varying year-to-year deviations in a
  real<lower=0> beta;                     // Ricker b
  vector<lower=0>[nRyrs] lnR;             // log recruitment states
  real<lower=0> sigma_R0;                 // process error for first a.max years with no spawner link
  real lnresid_0;                         // first residual for lag-1 correlation in process error
  real<lower=0> mean_ln_R0;               // "true" mean log recruitment in first a.max years with no spawner link
  vector<lower=0.01,upper=0.99>[nyrs] U;  // harvest rate
  vector<lower=0,upper=1>[3] prob;        // maturity schedule probs
  real<lower=0,upper=1> D_scale;          // governs variability of age proportion vectors across cohorts
  matrix<lower=0.01>[nRyrs, A] g;         // individual year/age class gamma variates for generating age at maturity proportions
}
transformed parameters{
  vector<lower=0>[nyrs] N;              // run size states
  vector<lower=0>[nyrs] S;              // spawner states
  vector[nyrs] lnS;                     // log spawner states
  vector<lower=0>[nyrs] C;              // catch states
  vector[nyrs] lnC;                     // log catch states
  vector<lower=0>[nRyrs] R;             // recruitment states
  vector[nRyrs] lnresid;                // log recruitment residuals
  vector[nRyrs] lnRm_1;                 // log recruitment states in absence of lag-one correlation //DO I NEED THIS
  matrix<lower=0>[nyrs, A] N_ta;        // returns by age matrix
  matrix<lower=0, upper=1>[nRyrs, A] p; // age at maturity proportions
  vector<lower=0,upper=1>[A] pi;        // maturity schedule probs
  real<lower=0> D_sum;                  // inverse of D_scale which governs variability of age proportion vectors across cohorts
  vector<lower=0>[A] Dir_alpha;         // Dirichlet shape parameter for gamma distribution used to generate vector of age-at-maturity proportions
  matrix<lower=0, upper=1>[nyrs, A] q;  // age composition by year/age classr matrix
  vector[nyrs] ln_alpha;               // ln_alpha in each year
  real<lower=0> sigma_alpha;              // amount alpha can deviate
  real<lower=0> sigma_R;                  // process error

  // Maturity schedule: use a common maturation schedule to draw the brood year specific schedules
  pi[1] = prob[1];
  pi[2] = prob[2] * (1 - pi[1]);
  pi[3] = 1 - pi[1] - pi[2];
  D_sum = 1/D_scale^2;
  
  for (a in 1:A) {
    Dir_alpha[a] = D_sum * pi[a];
    for (y in 1:nRyrs) {
      p[y,a] = g[y,a]/sum(g[y,]);
    }
  }
  
  sigma_alpha = F_rw*sigma_tot; //process variance - fraction of total attributed to random walk
  sigma_R = (1-F_rw)*sigma_tot; //obs variance - remaining fraction of total variance
  
  R = exp(lnR);
  // Calculate the numbers at age matrix as brood year recruits at age (proportion that matured that year)
  for (t in 1:nyrs) {
    for(a in 1:A){
      N_ta[t,a] = R[t+A-a] * p[t+A-a,a];
    }
  }
  // Calculate returns, spawners and catch by return year
  for(t in 1:nyrs) {
    N[t] = sum(N_ta[t,1:A]);
    S[t] = N[t] * (1 - U[t]);
    lnS[t] = log(S[t]);
    C[t] = N[t] * U[t];
    lnC[t] = log(C[t]);
  }
  // Calculate age proportions by return year
  for (t in 1:nyrs) {
    for(a in 1:A){
      q[t,a] = N_ta[t,a]/N[t];
    }
  }
  
  //time-varying alpha component
  ln_alpha[1] = ln_alpha0; //first gets the prior

  for(t in 2:nyrs){
    ln_alpha[t] = ln_alpha[t-1] + alpha_dev[t-1]*sigma_alpha; // random walk of log_a
  }
  
  // Ricker SR with with time-varying productivity for years with brood year spawners  
  for (i in 1:a_max) { //unobserved resids
    lnresid[i] = 0.0;
    lnRm_1[i] = 0.0;
  }
  for(y in (A+a_min):nRyrs){ //8:42
    lnRm_1[y] = lnS[y-a_max] + ln_alpha[y-a_max] - beta * S[y-a_max]; //ln_alpha[y-a_max] links it to the spawners
    lnresid[y] = lnR[y] - lnRm_1[y];
  }
}
model{
  // Priors
  sigma_tot ~ gamma(2,1);
  F_rw ~ beta(2,5); // HH changed from beta(2,4)
  ln_alpha0 ~ normal(1,4);
  alpha_dev ~ std_normal();                    //standardized (z-scales) deviances
  //sigma_alpha ~ normal(0,1);
  beta ~ normal(0,1);
  //sigma_R ~ normal(0,2);
  lnresid_0 ~ normal(0,20);
  mean_ln_R0 ~ normal(0,20);
  sigma_R0 ~ gamma(2,1); 
  prob ~ beta(1,1);
  D_scale ~ beta(1,1);

  // Likelihoods
  // Gamma variates for each year and age class which are used to determine age at maturity proportions
  for (y in 1:nRyrs) {
    for (a in 1:A) {
      //g[y,a] ~ gamma(Dir_alpha[a],1);
      target += gamma_lpdf(g[y,a]|Dir_alpha[a], 1);
    }
  }

  // Process model
  lnR[1:a_max] ~ normal(mean_ln_R0, sigma_R0);
  lnR[(A+a_min):nRyrs] ~ normal(lnRm_1[(A+a_min):nRyrs], sigma_R);

  // Observation model
  for(t in 1:nyrs){
  //A_obs[t,1:A]) ~ multinomial(q[t,1:A]);
    target += multinomial_lpmf(A_obs[t,1:A]|to_vector(q[t,1:A]));
    U[t] ~ beta(1,1);
    H_obs[t] ~ lognormal(lnC[t], sqrt(log((H_cv[t]^2) + 1)));
    S_obs[t] ~ lognormal(lnS[t], sqrt(log((S_cv[t]^2) + 1)));
  }
}

generated quantities {
  //benchmarks
  
  //prior and posterior predictive check
  //below is what I did for fraser pink, but I wonder if I should:
    //-add A_obs and 
    //-makle the CVs rng too
  array[nyrs] real H_rep;
  array[nyrs] real S_rep;

  H_rep = lognormal_rng(lnC, sqrt(log((H_cv^2) + 1))); 
  S_rep = lognormal_rng(lnS, sqrt(log((S_cv^2) + 1)));
}
