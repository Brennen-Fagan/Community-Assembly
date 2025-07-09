# following paulbuerkner.com/brms/articles/brms_customfamilies.html
# and paulbuerkner.com/brms/articles/brms_customfamily.html
library(brms)

skellam <- brms::custom_family(
  "Skellam", dpars = c("mu", "sigma"),
  links = c("identity", "log"),
  lb = c(NA, 0), ub = c(NA, NA),
  type = "int"
)

# lpmf => logged
stan_funs <- "
  real Skellam_lpmf(int k, real mu, real sigma) {
    real kk = k * 1.0;
    real lambda1 = (mu + sigma^2)/2;
    real lambda2 = (-mu + sigma^2)/2;
    real bes = log(modified_bessel_first_kind(k, 2 * sqrt(lambda1 * lambda2)));
    return((-lambda1-lambda2) + (kk/2) * (log(lambda1) - log(lambda2)) + bes);
  }
  int Skellam_rng(real mu, real sigma) {
    real lambda1 = (mu + sigma^2)/2;
    real lambda2 = (-mu + sigma^2)/2;
    return(poisson_rng(lambda1) - poisson_rng(lambda2));
  }
"

stanvars <- brms::stanvar(scode = stan_funs, block = "functions")

# test:
testSkellam <- data.frame(
  p1 = rpois(100, 20),
  p2 = rpois(100, 5)
)
testSkellam$s <- with(testSkellam, p1 - p2)
brms::brm(s ~ 1, family = skellam, stanvars = stanvars, data = testSkellam)

load("/shared/storage/biology/rsrch/lcab/research/btf500/CommunityAssembly/To Github/Community-Assembly/experiments/VariantExperiments/diversitiesFlattened9_142486-4929_359-359_2025-01-22.RData")
library(tidytable)
diversitiesFlattened %>% filter(is.na(Subset), NicheDistance == "5", Intervention %in% c("(0)", "(0.5)->(0)"), SpeciesAffinity == "runif", Metric == "Alpha Hill:0", Time >= 20000, Time <= 30000) %>% tidytable::group_by(Time) %>% tidytable::mutate(Value = ifelse(is.na(InterventionPatchType), -Value, Value)) %>% tidytable::summarise(Value = round(sum(Value))) -> tempdat
brms::brm(Value ~ 1 + arma(p = 1, q = 1, cov = TRUE), family = skellam, stanvars = stanvars, data = tempdat,
          cores = 8, init = 0, chains = 8, control = list(adapt_delta = 0.95), warmup = 2000, iter = 3000) -> temp
# cov must be TRUE because we're using a custom distribution, and
# cov = TRUE means p <= 1, q <= 1...
diversitiesFlattened %>% filter(is.na(Subset), NicheDistance == "5", Intervention %in% c("(0)"), SpeciesAffinity == "runif", Metric == "Alpha Hill:0", Time >= 20000, Time <= 30000) %>% mutate(Value = round(Value)) -> tempdatNB
brms::brm(Value ~ 1 + arma(p = 1, q = 1, cov = TRUE), family = negbinomial, stanvars = stanvars, data = tempdatNB,
          cores = 8, init = 0, chains = 8, control = list(adapt_delta = 0.95), warmup = 2000, iter = 3000) -> tempNB


library(rstan)
modelstan_mixNB <- "
data {
  int<lower=1> K; // mixture components
  int<lower=1> N; // data points
  int y[N]; // observations
}
parameters {
  simplex[K] theta; // mixture proportions
  positive_ordered[K] mu; // mixture means, ordered for identifiability
  vector<lower=0>[K] phi; // overdispersion parameters
}
model {
  vector[K] log_theta = log(theta); // recommended to cache log calculation
  // priors
  mu ~ lognormal(2.7, 0.3); // Median ~15, Std.Dev ~10
  phi ~ lognormal(1, 1); // Median ~2, Variance ~20
  for (n in 1:N) {
    vector[K] lps = log_theta;
    for (k in 1:K) {
      lps[k] += neg_binomial_2_lpmf(y[n] | mu[k], phi[k]);
    }
    target += log_sum_exp(lps);
  }
}
"
fit1 <- rstan::stan(
  model_code = modelstan_mixNB,
  data = list(
    K = 4,
    N = nrow(tempdatNB),
    y = round(tempdatNB$Value) # make sure no precision errors.
  ),
  chains = 6,
  warmup = 3000,
  iter = 4000,
  cores = 6
)

modelstan_mixarmaNB <- "
data {
  int<lower=1> K; // mixture components
  int<lower=1> N; // data points
  int y[N]; // observations
}
parameters {
  simplex[K] theta; // mixture proportions
  positive_ordered[K] mu; // mixture means, ordered for identifiability
  vector<lower=0>[K] phi; // overdispersion parameters
  real beta; // ar(1)
  real gamma; // ma(1)
}
transformed parameters {
  real epsilon[N, K]; // error terms, note now deprecated syntax (but 2.21)
  for (k in 1:K) {
    epsilon[1, k] = y[1]*1.0 - mu[k];
  }
  for (n in 2:N) {
    for (k in 1:K) {
      epsilon[n, k] = ( y[n]*1.0 - mu[k] - gamma * epsilon[n - 1, k] );
    }
  }
}
model {
  vector[K] log_theta = log(theta); // recommended to cache log calculation
  // priors
  mu ~ lognormal(2.7, 0.3); // Median ~15, Std.Dev ~10
  phi ~ lognormal(1, 1); // Median ~2, Variance ~20
  beta ~ normal(0, 0.75); // Probably should be |.| < 1 most of the time
  gamma ~ normal(0, 0.75); //
  for (n in 2:N) {
    vector[K] lps = log_theta;
    for (k in 1:K) {
      lps[k] += (
        neg_binomial_2_lpmf(y[n] | mu[k] + beta * (y[n-1] - mu[k]) + gamma * epsilon[n-1, k], phi[k])
      );
    }
    target += log_sum_exp(lps);
  }
}
"
fit2 <- rstan::stan(
  model_code = modelstan_mixarmaNB,
  data = list(
    K = 4,
    N = nrow(tempdatNB),
    y = round(tempdatNB$Value) # make sure no precision errors.
  ),
  chains = 6,
  warmup = 3000,
  iter = 4000,
  cores = 6
)
