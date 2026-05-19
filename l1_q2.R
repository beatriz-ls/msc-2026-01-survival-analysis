# =========================================================
# Modelos de Sobrevivência no Stan
# Banco: sobrevivencia::isoladores
# =========================================================

# =========================================================
# Pacotes
# =========================================================

library(rstan)
library(sobrevivencia)
library(survival)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# =========================================================
# Banco de dados
# =========================================================

data(isoladores, package = "sobrevivencia")

head(isoladores)

# =========================================================
# Preparação dos dados
# =========================================================

N <- nrow(isoladores)

stan_data <- list(
  N = N,
  t = isoladores$tempo,
  d = isoladores$status
)

# =========================================================
# MODELO EXPONENCIAL
# =========================================================

exp_stan <- '
data {
  int<lower=1> N;
  vector<lower=0>[N] t;
  array[N] int<lower=0,upper=1> d;
}

parameters {
  real<lower=0> lambda;
}

model {

  lambda ~ gamma(0.01, 0.01);

  for (i in 1:N) {

    if (d[i] == 1)
      target += exponential_lpdf(t[i] | lambda);

    else
      target += exponential_lccdf(t[i] | lambda);

  }

}
'

mod_exp <- stan_model(model_code = exp_stan)

fit_exp <- sampling(
  mod_exp,
  data = stan_data,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  seed = 123
)

print(fit_exp)

# =========================================================
# MODELO GAMMA
# =========================================================

gamma_stan <- '
data {
  int<lower=1> N;
  vector<lower=0>[N] t;
  array[N] int<lower=0,upper=1> d;
}

parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
}

model {

  alpha ~ gamma(0.01,0.01);
  beta ~ gamma(0.01,0.01);

  for (i in 1:N) {

    if (d[i] == 1)
      target += gamma_lpdf(t[i] | alpha, beta);

    else
      target += gamma_lccdf(t[i] | alpha, beta);

  }

}
'

mod_gamma <- stan_model(model_code = gamma_stan)

fit_gamma <- sampling(
  mod_gamma,
  data = stan_data,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  seed = 123
)

print(fit_gamma)

# =========================================================
# MODELO LOG-NORMAL
# =========================================================

lnorm_stan <- '
data {
  int<lower=1> N;
  vector<lower=0>[N] t;
  array[N] int<lower=0,upper=1> d;
}

parameters {
  real mu;
  real<lower=0> sigma;
}

model {

  mu ~ normal(0,100);
  sigma ~ cauchy(0,5);

  for (i in 1:N) {

    if (d[i] == 1)
      target += lognormal_lpdf(t[i] | mu, sigma);

    else
      target += lognormal_lccdf(t[i] | mu, sigma);

  }

}
'

mod_lnorm <- stan_model(model_code = lnorm_stan)

fit_lnorm <- sampling(
  mod_lnorm,
  data = stan_data,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  seed = 123
)

print(fit_lnorm)

# =========================================================
# MODELO LOG-LOGÍSTICO
# =========================================================

llogis_stan <- '
functions {

  real loglogistic_lpdf(real t, real alpha, real beta) {

    return log(beta) - log(alpha) +
           (beta - 1) * (log(t) - log(alpha)) -
           2 * log1p(pow(t / alpha, beta));

  }

  real loglogistic_lccdf(real t, real alpha, real beta) {

    return -log1p(pow(t / alpha, beta));

  }

}

data {

  int<lower=1> N;
  vector<lower=0>[N] t;
  array[N] int<lower=0,upper=1> d;

}

parameters {

  real<lower=0> alpha;
  real<lower=0> beta;

}

model {

  alpha ~ gamma(0.01,0.01);
  beta ~ gamma(0.01,0.01);

  for (i in 1:N) {

    if (d[i] == 1)
      target += loglogistic_lpdf(t[i] | alpha, beta);

    else
      target += loglogistic_lccdf(t[i] | alpha, beta);

  }

}
'

mod_llogis <- stan_model(model_code = llogis_stan)

fit_llogis <- sampling(
  mod_llogis,
  data = stan_data,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  seed = 123
)

print(fit_llogis)
