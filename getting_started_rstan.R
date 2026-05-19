# run the next line if you already have rstan installed
# remove.packages(c("StanHeaders", "rstan"))

install.packages("rstan", repos = c('https://stan-dev.r-universe.dev',
                                    getOption("repos")))

# Rodar isso primerio:
pkgbuild::has_build_tools(debug = TRUE)

# calmar biblioteca
library(rstan)

# To verify your installation, the RStan example/test model:

example(stan_model, package = "rstan", run.dontrun = TRUE)



rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
