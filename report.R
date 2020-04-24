library(ggplot2)
library(dplyr)
library(jsonlite)
library(cmdstanr)
library(posterior)
library(glue)
library(microbenchmark)
library(rhdf5)

set_cmdstan_path("~/cmdstan")
theme_set(theme_minimal())

branches <- c("develop", "proposal")
metrics <- c("dense", "diag")

dials <- tibble(cores = 4, chains = 4,
                warmup = 5000, samples = 2 * warmup,
                refresh = 2000, evals = 1)

fit_model <- function(mdl, metric, dial = dials) {
    dial$seed <- sample.int(10^5, 1)
    mod <- cmdstan_model(glue("{mdl}/{mdl}.stan"))
    data <- read_json(glue("{mdl}/{mdl}.json"), simplifyVector = TRUE)
    time <- microbenchmark("runtime" =
                               fit <- mod$sample(data = data,
                                                 num_warmup = dial$warmup,
                                                 num_samples = dial$samples,
                                                 num_chains = dial$chains,
                                                 num_cores = dial$cores,
                                                 metric = metric, refresh = dial$refresh),
                           times = dial$evals)
    list(time = time, fit = fit)
}



store_results <- function(output, branch, metric, model) {
    h5write(max(output$time$time / 1e9),
            "results.h5",
            glue("{branch}/{metric}/{model}_time")) # seconds
    h5write(output$fit$sampling_info()$seed,
            "results.h5",
            glue("{branch}/{metric}/{model}_seed"))

    draws <- output$fit$draws()
    h5write(draws,
            "results.h5",
            glue("{branch}/{metric}/{model}_draws"))
    h5write(summarise_draws(draws),
            "results.h5",
            glue("{branch}/{metric}/{model}_summary"))
    h5write(summarise_draws(draws^2),
            "results.h5",
            glue("{branch}/{metric}/{model}_summary2"))
}

models <- c(## "gauss_32d", "gauss_c25_32d", "gauss_c25_64d", "gauss_c25_128d",
            ## "gauss_c50_32d", "gauss_c50_64d", "gauss_c50_128d",
            ## "gauss_c75_32d", "gauss_c75_64d", "gauss_c75_128d",  "cauchy",
    "arK", "arma", "eight_schools", "garch", "gp_pois_regr", "gp_regr", "irt_2pl",
    "low_dim_gauss_mix", "low_dim_gauss_mix_collapse", #"pkpd",
    "sir")

dials$chains <- 1
dials$cores <- 1
dials$warmup <- 1000
dials$samples <- 1000

capture.output(out <- fit_model("arma", "dense_e", dial=dials))

dfs <- rep(NA, length(models))

for (m in 21) {
    ## essbulk_models <- matrix(NA, nrow = length(models), ncol=dials$samples)
    ## essbulk2_models <- matrix(NA, nrow = length(models), ncol=dials$samples)
    ## esstail_models <- matrix(NA, nrow = length(models), ncol=dials$samples)
    ## esstail2_models <- matrix(NA, nrow = length(models), ncol=dials$samples)
    ## rhat_models <- matrix(NA, nrow=)
    out <- fit_model(models[m], "diag_e", dial = dials)
    dfs[m] <- nrow(out$fit$summary())
}

## > dfs
##  [1]  33  33  65 129  33  65 129  33  65 129   8   5   4  19   5  25   4 145   6
## [20]   6  NA  85


## fh5 <- h5createFile("results.h5")


## for (b in branches) {
##     h5createGroup("results.h5", b)
##     for (m in metrics) {
##         h5createGroup("results.h5", glue("{b}/{m}"))
##         for (model in models) {
##             h5createGroup("results.h5", glue("{b}/{m}/{model}"))
##         }
##     }
## }
