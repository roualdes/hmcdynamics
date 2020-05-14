source("utilities.R")
set_cmdstan_path("~/cmdstan")
theme_set(theme_minimal())
replications <- 25
throwaway <- 5

## simulation
branch <- "proposal"
metric <- "dense_e"

## dim <- 8; rho <- 0; nu <- 3; rep <- 10; model <- "gaussian"

for (model in c("gaussian", "student_t")) {
    force_recompile("~/hmcdynamics", model)
    mod <- cmdstan_model(glue("{model}.stan"))

    for (dim in 2^(1:6)) {
        for (rho in seq(0, .75, by = 0.25)) {
            for (rep in 1:replications) {
                seed <- use_seed(model, branch, dim, rho, rep)

                fit <- mod$sample(data = gendata(dim, rho), num_cores = 4,
                                  num_warmup = 10000, num_samples = 10000,
                                  seed = seed, metric = metric, refresh = 10000)

                ## discard first throwaway runs
                if (rep > throwaway) {
                    write_output(fit, model, branch, dim, rho, rep - throwaway)
                }
            }
        }
    }
}

## pre analysis
rhos <- seq(0, 0.75, by=0.25)
models <- c("gaussian", "student_t")
branches <- c("proposal", "develop")
dims <- 2^(1:6)
reps <- 1:20

for (rho in rhos) {
    for (model in models) {
    df <- prepare_df(model, rho, branches, dims, reps)
    write_csv(df, glue("output/{model}_{rho}.csv.gz"))
    }
}

## analysis
rho <- 0
model <- "gaussian"
df <- read.csv(glue("output/{model}_{rho}.csv.gz"))

plot_ess_time(df, model, rho)
table_ess_time(df)

plot_ess_leapfrog(df, model, rho)
table_ess_leapfrog(df)

plot_stats(df, model, rho)
plot_rhat(df, model, rho)

check_stepsize(df)
check_leapfrog(df)

check_runtime(df)
plot_runtime(df)

check_divergence(df)
check_seed(df)
