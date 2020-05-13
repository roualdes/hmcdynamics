source("utilities.R")
set_cmdstan_path("~/cmdstan")
theme_set(theme_minimal())
replications <- 10

## simulation
branch <- "develop"

metric <- "dense_e"

## dim <- 8; rho <- 0; nu <- 3; rep <- 10; model <- "gaussian"

for (model in c("gaussian", "student_t")) {
    force_recompile("~/hmcdynamics", model)
    mod <- cmdstan_model(glue("{model}.stan"))

    for (dim in 2^(1:8)) {
        for (rho in seq(0, .75, by = 0.25)) {
            for (rep in 1:replications) {
                seed <- use_seed(model, branch, dim, rho, rep)

                fit <- mod$sample(data = gendata(dim, rho), num_cores = 4,
                                  seed = seed, metric = metric, refresh = 2000)

                write_output(fit, model, branch, dim, rho, rep)
            }
        }
    }
}


## analysis

rhos <- seq(0, .75, by = 0.25)
model <- "student_t"

branches <- c("proposal", "develop")
dims <- 2^(1:8)

df_ess <- expand.grid(branch = branches, dim = dims, rep = 1:replications)
df_ess$ess_bulk_sec <- NA
df_ess$ess_tail_sec <- NA
df_ess$ess2_bulk_sec <- NA
df_ess$ess2_tail_sec <- NA
df_ess$time <- NA

rho <- 0.0
for (b in branches) {
    for (d in dims) {
        for (r in 1:replications) {
            dfs <- read_csv(genpath(model, b, d, rho, r, "summary"))
            dft <- read_csv(genpath(model, b, d, rho, r, "time"))
            ess2 <- as.matrix(read_csv(genpath(model, b, d, rho, r, "draws")))^2 %>%
                array(dim=c(1000, 4, 10)) %>%
                summarise_draws

            idx <- with(df_ess,
                        which(branch == b & dim == d & rep == r))

            time <- sum(dft$total)
            df_ess$ess_bulk_sec[idx] <- min(dfs$ess_bulk / time)
            df_ess$ess_tail_sec[idx] <- min(dfs$ess_tail / time)
            df_ess$ess2_bulk_sec[idx] <- min(ess2$ess_bulk / time)
            df_ess$ess2_tail_sec[idx] <- min(ess2$ess_tail / time)
            df_ess$time[idx] <- time
        }
    }
}

df_ess %>%
    group_by(dim, branch) %>%
    summarise(essb = mean(ess_bulk_sec),
              esst = mean(ess_tail_sec),
              time = mean(time))

df_ess %>%
    ggplot(aes(factor(dim), log10(ess2_bulk_sec), color=branch)) +
    geom_point(position = position_jitterdodge()) +
    stat_summary(fun.data = "median_hilow", position = position_jitterdodge())

df_ess %>%
    ggplot(aes(time, color=branch)) +
    geom_density() +
    facet_wrap(~dim, scales="free")
