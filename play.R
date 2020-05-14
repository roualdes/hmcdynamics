source("utilities.R")
set_cmdstan_path("~/cmdstan")
theme_set(theme_minimal())
replications <- 25
throwaway <- 5

## simulation
branch <- "develop"
metric <- "dense_e"

## dim <- 8; rho <- 0; nu <- 3; rep <- 10; model <- "german_credit"

for (model in c("gaussian", "student_t")) {

    force_recompile("~/hmcdynamics", model)
    mod <- cmdstan_model(glue("{model}.stan"))

    for (dim in 2^(1:6)) {
        for (rho in seq(0, .75, by = 0.25)) {
            for (rep in 1:replications) {
                ## TODO something wrong with this seeding idea.
                seed <- use_seed(model, branch, dim, rho, rep)

                fit <- mod$sample(data = gendata(model, dim, rho), num_cores = 4,
                                  num_warmup = 5000, num_samples = 5000, adapt_delta = 0.95,
                                  seed = seed, metric = metric, refresh = 5000)

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
rho <- 0.75
model <- "student_t"
df <- read.csv(glue("output/{model}_{rho}.csv.gz"))

## ess
## time
df %>%
    select(branch, dim, ess_bulk_sec, ess_tail_sec, ess2_bulk_sec, ess2_tail_sec) %>%
    pivot_longer(cols = c(ess_bulk_sec, ess_tail_sec, ess2_bulk_sec, ess2_tail_sec)) %>%
    ggplot(aes(factor(dim), value, color = branch, shape = branch)) +
    geom_point(position = position_jitterdodge(), alpha = 0.25) +
    scale_y_log10() +
    stat_summary(fun.data = "median_hilow", position = position_jitterdodge()) +
    facet_wrap(~name, nrow=2, scales="free_y") +
    labs(title = bquote("Model" ~ .(model) ~ "with rho =" ~ .(rho)))

df %>%
    group_by(branch, dim) %>%
    summarise(mean_ess_bulk = mean(ess_bulk_sec),
              mean_ess_tail = mean(ess_tail_sec),
              mean_ess2_bulk = mean(ess2_bulk_sec),
              mean_ess2_tail = mean(ess2_tail_sec)) %>%
    arrange(dim)

## leapfrog
df %>%
    select(branch, dim, ess_bulk_leapfrog, ess2_bulk_leapfrog, ess_tail_leapfrog, ess2_tail_leapfrog) %>%
    pivot_longer(cols = c(ess_bulk_leapfrog, ess2_bulk_leapfrog, ess_tail_leapfrog, ess2_tail_leapfrog)) %>%
    ggplot(aes(factor(dim), value, color = branch, shape = branch)) +
    geom_point(position = position_jitterdodge(), alpha = 0.25) +
    stat_summary(fun.data = "median_hilow", position = position_jitterdodge()) +
    facet_wrap(~name, nrow=2, scales="free_y") +
    labs(title = bquote("Model" ~ .(model) ~ "with rho =" ~ .(rho)))

## stats
df %>%
    select(branch, dim, q5, median, q95, sd) %>%
    pivot_longer(cols = c(q5, median, q95, sd)) %>%
    ggplot(aes(factor(dim), value, color=branch)) +
    geom_point(position = position_jitterdodge(), alpha = 0.25) +
    stat_summary(fun.data = "median_hilow", position = position_jitterdodge()) +
    facet_wrap(~name, nrow=2, scale="free_y") +
    labs(title = bquote("Model" ~ .(model) ~ "with rho =" ~ .(rho)))

## rhat
df %>%
    ggplot(aes(factor(dim), rhat, color = branch, shape = branch)) +
    geom_point(position = position_jitterdodge(), alpha = 0.25) +
    stat_summary(fun.data = "median_hilow", position = position_jitterdodge()) +
    geom_hline(yintercept = 1.01, linetype = 2) +
    labs(x = "Dimension", y = "Rhat", title = bquote("Model" ~ .(model) ~ "with rho =" ~ .(rho)))

## stepsize
df %>%
    group_by(branch, dim) %>%
    summarise(m_eps = mean(stepsize)) %>%
    arrange(dim)

## leapfrog
df %>%
    group_by(branch, dim) %>%
    summarise(m_leapfrog = mean(leapfrog)) %>%
    arrange(dim)

## run time
df %>%
    group_by(branch, dim) %>%
    summarise(m_time = mean(time)) %>%
    arrange(dim)

df %>%
    ggplot(aes(factor(dim), time, color=branch)) +
    geom_point(position = position_jitterdodge(), alpha = 0.25) +
    stat_summary(fun.data = "median_hilow", position = position_jitterdodge())


## divergences
df %>%
    group_by(branch, dim) %>%
    summarise(m_divergent = mean(divergent)) %>%
    arrange(dim)

## seeds
df %>%
    group_by(branch, dim) %>%
    summarise(m_seed = mean(seed)) %>%
    arrange(dim)
