 models <- c("gauss_32d", "gauss_c25_32d", "gauss_c25_64d", "gauss_c25_128d",
             "gauss_c50_32d", "gauss_c50_64d", "gauss_c50_128d",
             "gauss_c75_32d", "gauss_c75_64d", "gauss_c75_128d",
             "cauchy", "arK", "eight_schools", "garch", "gp_pois_regr",
             "gp_regr", "irt_2pl", "low_dim_gauss_mix", "low_dim_gauss_mix_collapse",
             "pkpd", "sir")

cntrl <- list(num_warmup = 2000, num_samples = 4000, num_chains = 4, replicates = 20,
                num_cores = 10, refresh = 2000, adapt_delta = 0.8, max_depth = 10)

fit_model <- function(model, metric, branch, control = cntrl) {
    mod <- cmdstan_model(glue("{model}/{model}.stan"))

    fit <- mod$sample(data = glue("{model}/{model}.json"),
                      num_warmup = control$num_warmup, num_samples = control$num_samples,
                      num_chains = control$num_chains * control$replicates,
                      num_cores = control$num_cores,
                      metric = metric, refresh = control$refresh,
                      adapt_delta = control$adapt_delta, max_depth = control$max_depth)

    summary_diagnostics(fit, control$num_chains, control$replicates, branch)
}

force_recompile <- function(basepath, model_names = models) {
    paths <- sapply(model_names, function(m) glue("{basepath}/{m}/{m}"))
    file.remove(paths)
    invisible(NULL)
}

combine <- function(x, ...) {
    others <- list(...)
    for (nm in names(x)) {
        for (o in others) {
            x[[nm]] <- cbind(x[[nm]], o[[nm]])
        }
    }
    x
}

identify_df <- function(df, rownames, branch) {
    df %>%
        tibble %>%
        mutate(id = rownames, branch = branch) %>%
        pivot_longer(cols = -c(id, branch))
}

standardize_by <- function(ess, leapfrog_run) {
    data.frame(t(apply(ess, 1, function(ess_p) ess_p / leapfrog_run)))
}

summary_diagnostics <- function(fit, nchains, replicates, branch) {
    draws <- fit$draws()
    parameter_names <- fit$sampling_info()$model_params
    idx <- 1:nchains

    d <- foreach(i = 0:(replicates - 1), .combine=combine) %dopar% {
        x <- summarise_draws(draws[, i*nchains + idx, ])
        x2 <- summarise_draws(draws[, i*nchains + idx, ] ^ 2)

        list(essbulk = x$ess_bulk,esstail = x$ess_tail, rhat = x$rhat,
             ess2bulk = x2$ess_bulk, ess2tail = x2$ess_tail)
    }

    leapfrog <- colSums(fit$sampler_diagnostics()[,,"n_leapfrog__"])[,1]
    leapfrog_run <- sapply(0:(replicates - 1), function(i) sum(leapfrog[i*nchains + idx]))

    time <- fit$time()$chains$total
    time_run <- sapply(0:(replicates - 1), function(i) sum(time[i*nchains + idx]))

    list(essbulk_leapfrog = d$essbulk %>%
             standardize_by(leapfrog_run) %>%
             identify_df(parameter_names, branch),
         esstail_leapfrog = d$esstail %>%
             standardize_by(leapfrog_run) %>%
             identify_df(parameter_names, branch),
         ess2bulk_leapfrog = d$ess2bulk %>%
             standardize_by(leapfrog_run) %>%
             identify_df(parameter_names, branch),
         ess2tail_leapfrog = d$ess2tail %>%
             standardize_by(leapfrog_run) %>%
             identify_df(parameter_names, branch),

         essbulk_time = d$essbulk %>%
             standardize_by(time_run) %>%
             identify_df(parameter_names, branch),
         esstail_time = d$esstail %>%
             standardize_by(time_run) %>%
             identify_df(parameter_names, branch),
         ess2bulk_time = d$ess2bulk %>%
             standardize_by(time_run) %>%
             identify_df(parameter_names, branch),
         ess2tail_time = d$ess2tail %>%
             standardize_by(time_run) %>%
             identify_df(parameter_names, branch),

         rhat = d$rhat %>% data.frame %>% identify_df(parameter_names, branch),
         time = identify_df(fit$time()$chains[,"total",drop=FALSE],
                            fit$time()$chains$chain_id,
                            branch),
         leapfrog = identify_df(data.frame(leapfrog), fit$time()$chains$chain_id, branch))
}

plot_ess <- function(df1, df2) {
    bind_rows(df1, df2) %>%
        ggplot(aes(id, value, color=branch, shape=branch)) +
        geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
        stat_summary(fun = "mean", color="black", position = position_dodge(width=0.4)) +
        labs(x = "parameter") +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

plot_time <- function(df1, df2) {
    bind_rows(df1, df2) %>%
        ggplot(aes(value, color=branch)) +
        geom_density() +
        labs(x = "seconds")
}

plot_rhat <- function(df1, df2) {
    bind_rows(df1, df2) %>%
        ggplot(aes(id, value, color = branch, shape = branch)) +
        geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
        geom_hline(yintercept = 1.01) +
        labs(x = "parameter", y = "Rhat") +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

write_model_output <- function(output, model, branch) {
    for(n in names(output)) {
        write_csv(outp[[n]], path = glue("{model}/{branch}_{n}.csv.gz"))
    }
    invisible(NULL)
}

read_model_output <- function(model, branch) {
    files <- c("essbulk_leapfrog", "esstail_leapfrog", "ess2bulk_leapfrog", "ess2tail_leapfrog",
              "essbulk_time", "esstail_time", "ess2bulk_time", "ess2tail_time",
              "rhat", "time")
    out <- vector("list", length(files))
    names(out) <- files
    for (f in files) {
        out[[f]] <- read_csv(glue("{model}/{branch}_{f}.csv.gz"))
    }
    out
}
