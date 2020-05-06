library(jsonlite)
library(Matrix)
library(glue)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)
library(jsonlite)
library(cmdstanr)
library(posterior)
library(readr)

gendata <- function(dim, rho, nu) {
    m <- matrix(c(1, rho, rho, 1), nrow = 2)
    l <- lapply(1:(dim / 2), function(x) m)
    M <- as.matrix(bdiag(l))

    list(d = dim, S = M, nu = nu)
}


force_recompile <- function(basepath, model_names = models) {
    paths <- sapply(model_names, function(m) glue("{basepath}/{m}"))
    file.remove(paths)
    invisible(NULL)
}

write_output <- function(fit, model, branch, dim, rho, nu, rep) {

    fit$sampler_diagnostics() %>%
        as.data.frame %>%
        write_csv(glue("output/{model}_{branch}_diagnostics_{rep}rep_{dim}d_{rho}c_{nu}nu.csv.gz"))

    fit$draws() %>%
        as.data.frame %>%
        write_csv(glue("output/{model}_{branch}_draws_{rep}rep_{dim}d_{rho}c_{nu}nu.csv.gz"))

    fit$time()$chains %>%
                 as.data.frame %>%
                 write_csv(glue("output/{model}_{branch}_time_{rep}rep_{dim}d_{rho}c_{nu}nu.csv.gz"))

    fit$summary() %>%
        as.data.frame %>%
        write_csv(glue("output/{model}_{branch}_summary_{rep}rep_{dim}d_{rho}c_{nu}nu.csv.gz"))

    invisible(NULL)
}



## TODO look at everything below again, probably don't need much

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
    diagnostics <- fit$sampler_diagnostics()

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

    list(
        diagnostics = diagnostics %>%
            array(dim = c(prod(dim(diagnostics)[1:2]), dim(diagnostics)[3]),
                  dimnames = list(NULL, fit$sampling_info()$sampler_diagnostics)) %>%
            data.frame %>%
            tibble,

        draws = draws %>%
             array(dim = c(prod(dim(draws)[1:2]), length(parameter_names)),
                   dimnames = list(NULL, parameter_names)) %>%
             data.frame %>%
             tibble,

         essbulk = d$essbulk %>%
             data.frame %>%
             identify_df(parameter_names, branch),
         esstail = d$esstail %>%
             data.frame %>%
             identify_df(parameter_names, branch),
         ess2bulk = d$ess2bulk %>%
             data.frame %>%
             identify_df(parameter_names, branch),
         ess2tail = d$ess2tail %>%
             data.frame %>%
             identify_df(parameter_names, branch),

         essbulk_leapfrog = d$essbulk %>%
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
         leapfrog = identify_df(data.frame(leapfrog), fit$time()$chains$chain_id, branch),
         seed = tibble(seed = fit$sampling_info()$seed))
}

plot_ess <- function(df1, df2, add_overallmeans = FALSE) {

    df <- bind_rows(df1, df2)

    p <- df %>%
        ggplot(aes(id, value, color=branch, shape=branch)) +
        geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.25) +
        stat_summary(fun = "mean", color="black", position = position_dodge(width=0.4), alpha = 0.5) +
        labs(x = "parameter") +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))

    if (add_overallmeans) {
        means <- df %>%
            filter(id != "lp__") %>%
            group_by(branch) %>%
            summarise(m = mean(value))
        p <- p + geom_hline(data=means, aes(yintercept = m, linetype = branch))
    }

    p
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
               "essbulk", "esstail", "ess2bulk", "ess2tail",
               "rhat", "time")
    out <- vector("list", length(files))
    names(out) <- files
    for (f in files) {
        out[[f]] <- read_csv(glue("{model}/{branch}_{f}.csv.gz"))
    }
    out
}
