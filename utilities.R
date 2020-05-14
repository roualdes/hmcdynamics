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
library(stringr)
library(tidyr)

gendata <- function(model, dim, rho) {
    if (model %in% c("gaussian", "student_t")) {
        m <- matrix(c(1, rho, rho, 1), nrow = 2)
        l <- lapply(1:(dim / 2), function(x) m)
        M <- as.matrix(bdiag(l))

        return(list(d = dim, S = M, nu = 3))


    } else if (model == "german_credit") {
        df <- read_csv("german_credit.csv")

        N <- nrow(df)
        M <- ncol(df)

        x <- as.matrix(df[,-M])

        for (m in 1:(M - 1)) {
            x[,m] <- 2 * (x[,m] / (diff(range(x[,m])))) - 1
        }

        return(list(y = df %>% pull(X25) - 1,
                    x = cbind(1, x),
                    n = N,
                    d = M))
    }
}


force_recompile <- function(basepath, model_names = models) {
    paths <- sapply(model_names, function(m) glue("{basepath}/{m}"))
    file.remove(paths)
    invisible(NULL)
}

genpath <- function(model, branch, dim, rho, rep, thing = NULL) {
    dim <- stringr::str_pad(dim, width = 3, pad = 0)
    rho <- stringr::str_pad(rho * 100, width = 2, pad = 0)
    rep <- stringr::str_pad(rep, width = 3, pad = 0)

    path <- glue("output/{model}_{branch}_")
    if (!is.null(thing)) {
        path <- glue(path, "{thing}_")
    }
    glue(path, "{rep}rep_{dim}d_{rho}c.csv.gz")
}

use_seed <- function(model, branch, dim, rho, rep) {
    seed <- NULL
    b <- setdiff(c("proposal", "develop"), branch)

    path <- genpath(model, b, dim, rho, rep, "seed")
    exists <- file.exists(path)
    if (exists) {
        seed <- read_csv(path)$seed
    }

    seed
}

write_output <- function(fit, model, branch, dim, rho, rep) {
    fit$sampler_diagnostics() %>%
        as.data.frame %>%
        write_csv(genpath(model, branch, dim, rho, rep, "diagnostics"))

    fit$time()$chains %>%
                 as.data.frame %>%
                 write_csv(genpath(model, branch, dim, rho, rep, "time"))

    fit$summary() %>%
        as.data.frame %>%
        filter(str_detect(variable, "x")) %>%
        write_csv(genpath(model, branch, dim, rho, rep, "summary"))

    fit$draws()^2 %>%
        summarise_draws %>%
        as.data.frame %>%
        filter(str_detect(variable, "x")) %>%
        write_csv(genpath(model, branch, dim, rho, rep, "summary2"))

    data.frame(seed = fit$sampling_info()$seed) %>%
        write_csv(genpath(model, branch, dim, rho, rep, "seed"))

    invisible(NULL)
}

which_closest_mean <- function(x) {
    m <- mean(x)
    which.min(abs(x - m))
}

prepare_df <- function(model, rho, branches, dims, reps) {
    df <- expand.grid(branch = branches, dim = dims, rep = reps)
    df$ess_bulk_sec <- NA
    df$ess_tail_sec <- NA
    df$ess2_bulk_sec <- NA
    df$ess2_tail_sec <- NA
    df$ess_bulk_leapfrog <- NA
    df$ess_tail_leapfrog <- NA
    df$ess2_bulk_leapfrog <- NA
    df$ess2_tail_leapfrog <- NA
    df$time <- NA
    df$seed <- NA
    df$rhat <- NA
    df$q5 <- NA
    df$q95 <- NA
    df$median <- NA
    df$sd <- NA
    df$leapfrog <- NA
    df$stepsize <- NA
    df$divergent <- NA

    for (b in branches) {
        for (d in dims) {
            for (r in reps) {
                dfs <- suppressMessages(read_csv(genpath(model, b, d, rho, r, "summary")))
                dfs2 <- suppressMessages(read_csv(genpath(model, b, d, rho, r, "summary2")))

                dfseed <- suppressMessages(read_csv(genpath(model, b, d, rho, r, "seed")))
                seed <- dfseed$seed

                dft <- suppressMessages(read_csv(genpath(model, b, d, rho, r, "time")))
                time <- sum(dft$sampling)

                dfd <- suppressMessages(read_csv(genpath(model, b, d, rho, r, "diagnostics")))
                leapfrog <- dfd %>% select(contains("leapfrog")) %>% sum
                stepsize <- dfd %>% select(contains("stepsize")) %>% colMeans %>% mean
                divergent <- dfd %>% select(contains("divergent")) %>% sum

                idx <- with(df, which(branch == b & dim == d & rep == r))

                ## ess per time
                essbulks <- dfs$ess_bulk / time
                esstails <- dfs$ess_tail / time
                essbulk2s <- dfs2$ess_bulk / time
                esstail2s <- dfs2$ess_tail / time

                essbulks_idx <- which.min(essbulks)
                esstails_idx <- which.min(esstails)
                essbulk2s_idx <- which.min(essbulk2s)
                esstail2s_idx <- which.min(esstail2s)

                df$ess_bulk_sec[idx] <- essbulks[essbulks_idx]
                df$ess_tail_sec[idx] <- esstails[esstails_idx]
                df$ess2_bulk_sec[idx] <- essbulk2s[essbulk2s_idx]
                df$ess2_tail_sec[idx] <- esstail2s[esstail2s_idx]

                ## ess per leapfrog
                essbulkl <- dfs$ess_bulk / leapfrog
                esstaill <- dfs$ess_tail / leapfrog
                essbulk2l <- dfs2$ess_bulk / leapfrog
                esstail2l <- dfs2$ess_tail / leapfrog

                essbulkl_ldx <- which.min(essbulkl)
                esstaill_ldx <- which.min(esstaill)
                essbulk2l_ldx <- which.min(essbulk2l)
                esstail2l_ldx <- which.min(esstail2l)

                df$ess_bulk_leapfrog[idx] <- essbulkl[essbulkl_ldx]
                df$ess_tail_leapfrog[idx] <- esstaill[esstaill_ldx]
                df$ess2_bulk_leapfrog[idx] <- essbulk2l[essbulk2l_ldx]
                df$ess2_tail_leapfrog[idx] <- esstail2l[esstail2l_ldx]

                ## other
                df$rhat[idx] <- max(dfs$rhat)
                df$time[idx] <- time
                df$seed[idx] <- seed
                df$leapfrog[idx] <- leapfrog
                df$stepsize[idx] <- stepsize
                df$divergent[idx] <- divergent

                ## stats
                df$q5[idx] <- dfs$q5[esstails_idx]
                df$median[idx] <- dfs$median[essbulks_idx]
                df$q95[idx] <- dfs$q95[esstails_idx]
                df$sd[idx] <- dfs$sd[essbulk2s_idx]
            }
        }
    }
    df
}

plot_ess_time <- function(df, model, rho) {
    df %>%
        select(branch, dim, ess_bulk_sec, ess_tail_sec, ess2_bulk_sec, ess2_tail_sec) %>%
        pivot_longer(cols = c(ess_bulk_sec, ess_tail_sec, ess2_bulk_sec, ess2_tail_sec)) %>%
        ggplot(aes(factor(dim), value, color = branch, shape = branch)) +
        geom_point(position = position_jitterdodge(), alpha = 0.25) +
        scale_y_log10() +
        stat_summary(fun.data = "median_hilow", position = position_jitterdodge()) +
        facet_wrap(~name, nrow=2, scales="free_y") +
        labs(title = bquote("Model" ~ .(model) ~ "with rho =" ~ .(rho)))
}

table_ess_time <- function(df) {
    df %>%
        group_by(branch, dim) %>%
        summarise(mean_ess_bulk = mean(ess_bulk_sec),
                  mean_ess_tail = mean(ess_tail_sec),
                  mean_ess2_bulk = mean(ess2_bulk_sec),
                  mean_ess2_tail = mean(ess2_tail_sec)) %>%
        arrange(dim)
}


plot_ess_leapfrog <- function(df, model, rho) {
    df %>%
        select(branch, dim, ess_bulk_leapfrog, ess2_bulk_leapfrog, ess_tail_leapfrog, ess2_tail_leapfrog) %>%
        pivot_longer(cols = c(ess_bulk_leapfrog, ess2_bulk_leapfrog, ess_tail_leapfrog, ess2_tail_leapfrog)) %>%
        ggplot(aes(factor(dim), value, color = branch, shape = branch)) +
        geom_point(position = position_jitterdodge(), alpha = 0.25) +
        stat_summary(fun.data = "median_hilow", position = position_jitterdodge()) +
        facet_wrap(~name, nrow=2, scales="free_y") +
        labs(title = bquote("Model" ~ .(model) ~ "with rho =" ~ .(rho)))
}

table_ess_leapfrog <- function(df) {
    df %>%
        group_by(branch, dim) %>%
        summarise(mean_ess_bulk = mean(ess_bulk_leapfrog),
                  mean_ess_tail = mean(ess_tail_leapfrog),
                  mean_ess2_bulk = mean(ess2_bulk_leapfrog),
                  mean_ess2_tail = mean(ess2_tail_leapfrog)) %>%
        arrange(dim)
}


plot_stats <- function(df, model, rho) {
    df %>%
        select(branch, dim, q5, median, q95, sd) %>%
        pivot_longer(cols = c(q5, median, q95, sd)) %>%
        ggplot(aes(factor(dim), value, color=branch)) +
        geom_point(position = position_jitterdodge(), alpha = 0.25) +
        stat_summary(fun.data = "median_hilow", position = position_jitterdodge()) +
        facet_wrap(~name, nrow=2, scale="free_y") +
        labs(title = bquote("Model" ~ .(model) ~ "with rho =" ~ .(rho)))
}


plot_rhat <- function(df, model, rho) {
    df %>%
        ggplot(aes(factor(dim), rhat, color = branch, shape = branch)) +
        geom_point(position = position_jitterdodge(), alpha = 0.25) +
        stat_summary(fun.data = "median_hilow", position = position_jitterdodge()) +
        geom_hline(yintercept = 1.01, linetype = 2) +
        labs(x = "Dimension", y = "Rhat", title = bquote("Model" ~ .(model) ~ "with rho =" ~ .(rho)))
}


check_stepsize <- function(df) {
    df %>%
        group_by(branch, dim, rep) %>%
        summarise(m_eps = mean(stepsize)) %>%
        arrange(dim, rep)
}


check_leapfrog <- function(df) {
    df %>%
        group_by(branch, dim, rep) %>%
        summarise(m_leapfrog = mean(leapfrog)) %>%
        arrange(rep, dim)
}


check_runtime <- function(df) {
    df %>%
        group_by(branch, dim) %>%
        summarise(m_time = mean(time)) %>%
        arrange(dim)
}

plot_runtime <- function(df) {
    df %>%
        ggplot(aes(factor(dim), time, color=branch)) +
        geom_point(position = position_jitterdodge(), alpha = 0.25) +
        stat_summary(fun.data = "median_hilow", position = position_jitterdodge())
}


check_divergence <- function(df) {
    df %>%
        group_by(branch, dim) %>%
        summarise(m_divergent = mean(divergent)) %>%
        arrange(dim)
}


check_seed <- function(df) {
    df %>%
        group_by(branch, dim, rep) %>%
        summarise(m_seed = mean(seed)) %>%
        arrange(rep, dim)
}
