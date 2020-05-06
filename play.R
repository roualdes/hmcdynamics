library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)
library(jsonlite)
library(cmdstanr)
library(posterior)
library(glue)
library(readr)
library(doParallel)
registerDoParallel(10)
source("utilities.R")
set_cmdstan_path("~/cmdstan")
theme_set(theme_minimal())


controls <- list(ark = cntrl,
                 cauchy = modifyList(cntrl, list(adapt_delta = 0.95)),
                 eight_schools = modifyList(cntrl, list(adapt_delta = 0.95)),
                 garch = cntrl)

model <- "student_t"
metric <-  "dense_e"

## proposal, don't forget to switch branches.  Just double check.
force_recompile("~/hmcdynamics", model)
outp <- fit_model(model, metric, "proposal")
write_model_output(outp, model, "proposal")

## develop, don't forget to switch branches.  Just double check.
force_recompile("~/hmcdynamics", model)
outd <- fit_model(model, metric, "develop")
write_model_output(outd, model, "develop")


p1 <- plot_time(outd$time, outp$time)
p2 <- plot_rhat(outd$rhat, outp$rhat)
p3 <- plot_ess(outd$essbulk_time, outp$essbulk_time) + labs(y = "ESS_Bulk(x)")
p4 <- plot_ess(outd$esstail_time, outp$esstail_time) + labs(y = "ESS_Tail(x)")
p5 <- plot_ess(outd$ess2bulk_time, outp$ess2bulk_time) + labs(y = "ESS_Bulk(x^2)")
p6 <- plot_ess(outd$ess2tail_time, outp$ess2tail_time) + labs(y = "ESS_Tail(x^2)")
grid.arrange(p1, p2, p3, p4, p5, p6, ncol=2)

p1 <- plot_ess(outd$essbulk_leapfrog, outp$essbulk_leapfrog) + labs(y = "ESS_Bulk(x)")
p2 <- plot_ess(outd$esstail_leapfrog, outp$esstail_leapfrog) + labs(y = "ESS_Tail(x)")
p3 <- plot_ess(outd$ess2bulk_leapfrog, outp$ess2bulk_leapfrog) + labs(y = "ESS_Bulk(x^2)")
p4 <- plot_ess(outd$ess2tail_leapfrog, outp$ess2tail_leapfrog) + labs(y = "ESS_Tail(x^2)")
grid.arrange(p1, p2, p3, p4, ncol=2)




for (m in models) {
    force_recompile("~/hmcdynamics", m)
    out <- fit_model(model, metric, "proposal")
    write_model_output(out, m, "proposal")
}



###
model <- "gauss_c50_64d"
metric <- "dense_e"
force_recompile("~/hmcdynamics", model)
mod <- cmdstan_model(glue("{model}/{model}.stan"))

fitd <- mod$sample(data = glue("{model}/{model}.json"),
                  num_warmup = cntrl$num_warmup, num_samples = cntrl$num_samples,
                  num_chains = cntrl$num_chains * 1,
                  num_cores = 2,
                  metric = metric, refresh = cntrl$refresh,
                  adapt_delta = cntrl$adapt_delta, max_depth = cntrl$max_depth)

(smryd <- fitd$summary())



force_recompile("~/hmcdynamics", model)
mod <- cmdstan_model(glue("{model}/{model}.stan"))

fitp <- mod$sample(data = glue("{model}/{model}.json"),
                  num_warmup = cntrl$num_warmup, num_samples = 4000,
                  num_chains = cntrl$num_chains * 1,
                  num_cores = cntrl$num_cores,
                  metric = metric, refresh = cntrl$refresh,
                  adapt_delta = 0.9, max_depth = cntrl$max_depth)

samples <- tibble(data.frame(array(fitp$draws(), dim=c(100000*4, 6), dimnames = list(NULL, fitp$sampling_info()$model_params))))

(smryp <- fitp$summary())
