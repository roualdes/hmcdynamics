source("utilities.R")
set_cmdstan_path("~/cmdstan")
theme_set(theme_minimal())


controls <- list(gauss_128d = cntrl, gauss_256d = cntrl, gauss_32d = cntrl, gauss_64d = cntrl,
                 gauss_c25_128d = cntrl, gauss_c25_256d = cntrl, gauss_c25_32d = cntrl, gauss_c25_64d = cntrl,
                 gauss_c50_128d = cntrl, gauss_c50_256d = cntrl, gauss_c50_32d = cntrl, gauss_c50_64d = cntrl,
                 gauss_c75_128d = cntrl, gauss_c75_256d = cntrl, gauss_c75_32d = cntrl, gauss_c75_64d = cntrl
)

model <- "student_t"
metric <-  "dense_e"
replications <- 10
# force_recompile("~/hmcdynamics", model)
mod <- cmdstan_model(glue("{model}.stan"))


for (model in c("gaussian", "student_t")) {
    for (dim in 2^(1:8)) {
        for (rho in seq(0, .75, by = 0.25)) {
            for (nu in seq(3, 12, by = 3)) {
                for (rep in 1:10) {
                    fit <- mod$sample(data = gendata(2^3, 0, 100), num_cores = 4,
                                      metric = metric, refresh = 2000)
                    write_output(fit, model, "proposal", dim, rho, nu, rep)
                }
            }
        }
    }
}



model <- "gauss_c25_32d"
outp <- read_model_output(model, "proposal")
outd <- read_model_output(model, "develop")


p1 <- plot_time(outd$time, outp$time)
p2 <- plot_rhat(outd$rhat, outp$rhat)
p3 <- plot_ess(outd$essbulk_time, outp$essbulk_time) + labs(y = "ESS_Bulk(x)")
p4 <- plot_ess(outd$esstail_time, outp$esstail_time) + labs(y = "ESS_Tail(x)")
p5 <- plot_ess(outd$ess2bulk_time, outp$ess2bulk_time) + labs(y = "ESS_Bulk(x^2)")
p6 <- plot_ess(outd$ess2tail_time, outp$ess2tail_time) + labs(y = "ESS_Tail(x^2)")
grid.arrange(p1, p2, p3, p4, p5, p6, ncol=2)

bind_rows(outd$ess2tail_time, outp$ess2tail_time) %>%
    filter(id != "lp__") %>%
    group_by(branch) %>%
    summarise(m = mean(value))

bind_rows(outd$esstail, outp$esstail) %>%
    filter(id != "lp__") %>%
    group_by(branch) %>%
    summarise(m = mean(value))


p1 <- plot_ess(outd$essbulk_leapfrog, outp$essbulk_leapfrog) + labs(y = "ESS_Bulk(x)")
p2 <- plot_ess(outd$esstail_leapfrog, outp$esstail_leapfrog) + labs(y = "ESS_Tail(x)")
p3 <- plot_ess(outd$ess2bulk_leapfrog, outp$ess2bulk_leapfrog) + labs(y = "ESS_Bulk(x^2)")
p4 <- plot_ess(outd$ess2tail_leapfrog, outp$ess2tail_leapfrog) + labs(y = "ESS_Tail(x^2)")
grid.arrange(p1, p2, p3, p4, ncol=2)
