library(jsonlite)
library(Matrix)

dim <- 64
m <- matrix(c(1, 0.5, 0.5, 1), nrow=2)
l <- lapply(1:dim, function(x) m)
M <- as.matrix(bdiag(l))

data <- list(d = 2 * dim, mu = rep(0, 2 * dim), S = M)
write(toJSON(data, auto_unbox = TRUE), "gauss_c50_128d.json")
