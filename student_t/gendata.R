library(jsonlite)
library(Matrix)

dim <- 1
m <- matrix(c(1, 0.0, 0.0, 1), nrow=2)
l <- lapply(1:dim, function(x) m)
M <- as.matrix(bdiag(l))

data <- list(d = 2 * dim, S = M, nu = 1)
write(toJSON(data, auto_unbox = TRUE), "student_t.json")
