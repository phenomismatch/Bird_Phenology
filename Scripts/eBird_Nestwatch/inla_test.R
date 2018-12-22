n <- 10
A <- matrix(runif(n^2), n, n)
Q <- A %*% t(A)
print(mean(abs(inla.qinv(Q) - solve(Q))))

## sparse matrix example
rho <- 0.9
Q <- toeplitz(c(1+rho^2, -rho,  rep(0, n-3), -rho)) / (1-rho^2)
Q <- inla.as.dgTMatrix(Q)
Q.inv <- inla.qinv(Q)
