# Generating matrices

A <- matrix(c(1,2,3,4), nrow=2)
B <- matrix(c(4,5,6,7), nrow=2, byrow=TRUE)
c <- 10*A[,1]
x <- B[2,2]

D <- cbind(A,c)
K <- rbind(A,c)

# Addition/subtraction

A+B
A-B

# Matrix multiplication

A %*% B

K %*% B

B %*% K # incompatible dimensions

A*B # element-wise multiplication

A %x% B # Kronecker

# Additional operations

t(A) # transpose

det(A) # determinant

solve(A) # inverse

# generalized inverse
library(MASS)
G <- ginv(D)
D %*% G %*% D

# Cholesky Decomposition

W <- matrix(c(1,1,1,1,2,2,1,2,5), nrow=3)
L <- chol(W)
t(L) %*% L
