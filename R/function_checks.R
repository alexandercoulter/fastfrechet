n = 100  # number of samples - nrow(X) and nrow(Y).
p = 10   # number of covariates - ncol(X).
m = 100  # EQF grid density - ncol(Y).

set.seed(31)
mydata = fastfrechet::generate_zinbinom_qf(n = n,
                                           p = p,
                                           m = m)
X = mydata$X  # (n x p) matrix of covariates
Y = mydata$Y  # (n x m) matrix of EQFs, stored row-wise

# Check inputs:
n = 100 # shows warning statement when n = 0
p = 10 # shows a warning for any value less than 4
m = 150 # shows a warning for n = 0
set.seed(31)
mydata = fastfrechet::generate_zinbinom_qf(n = n,
                                           p = p,
                                           m = m)
X = mydata$X  # (n x p) matrix of covariates
Y = mydata$Y  # (n x m) matrix of EQFs, stored row-wise

# dimensions of return objects match the help file
dim(X)
dim(Y)



# frechetreg_univar2wass check --------------------------------------------

# Using the X and Y outputs from the previous function.

# i.e. Estimate conditional QFs:
Q = fastfrechet::frechetreg_univar2wass(X = X,
                                        Y = Y,
                                        Z = NULL,
                                        lambda = NULL,
                                        lower = 0,    # bounded from below by 0
                                        upper = Inf)

dim(Q)
# Plot these conditional QFs beside the EQFs:
mseq = seq(1 / (2 * m), 1 - 1 / (2 * m), len = m)
plot(x = c(), y = c(), xlim = c(0, 1), ylim = c(0, 60),
     main = 'Fréchet Regression QFs', xlab = 'p', ylab = 'Quantile')
for(i in 1:n) lines(mseq, Q[i, ], col = 'black', lwd = 1)


# Check optional Z input
z = 1
Z = matrix(seq(1,p), nrow = z, ncol = p) # Was able to take the optional Z input and when given a vector it did have the error warning
lambda = c(1:p) # ran with the optional lambda input
lambda = matrix(1:p, nrow = 1, ncol = p) # needs a check for data structure for lambda
lambda = c(1:(p-1))

# Try the X matrix as a character matrix:
# X = (matrix(data = as.character(1:(n*p)), nrow = n, ncol = p)) THIS WILL BREAK IT!

Q = fastfrechet::frechetreg_univar2wass(X = X,
                                        Y = Y,
                                        Z = NULL,
                                        lambda = NULL,
                                        lower = 0,    # bounded from below by 0
                                        upper = Inf)
Q
plot(x = c(), y = c(), xlim = c(0, 1), ylim = c(0, 60),
     main = 'Fréchet Regression QFs', xlab = 'p', ylab = 'Quantile')
for(i in 1:n) lines(mseq, Q[i, ], col = 'black', lwd = 1)
