#fixing the number of trials
n <- 10000 
# We use this to find the cross-product ratio, rho = p1(1-p2)/p2(1-p1)
# Setting the values for p1, q1 and p2, q2
p1 <- 0.8
p2 <- 0.8
q1 <- 1 - p1
q2 <- 1 - p2
# Setting seed
set.seed(200390872)

# Number of simulations
N <- 10000

# Actual value of rho
rho <- (p1*q2)/(p2*q1)
# Actual value of tau2
tau2 <- rho*(((p2*q1)+q2)/(p2*((q1)^2)))

# Creating empty vectors for the following
T = rep(0, N)
nu =rep(0, N)
rhohat = rep(0, N)
v = rep(0, N)
tausq = rep(0, N)
s2 = rep(0, N)
lower = rep(0,N)
upper = rep(0,N)
for (l in 1:N) {
  T[l] <- rbinom(1, n, p1) #random binomial number generation
  if(T[l]==0){l=l-1}
  else
  {nu[l] <- rnbinom(1,T[l],p2) + T[l]   # pascal number generation
  rhohat[l] <- (nu[l] - T[l])/(n+1-T[l])  # computing rho-hat values
  v[l] <- sqrt(n)*(rhohat[l] - rho)       # n multiplied by statistic
  # plus in estimators for tausq 
  tausq[l] <- ((nu[l]-T[l])/(n+1-T[l]))*((((n-T[l])/(T[l]+1))*((n+1)/(T[l]+1)))+
                                           (((nu[l]/T[l])-1)*((n+1)/(T[l]+1))^2))*((T[l]/(n+1-T[l]))^2)
  # Standard error for tausq
  s2[l] <- tausq[l]/n
  # lower limit
  lower[l] <- rhohat[l] - qnorm(0.975)*sqrt(s2[l])
  # upper limit
  upper[l] <- rhohat[l] + qnorm(0.975)*sqrt(s2[l])
  }
}

# Checking the condition to see how many of lower and upper lies in the
# 95% confidence interval range
cond <- numeric(N)
for(l in 1:N){
  if(rho >= lower[l] & rho <= upper[l]) cond[l]=1
}

# Finding the coverage probability by taking the mean.
mean(cond)

# Computing the interval width
width <- upper - lower
mean(width)

# Finding the standard error
se <- sqrt(s2) 
mean(se)

# Computing the variance of rhohat
var(rhohat)