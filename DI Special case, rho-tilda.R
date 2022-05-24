# Computing the values of Rho Tilda 

n <- 1000 #fixing the number of trials

# We use this to find the cross-product ratio, rho = p1(1-p2)/p2(1-p1)
# Setting the values for p1, q1 and p2, q2
p1 <- 0.8
p2 <- 0.8
q1 <- 1 - p1
q2 <- 1 - p2

# Seeting Seed
set.seed(200390872)

# Number of simulations
N <- 10000

# Cross-product ratio
rho <- (p1*q2)/(p2*q1)
# Estimate of tau^2
tau2 <- rho*(((p2*q1)+q2)/(p2*((q1)^2)))
# Asymptotic variance
A <-tau2/n
A

# Simulation asymptotic variance
# Creating empty vectors for variables in the study
T = rep(0, N)
nu =rep(0, N)
rhohat = rep(0, N)
v = rep(0, N)
tausq = rep(0, N)
s2 = rep(0, N)
for (l in 1:N) {
  T[l] <- rbinom(1, n, p1) #random binomial number generation
  if(T[l]==0){l=l-1}
  else
  {nu[l] <- rnbinom(1,T[l],p2) + T[l]   # pascal number generation
  rhohat[l] <- (nu[l] - T[l])/(n-T[l])   # computing rho-hat values
  v[l] <- sqrt(n)*(rhohat[l] - rho)      # n multiplied by statistic
  tausq[l] <- ((nu[l]-T[l])/(n+1-T[l]))*((((n-T[l])/(T[l]+1))*((n+1)/(T[l]+1)))+
                                           (((nu[l]/T[l])-1)*((n+1)/(T[l]+1))^2))*((T[l]/(n+1-T[l]))^2)
  # Variance
  s2[l] <- tausq[l]/n
  }
}

# Bias
B <- rhohat-rho
mean(B)

# True variance
T <- mean((rhohat - mean(rhohat))^2)
T

# Mean squared error
mse <- B^2 + T
mean(mse)

# Mean values for rhotilda
mean(rhotilda)
