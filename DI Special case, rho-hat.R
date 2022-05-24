# Computing the values of Rho Hat
n <- 100 #fixing the number of trials

# We use this to find the cross-product ratio, rho = p1(1-p2)/p2(1-p1)
# Setting the values for p1, q1 and p2, q2
p1 <- 0.2
p2 <- 0.2
q1 <- 1 - p1
q2 <- 1 - p2

# Setting seed
set.seed(200390872)

# Number of simulations
N <- 10000

# Cross-product ratio
rho <- (p1*q2)/(p2*q1)

# Estimate of tau^2 and tau^2/n
tau2 <- rho*(((p2*q1)+q2)/(p2*((q1)^2)))
tau2/n


for (l in 1:N) {
  T <- rbinom(1, n, p1)#random binomial number generation
  if(T==0){l=l-1}
  else
  {nu <- rnbinom(1,T,p2) + T # pascal number generation
  rhohat <- (nu - T)/(n+1-T)  # rhohat estimate 
  v <- sqrt(n)*(rhohat - rho)
  tausq <- ((nu-T)/(n+1-T))*((((n-T)/(T+1))*((n+1)/(T+1)))+
                               (((nu/T)-1)*((n+1)/(T+1))^2))*((T/(n+1-T))^2)
  gut <- n*tausq*((1-(T/n))^2) # Gut variance estimate
  viet <- n*tausq         # Vietnamese variance estimate
  }
}

# Bias
B <- rhohat-rho
mean(B)

# True variance of rhohat
Tr <- var(rhohat)
Tr

# Mean-squared error values for rhohat
MSE <- B^2 + Tr
MSE

# gut variance
gut

# vietnamese variance
viet

# Variance of rhohat
var(rhohat)

# Mean and variance of statistics V defined above
vhat <- mean(v)
sv2 <- var(v)
tauhatsq <- mean(tausq)


# Actual value of rho
rho

# Mean and bias of rhohat
mean(rhohat)
b <- mean(rhohat)-rho
b

# Mean squared errors and variance of rhohat
mse <- mean((rhohat - rho)^2)
mse
varr <- mean((rhohat - mean(rhohat))^2)
varr




