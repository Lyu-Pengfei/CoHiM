SimuData <- function(J     = 10000,
                     pi = c(0.25, 0.25, 0.25, 0.25),
                     A = 0.6 * diag(4) + 0.1,
                     muA   = 2,
                     muB   = 2,
                     sdA   = 1,
                     sdB   = 1){
  s <- c()
  s[1] <- sample(0:3, 1, prob = pi)
  for (j in 2:J){
    s[j] <- sample(0:3, 1, prob = A[s[j-1]+1,])
  }
  
  states1 = rep(0, J)
  states1[c(which(s == 2), which(s == 3))] = 1
  states2 = rep(0, J)
  states2[c(which(s == 1), which(s == 3))] = 1
  
  xa <- rnorm(J, mean = muA * states1, sd = sdA)
  xb <- rnorm(J, mean = muB * states2, sd = sdB)
  
  pa <- 1 - pnorm(xa)
  pb <- 1 - pnorm(xb)
  
  return(list(
    pa = pa,
    pb = pb,
    theta1 = states1,
    theta2 = states2
  ))
}

density_func <- function(t, mu){
  return(dnorm(qnorm(1-t) - mu) / dnorm(qnorm(1-t)))
}


trans.mat <- function(pi){
  a = c()
  a[1] = 1 - min(pi[2:4]/pi[1])*2/3
  A = matrix(NA, 4, 4)
  for (i in 2:4) {
    a[i] = 1 - pi[1]/pi[i] * (1-a[1])
  }
  for (i in 1:4) {
    A[i, ] = (1-a[i])/3
    A[i, i] = a[i]
  }
  return(A)
}



# Three studies -----------------------------------------------------------------
SimuData3 <- function(J     = 10000,
                      pi = rep(1/8, 8),
                      A = 0.65 * diag(8) + 0.05,
                      muA   = 2, muB   = 2, muC   = 2,
                      sdA   = 1, sdB   = 1, sdC   = 1){
  s <- c()
  s[1] <- sample(0:7, 1, prob = pi)
  for (j in 2:J){
    s[j] <- sample(0:7, 1, prob = A[s[j-1]+1,])
  }
  
  states1 = rep(0, J)
  states1[s %/% 4 == 1] = 1
  states2 = rep(0, J)
  states2[s %/% 2 %% 2 == 1] = 1
  states3 = rep(0, J)
  states3[s%%2 == 1] = 1
  
  xa <- rnorm(J, mean = muA * states1, sd = sdA)
  xb <- rnorm(J, mean = muB * states2, sd = sdB)
  xc <- rnorm(J, mean = muC * states3, sd = sdC)
  
  pa <- 1 - pnorm(xa)
  pb <- 1 - pnorm(xb)
  pc <- 1 - pnorm(xc)
  
  return(list(
    pa = pa,
    pb = pb,
    pc = pc,
    theta1 = states1,
    theta2 = states2,
    theta3 = states3
  ))
}

trans.mat3 <- function(pi){
  a = c()
  a[1] = 1 - min(pi[2:8]/pi[1])*2/3
  A = matrix(NA, 8, 8)
  for (i in 2:8) {
    a[i] = 1 - pi[1]/pi[i] * (1-a[1])
  }
  for (i in 1:8) {
    A[i, ] = (1-a[i])/7
    A[i, i] = a[i]
  }
  return(A)
}
