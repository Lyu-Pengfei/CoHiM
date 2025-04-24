#' Construct a 4-state transition matrix with given stationary distribution
#'
#' Given a 4-dimensional stationary distribution `pi`, this function constructs 
#' a symmetric-like transition matrix `A` that approximately preserves the stationary distribution.
#' The matrix is designed to have higher self-transition probabilities (diagonal entries)
#' and evenly distributes the remaining mass to the off-diagonal elements.
#'
#' @param J Number of hypotheses (default 10000)
#' @param pi Initial stationary distribution over 4 hidden states
#' @param A Transition matrix of size 4x4
#' @param muA Mean effect size for study A under alternatives
#' @param muB Mean effect size for study B under alternatives
#' @param sdA Standard deviation for study A
#' @param sdB Standard deviation for study B
#'
#' @return A list containing simulated p-values (pa, pb) and true hidden states (theta1, theta2)
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


#' Adjusted density ratio function
#'
#' Computes the ratio between the density of a shifted normal and the standard normal 
#' evaluated at the same quantile level. This is often used in local false discovery rate 
#' or rLIS-type calculations under effect size shift.
#'
#' @param t A numeric vector of (1 - p)-values
#' @param mu The effect size shift (mean under the alternative)
#'
#' @return A numeric vector of density ratios

density_func <- function(t, mu){
  return(dnorm(qnorm(1-t) - mu) / dnorm(qnorm(1-t)))
}

#' Construct a 4-state transition matrix with given stationary distribution
#'
#' @param pi A numeric vector of length 4 representing the stationary distribution.  Must be positive and sum to 1.
#' @return A 4x4 transition probability matrix.

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
#' Simulate three-study p-values under an 8-state HMM
#'
#' This function simulates p-values for three studies using a hidden Markov model
#' with 8 hidden states, where each state corresponds to different combinations of the null/alternative
#' hypotheses across three studies. It generates true latent states for each hypothesis and then
#' simulates z-scores and p-values based on study-specific mean and standard deviations.
#'
#' @param J Number of hypotheses to simulate (default is 10000)
#' @param pi Vector of stationary probabilities over 8 hidden states (length 8, sums to 1)
#' @param A 8x8 transition matrix (default: 0.65 * I + 0.05)
#' @param muA Mean effect size for study A under alternative
#' @param muB Mean effect size for study B under alternative
#' @param muC Mean effect size for study C under alternative
#' @param sdA Standard deviation for z-scores in study A
#' @param sdB Standard deviation for z-scores in study B
#' @param sdC Standard deviation for z-scores in study C
#'
#' @return A list containing:
#' \describe{
#'   \item{pa}{Vector of p-values from study A}
#'   \item{pb}{Vector of p-values from study B}
#'   \item{pc}{Vector of p-values from study C}
#'   \item{theta1}{Indicator vector for non-null status in study A}
#'   \item{theta2}{Indicator vector for non-null status in study B}
#'   \item{theta3}{Indicator vector for non-null status in study C}
#' }
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


#' Construct an 8-state transition matrix with given stationary distribution
#'
#' This function constructs a transition probability matrix `A` of size 8x8
#' that approximately preserves a specified stationary distribution `pi`.
#' It ensures diagonal dominance by assigning a higher self-transition probability `a[i]`
#' to each state and evenly distributes the remaining transition mass to off-diagonal entries.
#' This is typically used for modeling Markov chains with 8 latent states in three-study replicability simulations.
#'
#' @param pi A numeric vector of length 8 representing the stationary distribution (must sum to 1)
#'
#' @return A numeric matrix of size 8x8 representing the transition probabilities
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
