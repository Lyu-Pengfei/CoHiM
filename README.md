# CoHiM

**Co**mposite null hypotheses testing under **Hi**dden **M**arkov models.

## Overview
CoHiM is an empirical Bayesian method for testing compostie null hypotheses under dependence. It captures the depenednce structure by a hidden Markov model. 
CoHiM provides effective control of false discovery rate and higher power for the composite null hypotheses testing by borrowing information of the depenedence strucuture across different studies and different genes. It is scalable to high-diemnsional datasets with tens of thousands of genes measured on tens of thousands of spatial locations. 

## Installation
```R
## Install CoHiM
install.packages("devtools")
devtools::install_github("Lyu-Pengfei/CoHiM")

## Load CoHiM
library(CoHiM)
```

## A numeric example
We illustrate the use of CoHiM for testing composite null hypotheses in a three-study large-scale multiple testing setting using simulated data.

The data-generating functions are given below:
```R
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
```
The data is generated as shown below, and CoHiM is then applied to perform replicability analysis.
```R
## Pre-specify the number of hypotheses, m, the stationary probabilities of the joint hidden states, pi, the transition probabilites of the joint hidden states, A, and the alternative settings
J = 10000; pi = c(0.7, 0.1, 0.1, 0.1); A = trans.mat(pi)
muA = 2.5; muB = 3; sdA = 1; sdB = 1
data.obj <- SimuData(J = J,
                     pi = pi,
                     A  = A,
                     muA = muA, muB = muB,
                     sdA = sdA, sdB = sdB)

pa = data.obj$pa
pb = data.obj$pb
theta1 = data.obj$theta1
theta2 = data.obj$theta2
truth <- theta1 * theta2

## Apply CoHiM
library(CoHiM)
q = 0.05
p_matrix = rbind(pa, pb)
res.cohim = CoHiM(p_matrix, q)
rej_index = res.cohim$rej_index
FDP = sum(rej_index & !truth) / max(sum(rej_index), 1)
power = sum(rej_index & truth) / sum(truth)

## Display the results
cat(paste("There are totally", sum(rej_index), "discoveries.\nFDP =", round(FDP, 3),"\nPower =", round(power, 3)))
# There are totally 607 discoveries.
# FDP   = 0.046 
# Power = 0.572



```

