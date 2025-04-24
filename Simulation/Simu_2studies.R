library(tidyverse)
library(patchwork)
library(CoHiM)
source("./Simu_Func.R")



#-------------------------------
# 1. Parameter settings
#-------------------------------
pi1_values <- c(0.1, 0.15)
muA_values <- c(3, 2.5, 2)
muB_values <- c(3, 2.5, 2)
sdA = 1; sdB = 1
J = 10000
q = 0.05
methods <- c("CoHiM")
pi = rep(1/4, 4)
n_reps <- 3

#-------------------------------
# 2. Result container
#-------------------------------
results <- expand.grid(
  pi1   = pi1_values,
  muA   = muA_values,
  muB   = muB_values,
  method= methods,
  stringsAsFactors = FALSE
)
results$FDR_mean   <- NA_real_
results$FDR_sd     <- NA_real_
results$Power_mean <- NA_real_
results$Power_sd   <- NA_real_

#-------------------------------
# 3. Simulation and method evaluation
#-------------------------------
start_time = Sys.time()
for(i in seq_len(nrow(results))) {
  pi1_val    <- results$pi1[i]
  muA_val    <- results$muA[i]
  muB_val    <- results$muB[i]
  method_name<- results$method[i]
  cat(sprintf("pi1 = %.3f, muA = %.1f, muB = %.1f, method = %s\n", pi1_val, muA_val, muB_val, method_name))
  
  pi[4] = 0.1; pi[2:3] = pi1_val; pi[1] = 1 - sum(pi[2:4])
  A = trans.mat(pi)
  
  fdrs   <- numeric(n_reps)
  powers <- numeric(n_reps)
  
  for(rep in seq_len(n_reps)) {
    data.obj <- SimuData(J  = J,
                         pi = pi,
                         A  = A,
                         muA = muA_val, muB = muB_val,
                         sdA = sdA, sdB = sdB)
    
    pa = data.obj$pa
    pb = data.obj$pb
    theta1 = data.obj$theta1
    theta2 = data.obj$theta2
    truth <- theta1 * theta2
    
    if(method_name == "CoHiM"){
      res.hmm <- RepLIS(pa, pb)
      FDR <- sum(res.hmm$fdr <= q & !truth) / max(sum(res.hmm$fdr <= q), 1)
      Power <- sum(res.hmm$fdr <= q & truth) / sum(truth)
    }
    
    fdrs[rep]   <- FDR
    powers[rep] <- Power
  }
  
  results$FDR_mean[i]   <- mean(fdrs)
  results$FDR_sd[i]     <- sd(fdrs)
  results$Power_mean[i] <- mean(powers)
  results$Power_sd[i]   <- sd(powers)
}

#-------------------------------
# 4. Factor for facet layout
#-------------------------------
results$pi1 <- factor(results$pi1, levels = pi1_values)
results$muA <- factor(results$muA, levels = muA_values)
results$muB <- factor(results$muB, levels = muB_values)
results$method <- factor(results$method, levels = methods)

#-------------------------------
# 5. Plotting
#-------------------------------

fdr_plot <- ggplot(results, aes(x = pi1, y = FDR_mean, fill = method, shape = method)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = FDR_mean - FDR_sd, ymax = FDR_mean + FDR_sd),
                position = position_dodge(width = 0.9), width = 0.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 1) +
  xlab(expression(pi[1])) +
  ylab("Empirical FDR") +
  ylim(0, 0.125) +
  facet_grid(muB~muA, labeller = label_parsed, scales = "free_y") +
  theme(legend.position = 'none', 
        text = element_text(size = 28),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 20))

power_plot <- ggplot(results, aes(x = pi1, y = Power_mean, fill = method, shape = method)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = Power_mean - Power_sd, ymax = Power_mean + Power_sd),
                position = position_dodge(width = 0.9), width = 0.2) +
  xlab(expression(pi[1])) +
  ylab("Power") +
  ylim(0, 1) +
  facet_grid(muB~muA, labeller = label_parsed, scales = "free_y") +
  theme(legend.position = 'right',
        text = element_text(size = 28),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 20))

print(fdr_plot + power_plot)
