t0 <- Sys.time()

setwd("/home/projects/cu_10184/data/projects/groenbaek/data/data_raw/EVI2_EPICarray/")
library(tidyverse)
library(lme4)
library(vroom)
options(stringsAsFactors = F)

# Extract input arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 6) {
      samplesize <- as.integer(args[1]) # Number of patients in simulated dataset (0 = use original dataset)
      n_cpg      <- as.integer(args[2]) # Number of CpG sites in simulated dataset (0 = use original dataset)
      eff        <- as.integer(args[3]) # Simulated treatment effect to be used
      nsim       <- as.integer(args[4]) # Number of simulated datasets per simulation
      nperm      <- as.integer(args[5]) # Number of permutations per resampling procedure
      res.dist   <- as.integer(args[6]) # Residual distribution: 0: normal, 1: empirical residuals
} else {
      stop("Incorrect number of input arguments")
}

# For interactive testing:
# samplesize <- 0
# n_cpg      <- 720
# eff        <- 1
# nsim       <- 100
# nperm      <- 100
# res.dist   <- 1

runname <- paste0("N", samplesize, "_CpG", n_cpg, "_eff", eff,  "_nsim", nsim, "_nperm", nperm, "_resdist", res.dist)
print(runname)

# Functions to be used:
inverse_logit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

# Different g functions to transform ATE_k when calculating AATE
gs <- list(g1 = function(x) x,
           g2 = function(x) if_else(x > 0, 0, x),
           g3 = function(x) if_else(x > 0, 0, -(x^2)))

calculate_AATE <- function(data, varname, g_function = g){
      visit1_vitamin_aids <- data %>% filter(visit == 1 & treatment_arm == "vitaminC") %>% select(arrayID) %>% distinct() %>% pull(arrayID)
      visit5_vitamin_aids <- data %>% filter(visit == 5 & treatment_arm == "vitaminC") %>% select(arrayID) %>% distinct() %>% pull(arrayID)
      visit1_placebo_aids <- data %>% filter(visit == 1 & treatment_arm == "placebo") %>% select(arrayID) %>% distinct() %>% pull(arrayID)
      visit5_placebo_aids <- data %>% filter(visit == 5 & treatment_arm == "placebo") %>% select(arrayID) %>% distinct() %>% pull(arrayID)
      
      dwide <- data %>% select(all_of(c("probeID", "arrayID", varname))) %>% 
            pivot_wider(values_from = all_of(varname), names_from = probeID) %>% 
            column_to_rownames(var = "arrayID") %>% 
            as.matrix()
      
      Y_1_vitamin <- dwide[visit1_vitamin_aids, ]
      Y_1_placebo <- dwide[visit1_placebo_aids, ]
      Y_2_vitamin <- dwide[visit5_vitamin_aids, ]
      Y_2_placebo <- dwide[visit5_placebo_aids, ]
      
      Z_vitamin <- Y_2_vitamin - Y_1_vitamin
      Z_placebo <- Y_2_placebo - Y_1_placebo
      
      Z_vitamin_means <- apply(Z_vitamin, 2, mean)
      Z_placebo_means <- apply(Z_placebo, 2, mean)
      
      ATE_k <- Z_vitamin_means - Z_placebo_means
      ATE_k <- g_function(ATE_k)
      return(mean(ATE_k))
}

# Preparing datasets for simulation:
load("/home/projects/cu_10184/data/projects/groenbaek/data/data_raw/EVI2_EPICarray/trial_data_720sites.Rdata")
if(n_cpg == 0 & samplesize == 0){
      # If n_cpg and sample size both are 0, use the original dataset with 720 cpg sites
      dsim <- dlong
} 
if(n_cpg != 0 & samplesize == 0){
      # When investigating alternative CpG numbers but keeping the original patient set:
      mut_pts <- unique(dlong$patient_id[dlong$TET2_mut == "mut"])
      vitc_pts <- unique(dlong$patient_id[dlong$treatment_arm == "vitaminC"])
      dsim <- expand.grid(patient_id = rep(unique(dlong$patient_id), each = 2),
                          probeID = factor(1:n_cpg)) %>% 
            mutate(visit = rep(c(1,5), n_cpg*length(unique(dlong$patient_id))),
                   arrayID = paste("array", patient_id, visit, sep = "_"), 
                   TET2_mut = if_else(patient_id %in% mut_pts, "mut", "WT"),
                   treatment_arm = if_else(patient_id %in% vitc_pts, "vitaminC", "placebo"),
                   treatment_constrained = if_else(treatment_arm == "vitaminC" & visit == 5, "vitaminC", "baseline"),
                   patient_visit = paste(patient_id, visit, sep = "_"),
                   probe_mut = paste(probeID, TET2_mut, sep = "_"))
} 
if(n_cpg != 0 & samplesize != 0){
      # Otherwise create artifical dataset with 50% TET2 mut pts and 1:1 randomization
      dsim <- expand.grid(patient_id = rep(1:samplesize, each = 2),
                          probeID = factor(1:n_cpg)) %>% 
            mutate(visit = rep(c(1,5), n_cpg*samplesize),
                   arrayID = paste("array", patient_id, visit, sep = "_"), 
                   TET2_mut = if_else(patient_id %in% 1:(samplesize/2), "WT", "mut"),
                   treatment_arm = if_else(patient_id %in% seq(1, samplesize, 2), "vitaminC", "placebo"),
                   treatment_constrained = if_else(treatment_arm == "vitaminC" & visit == 5, "vitaminC", "baseline"),
                   patient_visit = paste(patient_id, visit, sep = "_"),
                   probe_mut = paste(probeID, TET2_mut, sep = "_"))
}


# Allowing for differential effects of treatment on CpG sites (including no effect)
if(n_cpg != 0){
      # Selecting positively associated probes with lowest P values in MNC dataset
      topmnc <- vroom("mnc_results_commonprobes.csv") %>% 
            arrange(desc(pval)) %>% 
            filter(estimate > 0) %>% 
            slice_tail(n = n_cpg) %>% 
            pull(probe)
      # Finding these probes' coefficients from EWAS with TET2 mutation status in the granulocyte dataset
      # These coefficients are then scaled with a factor 1/1.82, since 1.79 is the mean of the granulocyte coefficients
      # for the 500 top probes from the mnc dataset.
      # Estimate < 0 are set to 0
      site_effect <- vroom("granulocyte_results_commonprobes.csv") %>% 
            filter(probe %in% topmnc) %>% 
            transmute(coef_scaled = (estimate > 0) *estimate / 1.82,
                      probeID = factor(1:n_cpg))
      
      # Altering the treatment_constrained variable s.t. the design matrix will instead contain the scaled coefficients as treatment effects
      dsim <- dsim %>% 
            left_join(site_effect) %>% 
            mutate(treatment_constrained = if_else(treatment_constrained == "baseline", 0, coef_scaled))
}

### Simulation and permutation procedure:
# Model matrices for fixed and random effects:
X <- model.matrix(~ treatment_constrained + TET2_mut*visit, data = dsim)
Z <- model.matrix(~ -1 + patient_visit, data = dsim)
W <- model.matrix(~ -1 + probe_mut, data = dsim)

# Fixed effects: (entry 2 is treatment effect)
model <- lmer(Mval ~ TET2_mut * visit + (1|patient_visit) + (1|probe_mut), data = dlong)
nullcoefs <- summary(model)$coef[,1]

beta <- c(nullcoefs[1], 0, nullcoefs[2:4])
# Sequence of all possible effects that are simulated (only one used in each call to this script)
possible_effects <- 1-exp(seq(from = 0, to = 1, length.out = 10))
# Setting treatment effect to the selected value from input argument
beta[2] <- possible_effects[eff]

# Random effects:
raneffs <- as.data.frame(VarCorr(model))
tau1 <- raneffs[2,5] # SD for patient x visit product factor
tau2 <- raneffs[1,5] # SD for CpG site factor
sigma <- raneffs[3,5] # Residual SD

# Data objects for permutation sampling:
pids <- unique(dsim$patient_id)
N <- length(unique(dsim$patient_id))

# Matrix with all the permutations, stored as the patients assigned to active treatment in each permutation
# N <- length(unique(dsim$patient_id))
n_active <- length(unique(dsim$patient_id[dsim$treatment_arm == "vitaminC"]))

# Matrix with simulated P values
df.perm.pvalues <- data.frame(matrix(nrow = nsim, ncol = 6))
colnames(df.perm.pvalues) <- paste0("g", rep(1:3, each = 2), "_", rep(c("pval", "AATEhat"), 3))

# Vector of residuals to sample from:
if(res.dist == 1){
      resid <- residuals(model)
}

set.seed(555)
for(s in 1:nsim){
      dsim.s <- dsim
      print(paste("effect", eff, "; simulation step", s))
      write.table(paste("simulation step", s), quote = F, row.names = F, col.names = F,
                  file = paste0("simulations/progress_report_", runname, ".txt"))
      # Simulating random effects
      B1 <- rnorm(n = ncol(Z), mean = 0, sd = tau1) # Patient x visit random effects
      B2 <- rnorm(n = ncol(W), mean = 0, sd = tau2) # Site x TET2 random effects
      # Simulating residuals:
      if(res.dist == 0) epsilon <- rnorm(n = nrow(X), mean = 0, sd = sigma) # Normally distributed residuals
      if(res.dist == 1) epsilon <- sample(resid, size = nrow(X), replace = T)
      # Simulating M values:
      Mval <- X %*% beta + Z %*% B1 + W %*% B2 + epsilon 
      # Transforming to beta scale
      dsim.s$simulated_beta <- inverse_logit(Mval)
      
      AATE_hat_g1 <- calculate_AATE(data = dsim.s, "simulated_beta", g_function = gs[[1]])
      AATE_hat_g2 <- calculate_AATE(data = dsim.s, "simulated_beta", g_function = gs[[2]])
      AATE_hat_g3 <- calculate_AATE(data = dsim.s, "simulated_beta", g_function = gs[[3]])
      
      df.AATE_perm <- data.frame(AATE_perm_g1 = rep(NA, nperm),
                                 AATE_perm_g2 = rep(NA, nperm),
                                 AATE_perm_g3 = rep(NA, nperm))
      
      for(b in 1:nperm){
            sample_b <- sample(x = 1:N, size = n_active, replace = F)
            dsim.s$treatment_arm <- if_else(dsim.s$patient_id %in% pids[sample_b], "vitaminC", "placebo")
            df.AATE_perm$AATE_perm_g1[b] <- calculate_AATE(data = dsim.s, "simulated_beta", g_function = gs[[1]])
            df.AATE_perm$AATE_perm_g2[b] <- calculate_AATE(data = dsim.s, "simulated_beta", g_function = gs[[2]])
            df.AATE_perm$AATE_perm_g3[b] <- calculate_AATE(data = dsim.s, "simulated_beta", g_function = gs[[3]])
      }
      df.perm.pvalues$g1_pval[s] <- mean(df.AATE_perm$AATE_perm_g1 < AATE_hat_g1)
      df.perm.pvalues$g2_pval[s] <- mean(df.AATE_perm$AATE_perm_g2 < AATE_hat_g2)
      df.perm.pvalues$g3_pval[s] <- mean(df.AATE_perm$AATE_perm_g3 < AATE_hat_g3)
      
      df.perm.pvalues$g1_AATEhat[s] <- AATE_hat_g1
      df.perm.pvalues$g2_AATEhat[s] <- AATE_hat_g2
      df.perm.pvalues$g3_AATEhat[s] <- AATE_hat_g3
}

t1 <- Sys.time()
time_spent <- difftime(t1, t0, units = "mins")

write.csv(df.perm.pvalues, file = paste0("simulations/perm.results", runname, ".csv"),
          quote = F, row.names = F)
write.csv(time_spent, file = paste0("simulations/perm.comptime_", runname, ".csv"))

























