# Simulate datasets for fitting compositional metabarcoding data
library(tidyverse)
library(MCMCpack)
library(compositions)
library(rstan)
library(dplyr)
library(here)
#library(cmdstanr) #optional - only use if optimizing() crashes R session
  #check_cmdstan_toolchain()
  #install_cmdstan()

NAME <- "Test_A_1"

# Specifications for simulation
set.seed(122)
N_site    <- 40 # how many unique sites did you sample (assumes each site is sampled once)
N_species <- 25 # how many species are in your sample
tech_reps  <- 1 # how many technical replicates are you simulating
tech_reps_mock  <- 1 # how many technical replicates of the mock community are you simulating
PCR_cycles <- 35 # this is how many cycles for all of your field samples and mock communities

diri_alpha <- 0.2 # amount of variability in true species proportions at each Site
                  # This is a dirichlet parameter (smaller means more variability, bigger means more even)

# Are you using variable pcr cycle in your simulation?
vary_PCR <- TRUE
N_samp_vary_PCR <- max(round(0.1*N_site),2) #number of sites to choose to vary PCR cycles on
vary_PCR_cycles  <- c(26,29,32) # Number of cycles at which to run variable PCRs 
if(isTRUE(vary_PCR) ==FALSE){ # Catch if you specified vary_PCR = FALSe
  N_samp_vary_PCR <- NULL
  vary_PCR_cycles <- NULL
}

# Define amplification efficiency for each species
alpha_true <- rbeta(N_species, 18,12)

# Define overdispersion for each species
phi0_true  <- rnorm(N_species,10,1)

source(here("code","Simulating_data","quant_metabar_sim.R"))

########################################################################
#### Create data frames that can be read into Stan model
########################################################################

OPTIM = FALSE
CMDSTAN = TRUE
BAYES = TRUE
UNCERT = FALSE
# FALSE if you want only point estimate, (THIS WORKS)
# TRUE if you want approximate estimates of asymptotic errors (THIS HAS HESSIAN that is not pos-def)
# FALSE is faster than TRUE... sometime by a lot.

source(here("code","Simulating_data","Metabarcoding_make_design_files.R"))

# Define the parameters you'd like the model to save.
stan_pars <- c(
  "alpha",
  "beta",
  "eta_samp",
  "eta_mock",
  #"tau_base",
  "tau",
  "mu_samp",
  "mu_mock",
  "int_samp_small"
)

if(OPTIM ==TRUE){
  DRAWS = 1000 # only relevant if UNCERT == TRUE
  source(here("code","Simulating_data","Metabarcoding_optimizing.R")) #output is object point_est
  plot(point_est$alphas$value~alpha_true, main = "Likelihood, Optimizer")
  
  plot(point_est$int_samp_small ~ p_true, main = "Likelihood, Optimizer")
  abline(0,1)
}
if(CMDSTAN ==TRUE){
  source(here("code","Simulating_data","Metabarcoding_optimizing_cmdstanr.R")) #output is object point_est
  plot(stanOpt_summary$estimate[which(substr(stanOpt_summary$variable,1,6)=="alpha[")]~alpha_true, 
       main = "Likelihood, cmdstanr Optimizer", ylab="alpha Estimate", xlab="alpha Known")
  
  plot(stanOpt_summary$estimate[which(substr(stanOpt_summary$variable,1,14)=="int_samp_small")] ~ c(p_true), 
       main = "Likelihood, cmdstanr Optimizer", ylab="p true Estimate", xlab="p true Known")
  abline(0,1)
}
if(BAYES == TRUE){
  # Set up Preferences for MCMC.
  N_CHAIN = 3
  Warm = 50
  Iter = 50
  Treedepth = 11
  Adapt_delta = 0.7
  source(here("code","Simulating_data","Metabarcoding_bayesian.R")) #output is object stanMod
  
  plot(summary(stanMod, par = "alpha")$summary[,1] ~ alpha_true, main = "Bayesian")
  meanEstProp <- summary(stanMod, par = "int_samp_small")$summary[,1] %>% unlist %>% matrix(ncol = N_species, byrow = T)
  plot(meanEstProp[1,]~p_true[1,], main = "Bayesian")
  abline(0,1)
}



