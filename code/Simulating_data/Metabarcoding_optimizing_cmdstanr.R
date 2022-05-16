
# Initial values for some of the parameters.
stan_init_f1 <- function(N_species){#J_seb,K_seb){
  A <- list()
  A[[1]] <- list(
    tau = runif(N_species-1,0.1,0.5),
    alpha_raw = runif(N_species-1,-0.1,0.1)
  )
  return(A)
}

M <- cmdstan_model(here("code","Simulating_data","quant_metabar_multinom.stan"))

stanOpt <- M$optimize(data=stan_data, iter=30000,
                      algorithm="lbfgs",
                      init=stan_init_f1(N_species=N_species))

stanOpt_summary <- stanOpt$summary() %>% as.data.frame()

# Pull out the point estimates in more reasonable from than the original storage

####

# Write fit and raw data to file

Output <- list(
  #Simulation parameters
  N_species= N_species,
  N_site = N_site,
  alpha_true = alpha_true,
  PCR_cycles = PCR_cycles,
  tech_reps  = tech_reps,  
  tech_reps_mock  = tech_reps_mock,
  phi0_true  = phi0_true,
  vary_PCR   = vary_PCR,
  N_samp_vary_PCR = N_samp_vary_PCR, 
  vary_PCR_cycles = vary_PCR_cycles, 
  
  # Realizations of simulations (input data)
  p_true = p_true,
  p_samp_all = p_samp_all,
  p_mock_all = p_mock_all,
  
  # stan input objects
  stan_data = stan_data,
  
  # Fitted Objects
  #stanOpt = stanOpt, # Full stan Model fitted object
  stanOpt_summary = stanOpt_summary
)

save(Output,file=here("code","Simulating_data",paste0(NAME,"_optim_cmdstanr.Rdata")))

