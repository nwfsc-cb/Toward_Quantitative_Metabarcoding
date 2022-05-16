
# Initial values for some of the parameters.
stan_init_f3 <- function(N_species){#J_seb,K_seb){
  A <- list(
    tau = runif(N_species-1,0.1,0.5),
    alpha_raw = runif(N_species-1,-0.1,0.1)
  )
  return(A)
}

M <- stan_model(file=here("code","Simulating_data","quant_metabar_multinom.stan"))

if(UNCERT == FALSE){
  stanOpt <- optimizing(M, data=stan_data, iter=30000,draws=0,
                      verbose=T, 
                      #tol_param=1e-15, 
                      algorithm="LBFGS",
                      init=stan_init_f3(N_species=N_species))
}

if(UNCERT == TRUE){
  stanOpt <- optimizing(M, data=stan_data,  iter=400000,draws=10,
                        verbose=T,hessian=TRUE,
                        #output_dir=here("Output files","cmdstanr"),
                        tol_obj = 1e-20,
                        tol_rel_obj = 1e-15,
                        tol_grad = 1e-20,
                        tol_param=1e-20,
                        tol_rel_grad=1e-20,
                        algorithm="LBFGS",
                        history_size=5,
                        init=stan_init_f3(N_species=N_species))
}


# Pull out the point estimates in more reasonable from than the original storage
opt_fit <- data.frame(value = stanOpt$par) %>% mutate(par=rownames(.),id=1:nrow(.))
point_est <- list()

#parameters
point_est$alphas <- opt_fit %>% filter(grepl('alpha',par)) %>% filter(!grepl('raw',par)) 
point_est$betas <- opt_fit %>% filter(grepl('beta',par)) %>% filter(!grepl('raw',par)) 

point_est$betas_raw <- opt_fit %>% filter(grepl('beta_raw',par)) 
# make betas into the right size matrix each species with it's own column.
point_est$beta_mat <- matrix(point_est$betas$value,N_b_samp_col,N_species,byrow = TRUE)

# random effects for each observation
point_est$eta_samp <- opt_fit %>% filter(grepl('eta_samp',par)) %>% filter(!grepl('raw',par)) 
point_est$eta_mock <- opt_fit %>% filter(grepl('eta_mock',par)) %>% filter(!grepl('raw',par)) 
# make these into the right size matrix each species with it's own column.
point_est$eta_samp_mat <- matrix(point_est$eta_samp$value,N_obs_samp,N_species,byrow = FALSE)
point_est$eta_mock_mat <- matrix(point_est$eta_mock$value,N_obs_mock,N_species,byrow = FALSE)

# predictions for each observation
point_est$mu_samp <- opt_fit %>% filter(grepl('mu_samp',par)) 
point_est$mu_mock <- opt_fit %>% filter(grepl('mu_mock',par)) 
# make these into the right size matrix each species with it's own column.
point_est$mu_samp_mat <- matrix(point_est$mu_samp$value,N_obs_samp,N_species,byrow = FALSE)
point_est$mu_mock_mat <- matrix(point_est$mu_mock$value,N_obs_mock,N_species,byrow = FALSE)

# varaince parameters for each species
point_est$tau <- opt_fit %>% filter(grepl('tau',par))

# Predicted values for each unique site

A<- opt_fit %>% filter(grepl('int_samp_small',par))
point_est$int_samp_small <- opt_fit %>% filter(grepl('int_samp_small',par));
point_est$int_samp_small <- matrix(point_est$int_samp_small$value,N_obs_samp_small,N_species,byrow = FALSE)

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
  stanOpt = stanOpt, # Full stan Model fitted object
  point_est = point_est
)

save(Output,file=here("code","Simulating_data",paste0(NAME,"_optim.Rdata")))

