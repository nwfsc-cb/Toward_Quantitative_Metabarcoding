
###################################################
###################################################
###################################################
stan_init_f2 <- function(n.chain,N_species){#J_seb,K_seb){
  A <- list()
  for(i in 1:n.chain){
    A[[i]] <- list(
      tau = runif(N_species-1,0.1,0.5),
      alpha_raw = runif(N_species-1,-0.5,0.5)
    )
  }
  return(A)
}


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanMod = stan(file = here("code","Simulating_data","quant_metabar_multinom.stan") ,data = stan_data,
               verbose = FALSE, chains = N_CHAIN, thin = 1,
               warmup = Warm, iter = Warm + Iter,
               control = list(adapt_init_buffer = 175,
                              max_treedepth=Treedepth,
                              stepsize=0.01,
                              adapt_delta=Adapt_delta,
                              metric="diag_e"),
               pars = stan_pars,
               refresh = 10,
               boost_lib = NULL,
               init = stan_init_f2(n.chain=N_CHAIN,N_species=N_species),
               sample_file = paste0("./tmpB.csv")
)



# get_adaptation_info(stanMod)
pars <- rstan::extract(stanMod, permuted = TRUE)
samp_params <- get_sampler_params(stanMod)
#stanMod_summary <- summary(stanMod)$summary

stanMod_summary <- list()
stanMod_summary[["alpha"]] <- summary(stanMod,pars="alpha")$summary
stanMod_summary[["tau"]] <- summary(stanMod,pars=c("tau"))$summary
stanMod_summary[["beta"]] <- summary(stanMod,pars="beta")$summary
stanMod_summary[["eta_samp"]] <- summary(stanMod,pars="eta_samp")$summary
stanMod_summary[["eta_mock"]] <- summary(stanMod,pars="eta_mock")$summary
stanMod_summary[["mu_samp"]] <- summary(stanMod,pars="mu_samp")$summary
stanMod_summary[["mu_mock"]] <- summary(stanMod,pars="mu_mock")$summary
stanMod_summary[["int_samp_small"]] <- summary(stanMod,pars="int_samp_small")$summary

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
  Warm=Warm,
  Iter=Iter,
  # Fitted Objects
  stanMod = stanMod, # Full stan Model fitted object
  pars = pars, # MCMC output
  samp_params=samp_params, # Sampler informaion
  stanMod_summary = stanMod_summary # posterior summaries.
)

save(Output,file=here("code","Simulating_data",paste0(NAME,"_bayes.Rdata")))

#stanMod_summary[["mu_mock_small"]] <- summary(stanMod,pars="mu_mock_small")$summary
#stanMod_summary[["prop_samp"]] <- summary(stanMod,pars="prop_samp")$summary
# stanMod_summary[["eta_samp"]] <- summary(stanMod,pars="eta_samp")$summary
# stanMod_summary[["eta_mock"]] <- summary(stanMod,pars="eta_mock")$summary
#summary(stanMod,pars="prob_pred")

rstan::traceplot(stanMod,pars=c("lp__","alpha[1]","alpha[2]"),inc_warmup=FALSE)

##################################################
# make predicted-observed plots, proportion space (each point is an observation.)
##################################################
# field samples first
# pivot observations to long form
p_samp_trim <- p_samp_all %>% 
  dplyr::select(site,tech_rep,N_seq,N_pcr_samp,contains("obs_sp")) %>%
  mutate(id=1:nrow(p_samp_all))

p_samp_trim <- pivot_longer(p_samp_trim,cols = contains("obs_sp"),names_to="sp",values_to = "count") %>%
  mutate(prop_samp = count/N_seq)

PROBS <-c(0.025,0.05,0.25,0.5,0.75,0.95,0.975)
pred_samp <- apply(pars$mu_samp,c(2,3),mean) %>% as.data.frame() %>% mutate(id=1:nrow(.))
colnames(pred_samp)[grep("V",colnames(pred_samp))] <- paste0("obs_sp_",1:N_species)
pred_samp_long <- pivot_longer(pred_samp,
                               cols = -id,
                               names_to = "sp",
                               values_to = "mean.pred")

pred_samp_q <- apply(pars$mu_samp,c(2,3),quantile,probs=PROBS)
for(i in  1:length(PROBS)){
  A <- pred_samp_q[i,,] %>% as.data.frame()
  colnames(A)[grep("V",colnames(A))] <- paste0("obs_sp_",1:N_species)
  pred_samp_long <- A  %>% as.data.frame() %>% mutate(id=1:nrow(.)) %>%
    pivot_longer(. , cols = -id, names_to = "sp",values_to = paste0("q.",PROBS[i])) %>%
    left_join(pred_samp_long,.)
}

p_obs_pred <- left_join(p_samp_trim,pred_samp_long)

BREAKS <- c(1e-10,1e-8,1e-6,1e-4,1e-3,1e-2,1e-1,0.25,0.5,0.75,1.0)    
ggplot(p_obs_pred) +
  geom_point(aes(x=prop_samp,y=mean.pred,color=sp),alpha=0.5)+
  geom_errorbar(aes(x=prop_samp,ymin=q.0.05,ymax=q.0.95),alpha=0.05,width=0)+
  scale_x_continuous("Observed",trans="sqrt",breaks=BREAKS) +
  scale_y_continuous("Predicted",trans="sqrt",breaks=BREAKS) +
  geom_abline(intercept=0,slope=1,color="red") +
  theme_bw()

p_samp_trim_obs <- p_samp_trim

###############################################    
####### AMPLIFICATION 
###############################################

plot(stanMod, par = "alpha")    

alpha_summ <- stanMod_summary$alpha %>% as.data.frame() %>% mutate(param = rownames(stanMod_summary$alpha))
alpha_summ <- bind_cols(alpha_summ,alpha_true=1+alpha_true,alr_alpha_true=c(alr(1+alpha_true),0))

THESE <- grep("%",colnames(alpha_summ))
colnames(alpha_summ)[THESE] <- paste0("q.",colnames(alpha_summ)[THESE])
names(alpha_summ) = gsub(pattern = "%", replacement = "", x= names(alpha_summ))

ggplot(alpha_summ) +
  geom_point(aes(x=alr_alpha_true,y=mean),alpha=0.3)+
  geom_errorbar(aes(x=alr_alpha_true,ymin=q.2.5,ymax=q.97.5),alpha=0.2,width=0)+
  geom_abline(intercept=0,slope=1,color="red") +
  theme_bw()

plot(stanMod, par = "alpha")

###############################################    
####### True Proportions
###############################################

p_samp_trim <- p_samp_all %>% mutate(id=as.factor(site)) %>% filter(tech_rep==1)
p_samp_trim$site <- as.numeric(as.character(p_samp_trim$site))

p_true_all <- p_true %>% as.data.frame() 
colnames(p_true_all)[grep("V",colnames(p_true_all))] <- paste0("obs_sp_",1:N_species)
p_true_all <- p_true_all %>% mutate(site=1:nrow(.)) %>%
  pivot_longer(.,cols=contains("obs_sp"),names_to="sp",values_to="true_prop") %>%
  left_join(.,p_samp_trim %>% dplyr::select(id,site) )

PROBS <-c(0.025,0.05,0.25,0.5,0.75,0.95,0.975)
prop_samp <- apply(pars$int_samp_small,c(2,3),mean) %>% as.data.frame() %>% mutate(id=1:nrow(.))
colnames(prop_samp)[grep("V",colnames(prop_samp))] <- paste0("obs_sp_",1:N_species)
prop_samp_long <- pivot_longer(prop_samp,
                               cols = -id,
                               names_to = "sp",
                               values_to = "mean.pred")

prop_samp_q <- apply(pars$int_samp_small,c(2,3),quantile,probs=PROBS)
for(i in  1:length(PROBS)){
  A <- prop_samp_q[i,,] %>% as.data.frame()
  colnames(A)[grep("V",colnames(A))] <- paste0("obs_sp_",1:N_species)
  prop_samp_long <- A %>% as.data.frame() %>% mutate(id=1:nrow(.)) %>%
    pivot_longer(. , cols = -id, names_to = "sp",values_to = paste0("q.",PROBS[i])) %>%
    left_join(prop_samp_long,.)
}

prop_samp_long <- prop_samp_long %>% filter(id %in% p_samp_trim$id) %>% mutate(id = as.factor(id)) %>% 
  left_join(p_true_all,.)


BREAKS <- c(1e-10,1e-8,1e-6,1e-4,1e-3,1e-2,1e-1,0.25,0.5,0.75,1.0)    
ggplot(prop_samp_long) +
  geom_point(aes(x=true_prop,y=mean.pred,color=as.factor(sp)),alpha=0.5)+
  geom_errorbar(aes(x=true_prop,ymin=q.0.05,ymax=q.0.95),alpha=0.1,width=0)+
  scale_x_continuous(breaks=BREAKS,trans="sqrt") +
  scale_y_continuous(breaks=BREAKS,trans="sqrt") +
  geom_abline(intercept=0,slope=1,color="red") +
  theme_classic()
#facet_wrap(~site)

#####Make plot of true vs. raw observed.

raw_v_true <- left_join(p_samp_trim_obs %>% mutate(id=as.factor(id)),
                        p_true_all %>% dplyr::select(-id) %>% mutate(site=as.factor(site)))

ggplot(raw_v_true) +
  geom_point(aes(x=true_prop,y=prop_samp,color=as.factor(sp)),alpha=0.5)+
  #geom_errorbar(aes(x=true_prop,ymin=q.0.05,ymax=q.0.95),alpha=0.1,width=0)+
  scale_x_continuous("True Proportion",breaks=BREAKS,trans="sqrt") +
  scale_y_continuous("Observed proportion", breaks=BREAKS,trans="sqrt") +
  geom_abline(intercept=0,slope=1,color="red") +
  theme_classic()




