# This File pulls in data from Hanfling et al. 2016 and used the even numbered communities as 
# mock communities.

library(tidyverse)
library(compositions)
  select <- dplyr::select  #eliminate conflict w MASS::select vs. dplyr::select
library(rstan)
library(here)

  all <- read.csv(here("data/Hanfling2016/Mock_cytB.csv"), stringsAsFactors = F) %>%
    select(community, Species, nReads, start_conc_ng, Cycles, tech_rep) %>% 
    group_by(community) %>% 
    mutate(propReads = nReads/sum(nReads),
           totReads = sum(nReads)) %>% 
    filter(totReads > 10000) %>% 
    group_by(Species) %>% 
    mutate(totalSpeciesReads = sum(nReads)) %>% 
    filter(totalSpeciesReads > 0) %>% 
    ungroup() %>% 
    mutate(comm_idx = match(community, unique(community)))
  
  A<-all %>% distinct(community,Species) %>% arrange(community,Species) %>% mutate(val=1)
  wideform <- pivot_wider(A,id_cols="Species", names_from="community",values_from="val") %>% as.data.frame()
  
  #to cross-validate, treat some mock comm samples as unknown
  env <- read.csv(here("data/Hanfling2016/Mock_cytB.csv"), stringsAsFactors = F) %>%  #enviro samples
            filter(community %in% c("MC01", "MC03", "MC05", "MC07", "MC09")) %>%
            select(community, Species, nReads, start_conc_ng, Cycles, tech_rep) %>% 
            group_by(community) %>% 
            mutate(propReads = nReads/sum(nReads),
                 totReads = sum(nReads)) %>% 
            filter(totReads > 10000) %>% 
            group_by(Species) %>% 
            mutate(totalSpeciesReads = sum(nReads)) %>% 
            # filter(totalSpeciesReads > 0) %>% 
            ungroup() %>% 
            mutate(comm_idx = match(community, unique(community)))
  
  #assign most common species to be the reference species
  mostCommon <- env %>% 
    group_by(Species) %>% 
    tally(nReads > 0) %>%
    arrange(desc(n)) %>% 
    head(1) %>% 
    pull(Species)
  env$Species[env$Species == mostCommon] <- paste0("zRefSpecies_", mostCommon) #make the most common species the reference species
  env <- env %>% 
    arrange(Species, community)
  
  mc <- read.csv(here("data/Hanfling2016/Mock_cytB.csv"), stringsAsFactors = F) %>%  #mock comm samples
    filter(community %in% c("MC02", "MC04", "MC06", "MC08", "MC10")) %>%
    #filter(nReads > 0) %>% 
    select(community, Species, nReads, start_conc_ng, Cycles, tech_rep) %>% 
    mutate(comm_idx = match(community, unique(community)))
  # mc <- mc %>%
  #     filter(Species %in% env$Species)  #limit to species occurring in environmental dataset

  mc$Species[mc$Species == mostCommon] <- paste0("zRefSpecies_", mostCommon) #make the most common species the reference species
  
    # Filter so that you only keep species in the environment samples that are in the mock community.
    # It is ok to include species that are only in the mock community.
  env <- env %>%
      filter(Species %in% mc$Species)%>% #limit to species occurring in mock community dataset
      arrange(Species, community)  
    
    #double check
    sum(!mc$Species %in% unique(env$Species)) # This can be non-zero
    sum(!env$Species %in% unique(mc$Species)) # this had better be zero.
    
    # Make a single species list:
    sp.list   <- data.frame(Species = sort(unique(mc$Species)) ) %>% mutate(sp_idx =1:length(Species))
    N_species <- nrow(sp.list)
    
    comm.mock.list <- mc %>% group_by(community, tech_rep,Cycles) %>% summarise(n=length(tech_rep)) %>%
      ungroup() %>% mutate(id=1:length(n))
    comm.env.list   <- env %>% group_by(community, tech_rep,Cycles) %>% summarise(n=length(tech_rep)) %>%
      ungroup() %>% mutate(id=1:length(n))
    
    #make a list of species that are in mock community but not environment, 
    # expand grid to make it so the the environmental samples get padded with all the
    # missing species for all of the communities and all tech replicates.
    
    sp.comm.mc  <- expand_grid(Species = sp.list$Species, id = comm.mock.list$id) %>% 
      left_join(.,sp.list %>% select(Species,sp_idx)) %>%
      left_join(.,comm.mock.list %>% select(community,tech_rep,Cycles,id) ) %>% select(-id)
    sp.comm.env <- expand_grid(Species = sp.list$Species, id = comm.env.list$id) %>% 
      left_join(.,sp.list %>% select(Species,sp_idx)) %>%
      left_join(.,comm.env.list %>% select(community,tech_rep,Cycles,id) ) %>% select(-id)
    #convert to matrices
    # merge in species and indices first to make pivoting more efficient.
    
    mc  <- left_join(sp.comm.mc,mc) %>%
      mutate(nReads = ifelse(is.na(nReads),0,nReads),
             start_conc_ng = ifelse(is.na(start_conc_ng),0,start_conc_ng)) 
    env <- left_join(sp.comm.env,env) %>%
      mutate(nReads = ifelse(is.na(nReads),0,nReads),
             start_conc_ng = ifelse(is.na(start_conc_ng),0,start_conc_ng)) 
    
    sample_data <- env %>% 
      ungroup() %>% 
      dplyr::select(community, sp_idx, nReads,tech_rep, Cycles) %>% 
      arrange(sp_idx) %>% 
      pivot_wider(names_from = "sp_idx", values_from = "nReads", values_fill = 0) 
    
    sample_data_small <- sample_data %>% filter(tech_rep==1)
    
    mock_data <- mc %>% 
      ungroup() %>% 
      dplyr::select(community, sp_idx, nReads,tech_rep, Cycles) %>% 
      arrange(sp_idx) %>% 
      pivot_wider(names_from = "sp_idx", values_from = "nReads", values_fill = 0)
    
    mock_data_small <- mock_data %>% filter(tech_rep==1)
    
    #proportions
    p_mock <- mc %>% 
      select(community, tech_rep, sp_idx, start_conc_ng, Cycles) %>% 
      arrange(sp_idx) %>% 
      group_by(community, tech_rep, Cycles) %>% 
      mutate(prop_conc = start_conc_ng/sum(start_conc_ng)) %>% 
      select(-start_conc_ng) %>% #, -Species) %>% 
      pivot_wider(names_from = "sp_idx", values_from = "prop_conc", values_fill = 0) %>% 
      ungroup() %>% 
      arrange(community)
    #select(-community)
    
    p_mock_small <- p_mock %>% 
      filter(tech_rep == 1) 
    
    #calculate additive log ratios 
    alr_mock_true_prop <- p_mock[,4:(ncol(p_mock)-1)]*0
    for(i in 1:nrow(p_mock)){
      alr_mock_true_prop[i,] <- alr(p_mock[i,4:(ncol(p_mock))] + 1e-10)
    }
    alr_mock_true_prop[,N_species] <- 0 #adding explicit reference species column
    
    alr_mock_true_prop_small <- p_mock_small[,4:(ncol(p_mock_small)-1)]*0
    for(i in 1:nrow(p_mock_small)){
      alr_mock_true_prop_small[i,] <- alr(p_mock_small[i,4:(ncol(p_mock_small))] + 1e-10)
    }
    alr_mock_true_prop_small[,N_species] <- 0 
    
    #DESIGN MATRICES
    # mock communities first
    # species compositions (betas)
    # use mock_data  
    
    N_pcr_mock <- mock_data$Cycles
    
    if(length(unique(mock_data$community))==1){
      formula_b <- Cycles ~ 1  # what is on the left side of the equation doesn't matter.
    } else {
      formula_b <- Cycles ~ community # what is on the left side of the equation doesn't matter.
    }
    model_frame <- model.frame(formula_b, mock_data)
    model_matrix_b_mock <- model.matrix(formula_b, model_frame)
    
    #formula_b <- obs_sp_1 ~ community 
    model_frame <- model.frame(formula_b, mock_data_small)
    model_matrix_b_mock_small <- model.matrix(formula_b, model_frame)
    
    # efficiencies (alphas)
    formula_a <- community ~ Cycles -1
    model_frame <- model.frame(formula_a, mock_data)
    model_vector_a_mock <- model.matrix(formula_a, model_frame) %>% as.numeric()
    model_frame <- model.frame(formula_a, mock_data_small)
    model_vector_a_mock_small <- model.matrix(formula_a, model_frame) %>% as.numeric()
    
    N_obs_mock_small <- nrow(model_matrix_b_mock_small)
    N_obs_mock       <- nrow(mock_data)
    N_b_mock_col     <- ncol(model_matrix_b_mock)  
    
    # unknown communities second
    # species compositions (betas)
    
    # use sample_data
    
    N_pcr_samp <- sample_data$Cycles
    
    if(length(unique(sample_data$community))==1){
      formula_b <- Cycles ~ 1  
    } else {
      formula_b <- Cycles ~ community
    }
    model_frame <- model.frame(formula_b, sample_data)
    model_matrix_b_samp <- model.matrix(formula_b, model_frame)
    
    #formula_b <- obs_sp_1 ~ community 
    #p_samp_all$site <- as.factor(p_samp_all$site)
    model_frame <- model.frame(formula_b, sample_data_small)
    model_matrix_b_samp_small <- model.matrix(formula_b, model_frame)
    
    # efficiencies (alpha)
    formula_a <- community ~ Cycles -1
    model_frame <- model.frame(formula_a, sample_data)
    model_vector_a_samp <- model.matrix(formula_a, model_frame) %>% as.numeric()
    model_frame <- model.frame(formula_a, sample_data_small)
    model_vector_a_samp_small <- model.matrix(formula_a, model_frame) %>% as.numeric()
    
    #counters 
    N_obs_samp_small <- nrow(model_matrix_b_samp_small)
    N_obs_samp <- nrow(sample_data)
    N_b_samp_col <- ncol(model_matrix_b_samp)  
    

stan_data <- list(
  N_species = N_species,   # Number of species in data
  N_obs_samp = N_obs_samp, # Number of observed samples 
  N_obs_mock = N_obs_mock, # Number of observed mock samples
  N_obs_samp_small = N_obs_samp_small, # Number of observed samples 
  N_obs_mock_small = N_obs_mock_small, # Number of observed mock samples
  
  # Observed data of community matrices
  sample_data = sample_data[,4:ncol(sample_data)],
  mock_data   = mock_data[,4:ncol(mock_data)],
  
  # True proportions for mock community
  #mock_true_prop = p_mock_all %>% dplyr::select(contains("sp")),
  alr_mock_true_prop = alr_mock_true_prop,
  alr_mock_true_prop_small = alr_mock_true_prop_small,
  
  # vectors of PCR numbers
  N_pcr_samp = as.array(N_pcr_samp),
  N_pcr_mock = N_pcr_mock,
  
  # Design matrices: field samples
  N_b_samp_col = N_b_samp_col,
  model_matrix_b_samp = model_matrix_b_samp,
  model_matrix_b_samp_small = model_matrix_b_samp_small,
  model_vector_a_samp = as.array(model_vector_a_samp),
  model_vector_a_samp_small = as.array(model_vector_a_samp_small),
  
  # Design matrices: mock community samples
  N_b_mock_col = N_b_mock_col,
  model_matrix_b_mock = model_matrix_b_mock,
  model_matrix_b_mock_small = model_matrix_b_mock_small,
  model_vector_a_mock = model_vector_a_mock,
  model_vector_a_mock_small = model_vector_a_mock_small,
  
  # Priors
  alpha_prior = c(0,0.1),  # normal prior
  beta_prior = c(0,10),    # normal prior
  tau_prior = c(1.5,1.5)   # gamma prior
)


stan_pars <- c(
  "alpha",
  "beta",
  "mu_samp",
  "mu_mock",
  "int_samp_small"
)

stan_init_f2 <- function(n.chain,N_species){
  A <- list()
  for(i in 1:n.chain){
    A[[i]] <- list(
      # tau = runif(N_species-1,0.1,0.5),
      alpha_raw = runif(N_species-1,-0.5,0.5)
    )
  }
  return(A)
}

###################################################
###################################################
###################################################
N_CHAIN = 5
Warm = 500
Iter = 1000
Treedepth = 12
Adapt_delta = 0.7

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanMod = stan(file = "./code/Core_Stan_models/quant_metabar_no_overdispersion.stan" ,data = stan_data,
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
               init = stan_init_f2(n.chain=N_CHAIN,N_species=N_species)
               #sample_file = paste0("./tmpB.csv")
)

pars <- rstan::extract(stanMod, permuted = TRUE)
samp_params <- get_sampler_params(stanMod)
#stanMod_summary <- summary(stanMod)$summary

stanMod_summary <- list()
stanMod_summary[["alpha"]] <- summary(stanMod,pars="alpha")$summary
stanMod_summary[["beta"]] <- summary(stanMod,pars="beta")$summary
stanMod_summary[["mu_samp"]] <- summary(stanMod,pars="mu_samp")$summary
stanMod_summary[["mu_mock"]] <- summary(stanMod,pars="mu_mock")$summary
stanMod_summary[["int_samp_small"]] <- summary(stanMod,pars="int_samp_small")$summary


Output <- list(
  #Simulation parameters
  #N_species= N_species,
  #N_site = N_site,
  #alpha_true = alpha_true,
  #PCR_cycles = PCR_cycles,
  #tech_reps  = tech_reps,  
  #tech_reps_mock  = tech_reps_mock,
  #phi0_true  = phi0_true,
  #vary_PCR   = vary_PCR,
  #N_samp_vary_PCR = N_samp_vary_PCR, 
  #vary_PCR_cycles = vary_PCR_cycles, 
  env = env,
  mc = mc,
  
  # Realizations of simulations (input data)
  p_true = p_mock,
  p_samp_all = sample_data,
  p_mock_all = mock_data,
  
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

save(Output,file=paste0("./data/summarized_data/Hanfling_no_overdisp_even",".Rdata"))

## Plots
rstan::traceplot(stanMod,pars=c("lp__","alpha[1]","alpha[2]"),inc_warmup=FALSE)

plot(stanMod, par = "alpha")    

alphas <- data.frame(
  mean_log_alpha = stanMod_summary[["alpha"]][, 1],
  sp_idx = 1:length(stanMod_summary[["alpha"]][, 1])
)

posteriorProportions <- stanMod_summary[["int_samp_small"]][, c(1,4:8)] %>%
  as.data.frame() %>%
  mutate(community = rep(unique(env$community), each = N_species)) %>%
  mutate(sp_idx = rep(1:N_species, times = length(unique(env$community))))


(p11 <- env %>%
    group_by(tech_rep, community) %>%
    mutate(trueProp = start_conc_ng/sum(start_conc_ng)) %>%
    select(Species, sp_idx, community, tech_rep, Cycles, trueProp, propReads) %>% 
    right_join(posteriorProportions) %>%
    filter(mean > 0.005) %>% #omit rare things from plot
    ggplot(aes(x = Species, y = mean, col = Species)) +
    geom_point(aes(x = Species, y = propReads), color = "lightgrey") +
    geom_point(aes(x = Species, y = trueProp), color = "black", alpha = 0.5) +
    geom_point() +
    geom_segment(aes(x = Species, xend = Species, y = `25%`, yend = `75%`), size = .8) +
    geom_segment(aes(x = Species, xend = Species, y = `2.5%`, yend = `97.5%`), size = .4) +
    facet_grid(community~., drop = T)  +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("Proportion DNA") +
    guides(color = "none") +
    ggtitle("Hänfling Freshwater Fishes (mock community)"))

#ggsave(p11, filename = "Hanfling_Fish_Mock_even.pdf", width = 7, height = 5, units = "in")

