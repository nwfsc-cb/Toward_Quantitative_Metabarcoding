library(tidyverse)
library(compositions)
  select <- dplyr::select  #eliminate conflict w MASS::select vs. dplyr::select
library(rstan)
library(here)

ocean_dat <- readRDS("./data/Gold_Ocean/mifish_mock_community_data_North_Ocean.RDS") %>% 
  mutate(tech_rep = as.numeric(tech_rep),
         Cycles = as.numeric(Cycles)) %>% 
  filter(start_conc_ng > 0)
  

# species by community
sp.by.comm <- ocean_dat %>% filter(Cycles==39) %>% 
  group_by(species=ID_mifish,community) %>%
  summarise(conc = sum(start_conc_ng)) %>%
  mutate(conc2 = ifelse(conc>0,1,0)) %>%
  arrange(species,community) %>% group_by(community) %>%
  mutate(tot_conc = sum(conc),prop=conc/tot_conc)
  
sp.by.comm.slim <- sp.by.comm %>%
  dplyr::select(-conc) %>%
  pivot_wider(names_from = community,values_from = conc2,values_fill = 0) %>%
  pivot_longer(!species, names_to = "community",values_to= "present") %>% 
  mutate(present = as.factor(present))


# Make some grids of species to understand overlap.
A <- ggplot(sp.by.comm.slim) +
  geom_tile(aes(x=community,y=species,fill=present)) +
  scale_fill_viridis_d(option = "plasma",end = 0.75) +
  scale_y_discrete("") +
  theme(axis.text.x = element_text(angle=90) ) 
  

#REQUIRED FORMAT
#environmental (i.e., unknown) samples with the following column names:
  #community  -- the unique biological sample from which a set of sequences derived
  #Species   -- the biological species to which reads are assigned
  #nReads -- the number of reads assigned to that species in that community
  #tech_rep -- index of technical replicate, within community (e.g., if each community is sequenced 3 times, tech_rep will vary between 1 and 3)
  #Cycles -- number of PCR cycles the reaction underwent
  
  #mock community (i.e., known) samples with the following column names:
  #community  -- the unique biological sample from which a set of sequences derived
  #Species   -- the biological species to which reads are assigned
  #nReads -- the number of reads assigned to that species in that community
  #tech_rep -- index of technical replicate, within community (e.g., if each community is sequenced 3 times, tech_rep will vary between 1 and 3)
  #Cycles -- number of PCR cycles the reaction underwent
  #start_conc_ng -- starting concentration of each species in each community, in nanograms
  
  #set up environmental/unknown samples
  env <- ocean_dat %>%  
            filter(str_detect(community, "Skew", negate = F),Cycles==39) %>% 
            filter(str_detect(community, "Coastal", negate = T)) %>% 
            rename("Species" = "ID_mifish") %>% 
            group_by(community,tech_rep) %>% 
            mutate(propReads = nReads/sum(nReads), #calculate proportion of reads
                 totReads = sum(nReads)) %>%  #calculate total reads for community
            group_by(Species) %>% 
            mutate(totalSpeciesReads = sum(nReads)) %>%  
            add_tally(nReads > 0, name = "totalOccurrences") %>% 
            filter(totalSpeciesReads > 0) 
    
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
  
  #set up mock community samples
  mc <- ocean_dat %>%  
    filter(str_detect(community, "Even", negate = F),Cycles==39) %>%
    filter(str_detect(community, "Coastal", negate = T)) %>%
    rename("Species" = "ID_mifish")  #mock comm samples
    mc$Species[mc$Species == mostCommon] <- paste0("zRefSpecies_", mostCommon)   #make the most common species the reference species
  
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

p_mock_small <- mc %>% 
  filter(tech_rep == 1) %>% 
  select(community, sp_idx, start_conc_ng, Cycles) %>% 
  arrange(sp_idx) %>% 
  group_by(community) %>% 
  mutate(prop_conc = start_conc_ng/sum(start_conc_ng)) %>% 
  select(-start_conc_ng) %>%  # -Species) %>% 
  pivot_wider(names_from = "sp_idx", values_from = "prop_conc", values_fill = 0) %>% 
  ungroup() %>% 
  arrange(community)
  #select(-community, -Cycles)

  #calculate additive log ratios 
  alr_mock_true_prop <- p_mock[,4:(ncol(p_mock)-1)]*0
  for(i in 1:nrow(p_mock)){
    alr_mock_true_prop[i,] <- alr(p_mock[i,4:(ncol(p_mock))] + 1e-10)
  }
  alr_mock_true_prop[,N_species] <- 0 #adding explicit reference species column
  
  alr_mock_true_prop_small <- p_mock_small[,3:(ncol(p_mock_small)-1)]*0
  for(i in 1:nrow(p_mock_small)){
    alr_mock_true_prop_small[i,] <- alr(p_mock_small[i,3:(ncol(p_mock_small))] + 1e-10)
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
  sample_data = sample_data %>% select(-community,-Cycles,-tech_rep),
  mock_data   = mock_data  %>% select(-community,-Cycles,-tech_rep),
  
  # True proportions for mock community
  #mock_true_prop = p_mock_all %>% dplyr::select(contains("sp")),
  alr_mock_true_prop = alr_mock_true_prop,
  alr_mock_true_prop_small = alr_mock_true_prop_small,
  
  # vectors of PCR numbers
  N_pcr_samp = N_pcr_samp,
  N_pcr_mock = N_pcr_mock,
  
  # Design matrices: field samples
  N_b_samp_col = N_b_samp_col,
  model_matrix_b_samp = model_matrix_b_samp,
  model_matrix_b_samp_small = model_matrix_b_samp_small,
  model_vector_a_samp = model_vector_a_samp,
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
  "eta_samp",
  "eta_mock",
  "tau",
  "mu_samp",
  "mu_mock",
  "int_samp_small"
)

stan_init_f2 <- function(n.chain,N_species){#J_seb,K_seb){
  A <- list()
  for(i in 1:n.chain){
    A[[i]] <- list(
      # tau = runif(N_species-1,0.1,0.5),
      alpha_raw = runif(N_species-1,-0.5,0.5)
    )
  }
  return(A)
}

#########################################
#########################################
#Bayesian Estimation
N_CHAIN = 3
Warm = 1000
Iter = 1500
Treedepth = 13
Adapt_delta = 0.70

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanMod = stan(file = "./code/Core_stan_models/quant_metabar_multinom.stan" ,data = stan_data,
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
               sample_file = paste0("./tmpF.csv")
)

pars <- rstan::extract(stanMod, permuted = TRUE)
samp_params <- get_sampler_params(stanMod)

stanMod_summary <- list()
stanMod_summary[["alpha"]] <- summary(stanMod,pars="alpha")$summary
stanMod_summary[["tau"]] <- summary(stanMod,pars=c("tau"))$summary
stanMod_summary[["beta"]] <- summary(stanMod,pars="beta")$summary
stanMod_summary[["eta_samp_raw"]] <- summary(stanMod,pars="eta_samp")$summary
stanMod_summary[["eta_mock_raw"]] <- summary(stanMod,pars="eta_mock")$summary
stanMod_summary[["mu_samp"]] <- summary(stanMod,pars="mu_samp")$summary
stanMod_summary[["mu_mock"]] <- summary(stanMod,pars="mu_mock")$summary
stanMod_summary[["int_samp_small"]] <- summary(stanMod,pars="int_samp_small")$summary

Output <- list(
  ocean_dat = ocean_dat, # raw data from all observations, all communities.
  env = env,  #environmental data
  mc = mc, #mock data
  Species = unique(mc$Species),
  
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
  samp_params=samp_params, # Sampler information
  stanMod_summary = stanMod_summary # posterior summaries.
)

save(Output,file=paste0("./data/summarized_data/Fish_Oceanic_and_North_mock_2_crossValidate_Skew",".Rdata"))

posteriorProportions <- stanMod_summary[["int_samp_small"]][, c(1,4:8)] %>%
  as.data.frame() %>%
  mutate(comm_idx = rep(1:length(unique(env$community)), each = N_species)) %>%
  mutate(sp_idx = rep(1:N_species, times = length(unique(env$community))))

(p9 <- env %>%
    group_by(tech_rep, community) %>%
    mutate(trueProp = start_conc_ng/sum(start_conc_ng)) %>%
    right_join(posteriorProportions) %>%
    mutate(Species = gsub("zRefSpecies_", "", Species),
           BiolComm = gsub("_3.", "", community)) %>%
    filter(mean > 0.005) %>% #omit rare things from plot
    ggplot(aes(x = Species, y = mean, col = Species)) +
    geom_point(aes(x = Species, y = propReads), color = "lightgrey") +
    geom_point(aes(x = Species, y = trueProp), color = "black", alpha = 0.5) +
    geom_point() +
    geom_segment(aes(x = Species, xend = Species, y = `25%`, yend = `75%`), size = .6) +
    geom_segment(aes(x = Species, xend = Species, y = `2.5%`, yend = `97.5%`), size = .1) +
    #facet_grid(BiolComm~Cycles, drop = T)  +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("Proportion DNA") +
    guides(color = "none") +
    ggtitle("Pacific Ocean Fishes (mock community)"))

#ggsave(p9, filename = "Gold_Ocean_Mock.pdf", width = 7, height = 5, units = "in")




###################################################
###################################################
###################################################
#Likelihood

# M <- stan_model("../Core_stan_models/quant_metabar_multinom.stan")
# likMod <- optimizing(M, data=stan_data, iter=100000,draws=0,
#                      verbose=T,
#                      tol_obj = 1e-8,
#                      tol_rel_obj = 1e-8,
#                      tol_grad = 1e-8,
#                      tol_param=1e-8,
#                      tol_rel_grad=1e-8,
#                      algorithm="LBFGS")
# 
# 
# estProportions <- data.frame(mlEst = likMod$par[grepl("int_samp_small",names(likMod$par))]) %>%
#   mutate(community= rep(sample_data_small$community,N_species)) %>%
#   mutate(comm_idx = rep(1:length(unique(env$community)), times = N_species)) %>%
#   mutate(sp_idx = rep(1:N_species, each = length(unique(env$community))))
#   
# # merge in predictions to the env data.frame
# ppp <- env %>%
#   group_by(tech_rep, community) %>%
#   mutate(trueProp = start_conc_ng/sum(start_conc_ng)) %>%
#   left_join(.,estProportions) %>%
#   mutate(Species = gsub("zRefSpecies_", "", Species),
#          BiolComm = gsub("_3.", "", community)) %>%
#   mutate(propReads = ifelse(is.na(propReads),0,propReads))
#   
# p1 <-   ggplot(ppp) +
#   geom_point(aes(x = Species, y = propReads), color = "lightgrey") +
#   geom_point(aes(x = Species, y = trueProp), color = "black", alpha = 0.5) +
#   geom_point(data=ppp %>% filter(tech_rep==1),
#              aes(x = Species, y = mlEst, col = Species)) +
#   facet_grid(BiolComm~Cycles, drop = T) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   scale_y_continuous(trans="sqrt") +
#   ylab("Proportion DNA") +
#   guides(color = "none") +
#   ggtitle("ML Calibrating by Both Methods")
# 
# p1
# 
# # BREAKS <- c(0.001,0.01,0.05,0.10,0.20,0.30,0.40,0.50)
# # p2.raw <- ggplot(ppp) +
# #           geom_point(aes(x= trueProp,y=propReads,color=Species),alpha=0.5) +
# #           geom_abline(intercept=0,slope=1,color="red") +
# #           scale_y_continuous(trans="sqrt",limits=c(0,max(BREAKS)),breaks=BREAKS) +
# #           scale_x_continuous(trans="sqrt",limits=c(0,max(BREAKS)),breaks=BREAKS) +
# #           facet_wrap(~community)
# # p2.raw
# # 
# # p2.est <- ggplot(ppp) +
# #           geom_point(aes(x= trueProp,y=mlEst,color=Species),alpha=0.5) +
# #           geom_abline(intercept=0,slope=1,color="red") +
# #           scale_y_continuous(trans="sqrt",limits=c(0,max(BREAKS)),breaks=BREAKS) +
# #           scale_x_continuous(trans="sqrt",limits=c(0,max(BREAKS)),breaks=BREAKS) +
# #           facet_wrap(~community)  
# # p2.est  
# # 
# # (ppp <- env %>%
# #     group_by(tech_rep, community) %>%
# #     mutate(trueProp = start_conc_ng/sum(start_conc_ng)) %>%
# #     right_join(estProportions) %>%
# #     drop_na() %>% 
# #     mutate(Species = gsub("zRefSpecies_", "", Species),
# #            BiolComm = gsub("_3.", "", community)) %>%
# #     #filter(mlEst > 0.005) %>% #omit rare things from plot
# #     ggplot(aes(x = Species, y = mlEst, col = Species)) +
# #     geom_point(aes(x = Species, y = propReads), color = "lightgrey") +
# #     geom_point(aes(x = Species, y = trueProp), color = "black", alpha = 0.5) +
# #     geom_point() +
# #     #geom_segment(aes(x = Species, xend = Species, y = `25%`, yend = `75%`), size = .5) +
# #     facet_grid(BiolComm~Cycles, drop = T) +
# #     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
# #     ylab("Proportion DNA") +
# #     guides(color = "none") +
# #     ggtitle("ML Calibrating by Both Methods"))

###################################################
###################################################
###################################################


