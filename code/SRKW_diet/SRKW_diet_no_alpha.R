library(tidyverse)
library(compositions)
  select <- dplyr::select  #eliminate conflict w MASS::select vs. dplyr::select
library(rstan)
library(here)

#  YEAR <- 2021  
  
dat1 <- read.csv("./data/kw_diet/Quantmeta_species_model_mock.csv")
dat2 <- read.csv("./data/kw_diet/KWfecal_techrep_reads_allspecies.csv")
dat.true.mock <- read.csv("./data/kw_diet/Quantmeta_species_model_true_mock_prop.csv")

dat.mock <- dat1 %>% filter(type_samp=="mock")
dat.pcr.var <- dat1 %>% filter(type_samp=="PCR_var")

# Filter out all species which mock communities are non-zero in the experimental samples
# trim to include only non-labeling columns
A <- dat.true.mock[,8:ncol(dat.true.mock)]
B <- colSums(A) 

NOM <- c(colnames(A)[which(B>0)]) 

# Eliminate the reads associated with Orcas.
sp.list <- data.frame(sp=NOM[NOM != "Orcinus.orca"]) 
# Trim to include only Salmonids
# SAL <- sp.list$sp[c(grep("Oncor",sp.list$sp),grep("Salmo",sp.list$sp))]
# sp.list <- sp.list %>% filter(sp %in% SAL)

# Designate Oncorhynchus.tshawytscha as reference species.
sp.list$sp.mod <- sp.list$sp
sp.list <- sp.list %>% 
  mutate(sp.mod =ifelse(sp =="Oncorhynchus.tshawytscha","zzz_Ref_Oncorhynchus.tshawytscha",sp.mod)) %>% 
  arrange(sp.mod) %>%
  mutate(id=1:nrow(.)) %>%
  mutate(sp.lab = paste0("sp_",id))

N_species = nrow(sp.list)

#N_species <- nrow(sp.list)
sp.list$sp.lab <- factor(sp.list$sp.lab,levels=sp.list$sp.lab)
COLS <- c(sp.list$sp,"year","miseq_run","type_samp","community","pcr_cycles","sample_id","tech_rep")
         
# combine data into one file, pivot to long form, add in true proportions for mock communities.
dat1$community <- as.character(dat1$community)
dat.long <- dat1 %>% dplyr::select(all_of(COLS)) %>%
                  pivot_longer(all_of(sp.list$sp),names_to = "sp", values_to = "counts") %>%
          bind_rows(.,
                    dat2 %>% dplyr::select(all_of(COLS)) %>%
                    pivot_longer(all_of(sp.list$sp),names_to = "sp", values_to = "counts")
          )

dat.true.long <- dat.true.mock %>% dplyr::select(all_of(COLS)) %>% 
                  pivot_longer(all_of(sp.list$sp),names_to = "sp", values_to = "true_prop") %>%
                  mutate(community = as.character(community)) %>%
                  ungroup() %>% group_by(sample_id) %>%
                  mutate(sum_true_prop=sum(true_prop),
                         true_prop_mod = true_prop/sum_true_prop)

dat.long <- left_join(dat.long,dat.true.long)
dat.long <- dat.long %>% rename(Cycles=pcr_cycles) %>% 
                group_by(sample_id, Cycles, tech_rep) %>% 
                mutate(tot_count = sum(counts)) %>%
                ungroup() %>%
                filter(tot_count > 5000) %>%
                mutate(est_prop = counts / tot_count) %>%
                left_join(.,sp.list)


# This removes non-zero counts from mock communities
# where the true value is zero
dat.long.trim <- dat.long %>%
                mutate(counts = ifelse(type_samp=="mock" & 
                                         true_prop == 0, 0, counts)) %>%
                group_by(sample_id, Cycles, tech_rep) %>% 
                mutate(tot_count = sum(counts)) %>%
                mutate(est_prop = counts / tot_count)


####### Use only the terminal PCR cycle. (32) and samples that have replicates.
# pivot back to wide form.
sample_data <- dat.long.trim %>% 
                    filter(!type_samp %in% c("mock","PCR_var"),
                            !community %in% c("NGOS_2021_17","2021Sep11-02"),  # this is the only sample that has a singleton observation.
                            Cycles==32) %>% 
                dplyr::select(community,sample_id,Cycles,tech_rep,sp.lab,counts) %>% 
                pivot_wider(names_from = sp.lab,values_from = counts)

sample_data_small <- sample_data %>% filter(tech_rep == 1)

mock_data <- dat.long.trim %>% filter(type_samp == "mock",
                                      Cycles==32) %>% 
                    dplyr::select(community,sample_id,Cycles,tech_rep,sp.lab,counts) %>% 
                    pivot_wider(names_from = sp.lab,values_from = counts)

mock_data_small <- mock_data %>% filter(tech_rep == 1)
dat.mock.true.wide <- dat.long.trim %>% filter(type_samp == "mock",
                                               Cycles==32) %>%  
                    dplyr::select(community,sample_id,Cycles,tech_rep,sp.lab,true_prop ) %>%
                    #mutate(true_prop_mod = true_prop + 1e-12) %>% dplyr::select(-true_prop) %>%
                    pivot_wider(names_from = sp.lab,values_from = true_prop)

# Do ALR transform

#calculate additive log ratios 
alr_mock_true_prop <- dat.mock.true.wide[,5:(ncol(dat.mock.true.wide)-1)]*0
for(i in 1:nrow(dat.mock.true.wide)){
  alr_mock_true_prop[i,] <- alr(dat.mock.true.wide[i,5:(ncol(dat.mock.true.wide))] + 1e-10)
}
alr_mock_true_prop[,N_species] <- 0 #adding explicit reference species column

# FOR SMALL IF I WANT TO USE THEM.
# alr_mock_true_prop_small <- p_mock_small[,3:(ncol(p_mock_small)-1)]*0
# for(i in 1:nrow(p_mock_small)){
#   alr_mock_true_prop_small[i,] <- alr(p_mock_small[i,3:(ncol(p_mock_small))] + 1e-10)
# }
# alr_mock_true_prop_small[,N_species] <- 0 

# DESIGN MATRICES
# mock communities first
# species compositions (betas)
# use mock_data  

N_pcr_mock <- mock_data$Cycles

if(length(unique(mock_data$community))==1){
  formula_b <- Cycles ~ 1  # what is on the left side of the equation doesn't matter.
} else {
  formula_b <- Cycles ~ as.factor(community) # what is on the left side of the equation doesn't matter.
}
model_frame <- model.frame(formula_b, mock_data)
model_matrix_b_mock <- model.matrix(formula_b, model_frame)

model_frame <- model.frame(formula_b, mock_data_small)
model_matrix_b_mock_small <- model.matrix(formula_b, model_frame)

# efficiencies (alphas)
formula_a <- community ~ Cycles -1
model_frame <- model.frame(formula_a, mock_data)
model_vector_a_mock <- model.matrix(formula_a, model_frame) %>% as.numeric()

model_frame <- model.frame(formula_a, mock_data_small)
model_vector_a_mock_small <- model.matrix(formula_a, model_frame) %>% as.numeric()

N_obs_mock       <- nrow(mock_data)
N_b_mock_col     <- ncol(model_matrix_b_mock)  

# unknown communities second
# species compositions (betas)

# use sample_data
N_pcr_samp <- sample_data$Cycles

if(length(unique(sample_data$community))==1){
  formula_b <- Cycles ~ 1  
} else {
  formula_b <- Cycles ~ as.factor(community)
}
model_frame <- model.frame(formula_b, sample_data)
model_matrix_b_samp <- model.matrix(formula_b, model_frame)

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

# Put all of the data into a list.
stan_data <- list(
  N_species = N_species,   # Number of species in data
  N_obs_samp = N_obs_samp, # Number of observed samples 
  N_obs_mock = N_obs_mock, # Number of observed mock samples
  N_obs_samp_small = N_obs_samp_small, # Number of observed samples 
  
  # Observed data of community matrices
  sample_data = sample_data %>% ungroup() %>% dplyr::select(contains("sp_")),
  mock_data   = mock_data %>% ungroup() %>% dplyr::select(contains("sp_")),
  
  # True proportions for mock community
  #mock_true_prop = p_mock_all %>% dplyr::select(contains("sp")),
  alr_mock_true_prop = alr_mock_true_prop,
  
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
  #"alpha",
  "beta",
  "tau",
  "mu_samp",
  #"mu_mock",
  "int_samp_small"
)

stan_init_f2 <- function(n.chain,N_species){#J_seb,K_seb){
  A <- list()
  for(i in 1:n.chain){
    A[[i]] <- list(
      tau = runif(N_species-1,0.1,0.5)
      #alpha_raw = runif(N_species-1,-0.5,0.5)
    )
  }
  return(A)
}

### Bayesian Estimation.
#########################################
#########################################
#Bayesian Estimation
N_CHAIN = 3
Warm = 1000
Iter = 1500
Treedepth = 13
Adapt_delta = 0.7

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanMod = stan(file = "./code/Core_stan_models/quant_metabar_no_mock_no_alpha.stan" ,data = stan_data,
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

pars <- rstan::extract(stanMod, permuted = TRUE)
samp_params <- get_sampler_params(stanMod)


stanMod_summary <- list()
#stanMod_summary[["beta"]] <- summary(stanMod,pars="beta")$summary
stanMod_summary[["tau"]] <- summary(stanMod,pars="tau")$summary
stanMod_summary[["mu_samp"]] <- summary(stanMod,pars="mu_samp")$summary
#stanMod_summary[["mu_mock"]] <- summary(stanMod,pars="mu_mock")$summary
stanMod_summary[["int_samp_small"]] <- summary(stanMod,pars="int_samp_small")$summary



# Merge in predictions for unknown communities.

orca_poo <- list()

orca_poo$sp.list <- sp.list
orca_poo$sample_data <- sample_data
orca_poo$field_data_basic <- dat2
orca_poo$dat.long.trim <- dat.long.trim
orca_poo$stan_data <- stan_data
orca_poo$stanMod <- stanMod
orca_poo$pars <- pars
orca_poo$stanMod_summary <- stanMod_summary

save(orca_poo, file="./data/summarized_data/kw_diet_output_raw.RData")












##
# Mock communities only
##
###
ggplot(dat.long %>% filter(type_samp=="mock",true_prop>0)) +
    #geom_point(aes(x=true_prop,y=est_prop,color=sp,shape=as.factor(community)),fill="white",alpha=0.5) + 
    geom_point(data=dat.m %>% filter(type_samp=="mock",true_prop>0),
                  aes(x=true_prop,y=mean_est_prop,color=sp,fill=sp,shape=as.factor(community)),alpha=1,size=2) + 
    geom_errorbar(data=dat.m %>% filter(type_samp=="mock",true_prop>0),
             aes(x=true_prop,ymin=min_est,ymax=max_est,color=sp),alpha=0.5,width=0) + 
    geom_abline(intercept=0,slope=1,color="red",linetype="dashed") +
    scale_color_viridis_d("Species",begin=0,end=0.9,option = "plasma") +
    scale_fill_viridis_d("Species",begin=0,end=0.9,option = "plasma") +
    scale_shape_manual("Community",values=c(21,22)) +
    scale_x_continuous("True Proportion",limits = c(0,0.42),expand=c(0,0)) +
    scale_y_continuous("Estimated Proportion",limits = c(0,0.45),expand=c(0,0)) +
    #facet_wrap(~community) +
    theme_bw()

## Log-ratio of observed-predicted for mock communities.


ggplot(dat.long %>% filter(type_samp=="mock",true_prop>0)) +
    #geom_point(aes(y=log_obs_true,x=true_prop,,color=sp,fill=sp,shape=as.factor(community)),alpha=1,size=1,fill="white") +
    geom_point(data=dat.m %>% filter(type_samp=="mock",true_prop>0),
             aes(x=true_prop,y=mean_log_obs_true,color=sp,fill=sp,shape=as.factor(community)),alpha=1,size=2) + 
    geom_errorbar(data=dat.m %>% filter(type_samp=="mock",true_prop>0),
                aes(x=true_prop,ymin=min_log_obs_true,ymax=max_log_obs_true,color=sp),alpha=0.5,width=0) +
    geom_abline(intercept=0,slope=0,color="red",linetype="dashed") +
    scale_color_viridis_d("Species",begin=0,end=0.75,option = "plasma") +
    scale_fill_viridis_d("Species",begin=0,end=0.75,option = "plasma") +
    scale_shape_manual("Community",values=c(21,22)) +
    scale_x_continuous("True Proportion",limits = c(0,0.42),expand=c(0,0)) +
    scale_y_continuous("Fold error") + #,limits = c(0,0.45),expand=c(0,0)) +
    #facet_wrap(~community) +
    theme_bw()



# Raw observations vs. modified model estimates.


