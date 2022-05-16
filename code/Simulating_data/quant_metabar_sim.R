

# This script takes input from
p_true <- gtools::rdirichlet(N_site,alpha=rep(diri_alpha,N_species))

p_mock <- bind_cols(
  m_1=c(rep(1, N_species)),
  m_2=c(rep(1,floor(N_species/2)),rep(0,N_species- floor(N_species/2))),
  m_3=c(rep(0,floor(N_species/4)),
        rep(1,floor(N_species/2)),rep(0,N_species-floor(N_species/4)-floor(N_species/2))),
  m_4=c(rep(0,floor(N_species/2)),rep(1,N_species- floor(N_species/2))))
p_mock <- t(p_mock)
for(i in 1:nrow(p_mock)){
  p_mock[i,] <- p_mock[i,] / sum(p_mock[i,])
}

N_mock <- nrow(p_mock)

alr_mock_true_prop = p_mock[,1:(ncol(p_mock)-1)]*0
for(i in 1:nrow(p_mock)){
  alr_mock_true_prop[i,] = alr(p_mock[i,] + 1e-10)
}

###########################################################
# Make total amplicons (in log space)
###########################################################
p_amp <- t(t(log(p_true)) + PCR_cycles * log(1 + alpha_true))
p_exp_amp <- 0* p_amp
for(i in 1:ncol(p_amp)){
  p_exp_amp[,i] <- rnbinom(nrow(p_amp),mu=exp(p_amp[,i]),size=exp(phi0_true[i]))
}
max_p_exp_amp <- exp(log(p_true) + PCR_cycles * log(2))

# Check to make sure the observed do not exceed perfect amplification stochastically.
if(min(max_p_exp_amp - p_exp_amp ) < 0){ 
  STOP.CHECK<- "STOP, too many amplicons"
  print(STOP.CHECK)
  THESE <- which(p_exp_amp - max_p_exp_amp < 0,arr.ind=T)
}


#### For a constant community, visualize the proportion of amplicons for each species 
#### Across many replicates
# par(mfrow=c(4,4))
# for(i in 1:N_species){
#   hist(prop_amp[,i],breaks=100,main="")
#   title(paste("Species",i,"True prop =",p_true_vec[i]))
# }

# Declare a number of reads for each sample 
numb.samp.mid <- 100000
range.samp.frac <- 0.3 # what proportion (+/-) of numb.samp.mid do you observe for each sample

p_all <- NULL
for(j in 1:tech_reps){
  p_exp_amp <- 0* p_amp
  for(i in 1:ncol(p_amp)){
    p_exp_amp[,i] <- rnbinom(nrow(p_amp),mu=exp(p_amp[,i]),size=exp(phi0_true[i]))
  }
  # max_p_exp_amp <- exp(log(p_true) + PCR_cycles * log(2))
  # # Check to make sure the observed do not exceed perfect amplification stochastically.
  # if(min(max_p_exp_amp - p_exp_amp ) < 0){ 
  #   STOP.CHECK<- "STOP, too many amplicons"
  #   print(STOP.CHECK)
  #   THESE <- which(p_exp_amp - max_p_exp_amp < 0,arr.ind=T)
  # }
  # calculate proportion of reads in each sample for each species.
  prop_amp <- p_true * 0 
  for(i in 1:nrow(p_exp_amp)){
    prop_amp[i,] <- p_exp_amp[i,]/sum(p_exp_amp[i,])
  }
  p_samp <- prop_amp * 0 
  # number of observed reads among all species.
  n.samp <- round(runif(N_site,numb.samp.mid*(1-range.samp.frac),numb.samp.mid*(1+range.samp.frac)))
  for(i in 1:N_site){
    p_samp[i,] <- rmultinom(1,size=n.samp[i],prob = prop_amp[i,])
  }
  p_samp <- as.data.frame(p_samp)
  p_samp$site <- 1:N_site
  #A <- pivot_longer(p_samp,cols=contains(c("V")),names_to="sample",values_to="count")
  p_samp$tech_rep <- j
  if(j==1){p_all <- p_samp}
  if(j>1){p_all <- bind_rows(p_all,p_samp)}
}

p_all$N_pcr_samp <- PCR_cycles

# add a second batch of PCRs in "vary_PCR == TRUE" 
if(vary_PCR == TRUE){
  # Select the specified number of sites to run variable PCRs on
  rand_samp <- sample(1:N_site,N_samp_vary_PCR) %>% sort()
  p_all_2 <- NULL
  for(k in 1:length(vary_PCR_cycles)){

    p_amp_trim <- t(t(log(p_true[rand_samp,])) + vary_PCR_cycles[k] * log(1 + alpha_true))
    for(j in 1:tech_reps){
      p_exp_amp <- 0* p_amp_trim
      for(i in 1:ncol(p_amp_trim)){
        p_exp_amp[,i] <- rnbinom(nrow(p_amp_trim),mu=exp(p_amp_trim[,i]),size=exp(phi0_true[i]))
      }
      # max_p_exp_amp <- exp(log(p_true) + PCR_cycles * log(2))
      # # Check to make sure the observed do not exceed perfect amplification stochastically.
      # if(min(max_p_exp_amp - p_exp_amp ) < 0){ 
      #   STOP.CHECK<- "STOP, too many amplicons"
      #   print(STOP.CHECK)
      #   THESE <- which(p_exp_amp - max_p_exp_amp < 0,arr.ind=T)
      # }
      # calculate proportion of reads in each sample for each species.
      prop_amp <- p_true[rand_samp,] * 0 
      for(i in 1:nrow(p_exp_amp)){
        prop_amp[i,] <- p_exp_amp[i,]/sum(p_exp_amp[i,])
      }
      p_samp <- prop_amp * 0 
      # number of observed reads among all species.
      n.samp <- round(runif(N_site,numb.samp.mid*(1-range.samp.frac),numb.samp.mid*(1+range.samp.frac)))
      for(i in 1:N_samp_vary_PCR){
        p_samp[i,] <- rmultinom(1,size=n.samp[i],prob = prop_amp[i,])
      }
      p_samp <- as.data.frame(p_samp)
      p_samp$site <- rand_samp
      #A <- pivot_longer(p_samp,cols=contains(c("V")),names_to="sample",values_to="count")
      p_samp$tech_rep <- j + tech_reps
      p_samp$N_pcr_samp <- vary_PCR_cycles[k]
      if(j==1 & k==1){p_all_2 <- p_samp}
      if(j>1 | k > 1){p_all_2 <- bind_rows(p_all_2,p_samp)}
    } # end j loop
    
  }# end k loop
} # end vary PCR if statement.

p_all <- bind_rows(p_all,p_all_2)
p_all <- p_all %>% arrange(site,tech_rep)
head(p_all)
# print(STOP.CHECK)

p_samp_all <- p_all
p_samp_all <- p_samp_all %>% unite("site_tech",site:tech_rep,remove=F)
p_samp_all$N_seq <- rowSums(p_samp_all %>% dplyr::select(contains("V")))
colnames(p_samp_all)[grep("V",colnames(p_samp_all))] <- paste0("obs_sp_",1:N_species)

###########################################################
# Repeat for the mock communities
###########################################################
p_mock_amp <- t(t(log(p_mock + 1e-12)) + PCR_cycles * log(1 + alpha_true))
p_exp_mock_amp <- 0* p_mock_amp
for(i in 1:ncol(p_mock_amp)){
  p_exp_mock_amp[,i] <- rnbinom(nrow(p_mock_amp),mu=exp(p_mock_amp[,i]),size=exp(phi0_true[i]))
}
# max_p_exp_amp <- exp(log(p_true) + PCR_cycles * log(2))
# 
# # Check to make sure the observed do not exceed perfect amplification stochastically.
# if(min(max_p_exp_amp - p_exp_amp ) < 0){ 
#   STOP.CHECK<- "STOP, too many amplicons"
#   print(STOP.CHECK)
#   THESE <- which(p_exp_amp - max_p_exp_amp < 0,arr.ind=T)
# }

p_mock_all <- NULL
for(j in 1:tech_reps_mock){
  p_exp_mock_amp <- 0* p_mock_amp
  for(i in 1:ncol(p_mock_amp)){
    p_exp_mock_amp[,i] <- rnbinom(nrow(p_mock_amp),mu=exp(p_mock_amp[,i]),size=exp(phi0_true[i]))
  }
  # calculate proportion of reads in each sample for each species.
  prop_mock_amp <- p_mock_amp * 0 
  for(i in 1:nrow(p_exp_mock_amp)){
    prop_mock_amp[i,] <- p_exp_mock_amp[i,]/sum(p_exp_mock_amp[i,])
  }
  p_mock_samp <- prop_mock_amp * 0 
  # number of observed reads among all species.
  n.mock <- round(runif(N_mock,numb.samp.mid*(1-range.samp.frac),numb.samp.mid*(1+range.samp.frac)))
  for(i in 1:N_mock){
    p_mock_samp[i,] <- rmultinom(1,size=n.mock[i],prob = prop_mock_amp[i,])
  }
  p_mock_samp <- as.data.frame(p_mock_samp)
  p_mock_samp$site <- 1:N_mock
  #A <- pivot_longer(p_samp,cols=contains(c("V")),names_to="sample",values_to="count")
  
  p_mock_samp$tech_rep <- j
  if(j==1){p_mock_all <- p_mock_samp}
  if(j>1){p_mock_all <- bind_rows(p_mock_all,p_mock_samp)}
}

p_mock_all <- p_mock_all %>% arrange(site,tech_rep)
head(p_mock_all)

p_mock_all  <- p_mock_all %>% mutate(mock=paste0("m_",site))
p_mock_true <- p_mock %>% as.data.frame()
colnames(p_mock_true) <- paste0("sp_",1:N_species)
p_mock_true$mock <- as.character(rownames(p_mock_true))
alr_mock_true_prop <- alr_mock_true_prop %>% as.data.frame()

colnames(alr_mock_true_prop) <- paste0("alr_",1:(N_species-1))
alr_mock_true_prop[,paste0("alr_",N_species)] <- 0 
alr_mock_true_prop$mock <- as.character(rownames(alr_mock_true_prop))

p_mock_all  <- left_join(p_mock_all,p_mock_true) %>% left_join(.,alr_mock_true_prop)
p_mock_all <- p_mock_all %>% unite("mock_tech",mock:tech_rep,remove=F)
p_mock_all$N_pcr_mock <- PCR_cycles
p_mock_all$N_seq <- rowSums(p_mock_all %>% dplyr::select(contains("V")))
colnames(p_mock_all)[grep("V",colnames(p_mock_all))] <- paste0("obs_sp_",1:N_species)