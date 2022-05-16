# Compile fits for Ocean fish data and make figures of Ocean data for the publication.
library(tidyverse)
library(data.table)
library(gridExtra)
library(robCompositions)
library(grid)
library(cowplot)
library(ggsci)

# Read in the posteriors for the 
load("./data/summarized_data/kw_diet_output_raw.RData")
raw <- orca_poo  # Model without amplification adjustment
load("./data/summarized_data/kw_diet_output_fullmodel.RData")
mock <- orca_poo # Model with amplification adjustment

#########################################################
# Raw estimates, no amp bias adjust.
#########################################################

comm.sp.expand <- expand_grid(community=unique(mock$sample_data$community),
                    raw$sp.list %>% dplyr::select(id))

# This makes point estimates for each SRKW sample, each diet species.
raw.prop <- raw$stanMod_summary$int_samp_small %>% # rows are communities, columns are species. 
              as.data.frame() %>%
              bind_cols(.,comm.sp.expand) %>%
              left_join(.,raw$sp.list %>% dplyr::select(id,sp)) %>%
              left_join(.,raw$field_data_basic %>% dplyr::distinct(community,Pop)) %>% 
              filter(Pop != "ARKW")

COMM <- raw.prop %>% distinct(community)
# This uses the raw posterior, subsets for salmonids and calculates averages across samples:

# Make a general function for drawing from the posterior and calculating communities proportions.
ppd_prop <- function(model_matrix_b_samp_small,beta_matrix){ 
  prop_out <- NULL
  logit_val_samp_ppd <- matrix(0,dim(beta_matrix)[3],dim(beta_matrix)[2])
  for(i in 1:nrow(beta_matrix)){ # loop over posterior samples.
    for(n in 1:ncol(beta_matrix)){ # loop over retained species
      logit_val_samp_ppd[,n] = model_matrix_b_samp_small %*% beta_matrix[i,n,] ; 
    }
    temp <- (exp(logit_val_samp_ppd) / rowSums(exp(logit_val_samp_ppd)))%>% as.data.frame()
    temp$MCMC = i
    temp$comm.id <- 1:nrow(temp)
    
    prop_out <- bind_rows(prop_out,temp)
  }
  return(prop_out)
}

# apply function for all species in the mock community
raw_all_prop <- ppd_prop(raw$stan_data$model_matrix_b_samp_small,
                            beta_matrix=raw$pars$beta)
# add in species names as columns.
colnames(raw_all_prop)[1:max(raw$sp.list$id,na.rm=T)] <- raw$sp.list$sp[is.na(raw$sp.list$id)==F]
# merge in community identifiers.
comm.names <- data.frame(community=unique(raw$sample_data$community))
comm.names$comm.id <- 1:nrow(comm.names)
COLS <- raw$sp.list$sp
raw_all_prop_long <- left_join(raw_all_prop,comm.names) %>%
  filter(community %in% COMM$community) %>% # Include only communities for SRKW (derived above)
  pivot_longer(.,cols=all_of(COLS),names_to="species",values_to = "prop")

# Use the design matrices to generate the predicted proportions for each postrior draw
# but including only the salmon species.

# Identify which species are not salmonids
salmon.id <- sort(c(grep("Oncor",raw$sp.list$sp),grep("Salmo",raw$sp.list$sp)))
raw$sp.list$salmon.id[salmon.id]   <- 1:length(salmon.id)
                
# Subset the posterior samples to only include the salmon species.
raw_beta_salmon <- raw$pars$beta[,salmon.id,]

# apply function for just salmon species 
  raw_salmon_prop <- ppd_prop(raw$stan_data$model_matrix_b_samp_small,
                              beta_matrix=raw_beta_salmon)
  # add in species names as columns.
  colnames(raw_salmon_prop)[1:max(raw$sp.list$salmon.id,na.rm=T)] <- raw$sp.list$sp[is.na(raw$sp.list$salmon.id)==F]
  # merge in community identifiers.
  comm.names <- data.frame(community=unique(raw$sample_data$community))
  comm.names$comm.id <- 1:nrow(comm.names)
  COLS <- c(raw$sp.list$sp[is.na(raw$sp.list$salmon.id)==F])
  raw_salmon_prop_long <- left_join(raw_salmon_prop,comm.names) %>%
                        filter(community %in% COMM$community) %>% # Include only communities for SRKW (derived above)
                        pivot_longer(.,cols=all_of(COLS),names_to="species",values_to = "prop")
  
  #########################################################
  # Full model estimates with amp bias adjust.
  #########################################################

  # This makes point estimates for each SRKW sample, each diet species.
  mock.prop <- mock$stanMod_summary$int_samp_small %>% # rows are communities, columns are species. 
    as.data.frame() %>%
    bind_cols(.,comm.sp.expand) %>%
    left_join(.,mock$sp.list %>% dplyr::select(id,sp)) %>%
    left_join(.,mock$field_data_basic %>% dplyr::distinct(community,Pop)) %>% 
    filter(Pop != "ARKW")
  
  COMM <- mock.prop %>% distinct(community)
  # This uses the posterior, subsets for salmonids and calculates averages across samples:
  
    # apply function for all species in the mock community
  mock_all_prop <- ppd_prop(mock$stan_data$model_matrix_b_samp_small,
                           beta_matrix=mock$pars$beta)
  # add in species names as columns.
  colnames(mock_all_prop)[1:max(mock$sp.list$id,na.rm=T)] <- mock$sp.list$sp[is.na(mock$sp.list$id)==F]
  # merge in community identifiers.
  comm.names <- data.frame(community=unique(mock$sample_data$community))
  comm.names$comm.id <- 1:nrow(comm.names)
  COLS <- mock$sp.list$sp
  mock_all_prop_long <- left_join(mock_all_prop,comm.names) %>%
    filter(community %in% COMM$community) %>% # Include only communities for SRKW (derived above)
    pivot_longer(.,cols=all_of(COLS),names_to="species",values_to = "prop")
  
  # Use the design matrices to generate the predicted proportions for each postrior draw
  # but including only the salmon species.
  
  # Identify which species are not salmonids
  salmon.id <- sort(c(grep("Oncor",mock$sp.list$sp),grep("Salmo",mock$sp.list$sp)))
  mock$sp.list$salmon.id[salmon.id]   <- 1:length(salmon.id)
  
  # Subset the posterior samples to only include the salmon species.
  mock_beta_salmon <- mock$pars$beta[,salmon.id,]
  
  # apply function for just salmon species 
  mock_salmon_prop <- ppd_prop(mock$stan_data$model_matrix_b_samp_small,
                              beta_matrix=mock_beta_salmon)
  # add in species names as columns.
  colnames(mock_salmon_prop)[1:max(mock$sp.list$salmon.id,na.rm=T)] <- mock$sp.list$sp[is.na(mock$sp.list$salmon.id)==F]
  # merge in community identifiers.
  comm.names$comm.id <- 1:nrow(comm.names)
  COLS <- c(mock$sp.list$sp[is.na(mock$sp.list$salmon.id)==F])
  mock_salmon_prop_long <- left_join(mock_salmon_prop,comm.names) %>%
    filter(community %in% COMM$community) %>% # Include only communities for SRKW (derived above)
    pivot_longer(.,cols=all_of(COLS),names_to="species",values_to = "prop")
  
  ################################################################
  # These are the files of interest that can be used to create figures
  ################################################################
  
  # raw.prop # This is the summarized file with the estimates for each community (including all species)
  # raw_all_prop_long # file with all species from mock community
  # raw_salmon_prop_long # File with only salmonids
  
  # mock.prop # This is the summarized file with the estimates for each community (including all species)
  # mock_all_prop_long # file with all species from mock community
  # mock_salmon_prop_long # File with only salmonids
  
  raw_all_prop_long <- raw_all_prop_long %>%
    mutate(Species = "Other") %>%
    mutate(Species = ifelse(species=="Oncorhynchus.mykiss","O. mykiss",Species)) %>%
    mutate(Species = ifelse(species=="Oncorhynchus.keta","O. keta",Species)) %>%
    mutate(Species = ifelse(species=="Oncorhynchus.kisutch","O. kisutch",Species)) %>%
    mutate(Species = ifelse(species=="Oncorhynchus.tshawytscha","O. tshawytscha",Species)) %>%
    mutate(Species = ifelse(species=="Hippoglossus.stenolepis","H. stenolepis",Species))
  raw_all_prop_long$Species <- factor(raw_all_prop_long$Species,
                                      levels=
                                        c("H. stenolepis","O. keta","O. kisutch","O. mykiss",
                                          "O. tshawytscha","Other"))
  

  mock_all_prop_long <- mock_all_prop_long %>%
    mutate(Species = "Other") %>%
    mutate(Species = ifelse(species=="Oncorhynchus.mykiss","O. mykiss",Species)) %>%
    mutate(Species = ifelse(species=="Oncorhynchus.keta","O. keta",Species)) %>%
    mutate(Species = ifelse(species=="Oncorhynchus.kisutch","O. kisutch",Species)) %>%
    mutate(Species = ifelse(species=="Oncorhynchus.tshawytscha","O. tshawytscha",Species)) %>%
    mutate(Species = ifelse(species=="Hippoglossus.stenolepis","H. stenolepis",Species))
  mock_all_prop_long$Species <- factor(mock_all_prop_long$Species,
                                      levels=
                                        c("H. stenolepis","O. keta","O. kisutch","O. mykiss",
                                          "O. tshawytscha","Other"))
  
  # This should be identical to mock.prop... just good to spot check.
  mock_all_summary  <- mock_all_prop_long %>% 
                        group_by(species,community) %>% 
                        summarise(Mean=mean(prop),
                          Median=median(prop),
                          q.0.025 = quantile(prop,probs=c(0.025)),
                          q.0.05 = quantile(prop,probs=c(0.05)),
                          q.0.25 = quantile(prop,probs=c(0.25)),
                          q.0.75 = quantile(prop,probs=c(0.75)),
                          q.0.95 = quantile(prop,probs=c(0.95)),
                          q.0.975 = quantile(prop,probs=c(0.975))) %>%
                        mutate(model="mock")
  
  raw_all_summary  <- raw_all_prop_long %>% 
    group_by(species,community) %>% 
    summarise(Mean=mean(prop),
              Median=median(prop),
              q.0.025 = quantile(prop,probs=c(0.025)),
              q.0.05 = quantile(prop,probs=c(0.05)),
              q.0.25 = quantile(prop,probs=c(0.25)),
              q.0.75 = quantile(prop,probs=c(0.75)),
              q.0.95 = quantile(prop,probs=c(0.95)),
              q.0.975 = quantile(prop,probs=c(0.975))) %>%
    mutate(model="raw")
  
  mock_all_summary_groups  <- mock_all_prop_long %>% 
    group_by(Species,community) %>% 
    summarise(Mean=mean(prop),
              Median=median(prop),
              q.0.025 = quantile(prop,probs=c(0.025)),
              q.0.05 = quantile(prop,probs=c(0.05)),
              q.0.25 = quantile(prop,probs=c(0.25)),
              q.0.75 = quantile(prop,probs=c(0.75)),
              q.0.95 = quantile(prop,probs=c(0.95)),
              q.0.975 = quantile(prop,probs=c(0.975))) %>%
    mutate(model="mock")
  
  raw_all_summary_groups  <- raw_all_prop_long %>% 
    group_by(Species,community) %>% 
    summarise(Mean=mean(prop),
              Median=median(prop),
              q.0.025 = quantile(prop,probs=c(0.025)),
              q.0.05 = quantile(prop,probs=c(0.05)),
              q.0.25 = quantile(prop,probs=c(0.25)),
              q.0.75 = quantile(prop,probs=c(0.75)),
              q.0.95 = quantile(prop,probs=c(0.95)),
              q.0.975 = quantile(prop,probs=c(0.975))) %>%
    mutate(model="raw")
  
  
  
  # Use species groups
  mock_all_pop_avg_summary <- mock_all_prop_long %>% 
                                group_by(MCMC,Species) %>% 
                                summarise(Mean=mean(prop)) %>%
                                ungroup() %>% group_by(Species) %>%
                                summarise(grand.mean=mean(Mean),
                                          Median=median(Mean),
                                          q.0.025 = quantile(Mean,probs=c(0.025)),
                                          q.0.05 = quantile(Mean,probs=c(0.05)),
                                          q.0.25 = quantile(Mean,probs=c(0.25)),
                                          q.0.75 = quantile(Mean,probs=c(0.75)),
                                          q.0.95 = quantile(Mean,probs=c(0.95)),
                                          q.0.975 = quantile(Mean,probs=c(0.975))) %>%
                                mutate(model="mock")
  
  raw_all_pop_avg_summary <- raw_all_prop_long %>% 
    group_by(MCMC,Species) %>% 
    summarise(Mean=mean(prop)) %>%
    ungroup() %>% group_by(Species) %>%
    summarise(grand.mean=mean(Mean),
              Median=median(Mean),
              q.0.025 = quantile(Mean,probs=c(0.025)),
              q.0.05 = quantile(Mean,probs=c(0.05)),
              q.0.25 = quantile(Mean,probs=c(0.25)),
              q.0.75 = quantile(Mean,probs=c(0.75)),
              q.0.95 = quantile(Mean,probs=c(0.95)),
              q.0.975 = quantile(Mean,probs=c(0.975))) %>%
    mutate(model="raw")

  # mock_salmon_pop_avg_summary <- mock_salmon_prop_long %>%
  #   group_by(MCMC,Species) %>%
  #   summarise(Mean=mean(prop)) %>%
  #   ungroup() %>% group_by(species) %>%
  #   summarise(grand.mean=mean(Mean),
  #             Median=median(Mean),
  #             q.0.025 = quantile(Mean,probs=c(0.025)),
  #             q.0.05 = quantile(Mean,probs=c(0.05)),
  #             q.0.25 = quantile(Mean,probs=c(0.25)),
  #             q.0.75 = quantile(Mean,probs=c(0.75)),
  #             q.0.95 = quantile(Mean,probs=c(0.95)),
  #             q.0.975 = quantile(Mean,probs=c(0.975))) %>%
  #   mutate(model="mock")

  # raw_salmon_pop_avg_summary <- raw_salmon_prop_long %>%
  #   group_by(MCMC,Species) %>%
  #   summarise(Mean=mean(prop)) %>%
  #   ungroup() %>% group_by(Species) %>%
  #   summarise(grand.mean=mean(Mean),
  #             Median=median(Mean),
  #             q.0.025 = quantile(Mean,probs=c(0.025)),
  #             q.0.05 = quantile(Mean,probs=c(0.05)),
  #             q.0.25 = quantile(Mean,probs=c(0.25)),
  #             q.0.75 = quantile(Mean,probs=c(0.75)),
  #             q.0.95 = quantile(Mean,probs=c(0.95)),
  #             q.0.975 = quantile(Mean,probs=c(0.975))) %>%
  #   mutate(model="raw")

  # Use only September samples.
  mock_all_pop_avg_summary_sep <- mock_all_prop_long %>%
    filter(grepl("Sep|SEP",community)) %>%
    group_by(MCMC,Species) %>%
    summarise(Mean=mean(prop)) %>%
    ungroup() %>% group_by(Species) %>%
    summarise(grand.mean=mean(Mean),
              Median=median(Mean),
              q.0.025 = quantile(Mean,probs=c(0.025)),
              q.0.05 = quantile(Mean,probs=c(0.05)),
              q.0.25 = quantile(Mean,probs=c(0.25)),
              q.0.75 = quantile(Mean,probs=c(0.75)),
              q.0.95 = quantile(Mean,probs=c(0.95)),
              q.0.975 = quantile(Mean,probs=c(0.975))) %>%
    mutate(model="mock")
  
  raw_all_pop_avg_summary_sep <- raw_all_prop_long %>%
    filter(grepl("Sep|SEP",community)) %>%
    group_by(MCMC,Species) %>%
    summarise(Mean=mean(prop)) %>%
    ungroup() %>% group_by(Species) %>%
    summarise(grand.mean=mean(Mean),
              Median=median(Mean),
              q.0.025 = quantile(Mean,probs=c(0.025)),
              q.0.05 = quantile(Mean,probs=c(0.05)),
              q.0.25 = quantile(Mean,probs=c(0.25)),
              q.0.75 = quantile(Mean,probs=c(0.75)),
              q.0.95 = quantile(Mean,probs=c(0.95)),
              q.0.975 = quantile(Mean,probs=c(0.975))) %>%
    mutate(model="raw")
  
  #######################
  ###### OK.  Make stacked bar charts
  #######################
  pop_avg_summary <- bind_rows(mock_all_pop_avg_summary,
                               raw_all_pop_avg_summary )
  pop_avg_summary_sep <- bind_rows(mock_all_pop_avg_summary_sep,
                               raw_all_pop_avg_summary_sep )
  # pop_avg_summary_salmon <- bind_rows(mock_salmon_pop_avg_summary,
  #                              raw_salmon_pop_avg_summary )
  
  sample_avg_summary <- bind_rows(mock_all_summary,raw_all_summary)
  sample_avg_summary_sep <- bind_rows(mock_all_summary %>% filter(grepl("Sep|SEP",community)),
                                      raw_all_summary%>% filter(grepl("Sep|SEP",community)))
  sample_avg_summary_groups <- bind_rows(mock_all_summary_groups,raw_all_summary_groups)
  sample_avg_summary_groups_sep <- bind_rows(mock_all_summary_groups %>% filter(grepl("Sep|SEP",community)),
                                      raw_all_summary_groups %>% filter(grepl("Sep|SEP",community)))
  
    ### Make stacked bar charts to accompany the above figures.
  # This is the population average (all samples)
  pop.avg.stack <-  
    ggplot(pop_avg_summary) + 
    geom_bar(aes(x=model,fill=Species,y=grand.mean), position="fill", stat="identity") +
    #geom_bar(aes(x=2,fill=Species,y=est_prop), position="fill", stat="identity") +
    scale_fill_jco(palette = c("default"),"Species",alpha=1) +
    #scale_x_continuous("",labels=LAB.stack,breaks=c(1,2)) +
    scale_y_continuous(expand=c(0,0)) +
    ylab("Proportion") +
    theme_bw() # +
  #theme(legend.position="none")
  
  pop.avg.stack
  
  # pop.avg.salmon.stack <-  
  #   ggplot(pop_avg_summary_salmon) + 
  #   geom_bar(aes(x=model,fill=species,y=grand.mean), position="fill", stat="identity") +
  #   #geom_bar(aes(x=2,fill=Species,y=est_prop), position="fill", stat="identity") +
  #   scale_fill_jco(palette = c("default"),"species",alpha=1) +
  #   #scale_x_continuous("",labels=LAB.stack,breaks=c(1,2)) +
  #   scale_y_continuous(expand=c(0,0)) +
  #   ylab("Proportion") +
  #   theme_bw() # +
  # #theme(legend.position="none")
  
  sample.stack.community <-  
    ggplot(sample_avg_summary) + 
    geom_bar(aes(x=model,fill=species,y=Mean), position="fill", stat="identity") +
    #geom_bar(aes(x=2,fill=Species,y=est_prop), position="fill", stat="identity") +
    scale_fill_jco(palette = c("default"),"species",alpha=1) +
    #scale_x_continuous("",labels=LAB.stack,breaks=c(1,2)) +
    scale_y_continuous(expand=c(0,0)) +
    ylab("Proportion") +
    facet_wrap(~community) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))# +

  sample.stack.community.sep <-  
    ggplot(sample_avg_summary_sep) + 
    geom_bar(aes(x=model,fill=Species,y=Mean), position="fill", stat="identity") +
    #geom_bar(aes(x=2,fill=Species,y=est_prop), position="fill", stat="identity") +
    scale_fill_jco(palette = c("default"),"species",alpha=1) +
    #scale_x_continuous("",labels=LAB.stack,breaks=c(1,2)) +
    scale_y_continuous(expand=c(0,0)) +
    ylab("Proportion") +
    facet_wrap(~community) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))# +
  
  sample.stack.community.groups.sep <-  
    ggplot(sample_avg_summary_groups_sep) + 
    geom_bar(aes(x=model,fill=Species,y=Mean), position="fill", stat="identity") +
    #geom_bar(aes(x=2,fill=Species,y=est_prop), position="fill", stat="identity") +
    scale_fill_jco(palette = c("default"),"Species",alpha=1) +
    scale_x_discrete("",labels=c("Mock","None")) +
    scale_y_continuous(expand=c(0,0)) +
    ylab("Proportion") +
    facet_wrap(~community,nrow=2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))# +
  
  sample.stack.community.groups.sep
  
  ### Amplification Efficiencies.
  
  
  mock$alpha_clr$Species <- mock$alpha_clr$sp
  mock$alpha_clr$Species <-  gsub('[.]', " ", mock$alpha_clr$Species)
  mock$alpha_clr <- mock$alpha_clr %>%
    mutate(Group = "Other") %>%
    mutate(Group = ifelse(Species=="Oncorhynchus mykiss","O. mykiss",Group)) %>%
    mutate(Group = ifelse(Species=="Oncorhynchus keta","O. keta",Group)) %>%
    mutate(Group = ifelse(Species=="Oncorhynchus kisutch","O. kisutch",Group)) %>%
    mutate(Group = ifelse(Species=="Oncorhynchus tshawytscha","O. tshawytscha",Group)) %>%
    mutate(Group = ifelse(Species=="Hippoglossus stenolepis","H. stenolepis",Group))
    
  mock$alpha_clr$Group <- factor(mock$alpha_clr$Group,
                                  levels=c("H. stenolepis","O. keta","O. kisutch","O. mykiss",
                                          "O. tshawytscha","Other")) 
  mock$alpha_clr <-   mock$alpha_clr %>% arrange(Mean)
  mock$alpha_clr$Species <- factor(mock$alpha_clr$Species,
                                levels=mock$alpha_clr$Species) 
    
  p_clr_srkw <- ggplot(mock$alpha_clr) +
    geom_errorbarh(aes(xmin=q.0.25,xmax=q.0.75,y=Species),size=2,height=0) +
    geom_errorbarh(aes(xmin=q.0.025,xmax=q.0.975,y=Species),size=0.8,height=0) +
    geom_point(aes(x=Mean,y=Species,fill=Group,),size=3,shape=21) +
    geom_vline(xintercept=0,linetype="dashed") +
    scale_fill_jco(palette = c("default"),alpha=1) +
    scale_x_continuous("Amplification Efficiency\n(CLR)",limits=c(NA,0.06)) +
    scale_y_discrete(NULL) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.y = element_text(size=7))
  
  
  ###########################################
  ## MAKE FIGURES FOR PAPER.
  ###########################################333
  # Make a pretty plot for a subset of samples.
  FOCAL <- c("2017SEP10-12","2017SEP24-26","2017SEP22-06","2021Sep12-02")
  focal <- data.frame(community =FOCAL,numb=1:length(FOCAL))
  
  sample_avg_summary_groups_sep_mod <- left_join(sample_avg_summary_groups_sep,focal) %>%
          mutate(numb = ifelse(model=="mock",numb+0.2,numb),
                 numb = ifelse(model=="raw",numb-0.2,numb))
  
  BREAKS <- unique(sample_avg_summary_groups_sep_mod$numb) %>% sort()
  sample.stack.community.groups.sep.fin <-  
    ggplot(sample_avg_summary_groups_sep_mod %>% filter(community %in% focal$community)) + 
    geom_bar(aes(x=numb,fill=Species,y=Mean), position="fill", stat="identity",width=0.35) +
    #geom_bar(aes(x=2,fill=Species,y=est_prop), position="fill", stat="identity") +
    scale_fill_jco(palette = c("default"),"Species",alpha=1) +
    scale_x_continuous("",breaks = BREAKS, labels= rep(c("None","Mock"),4)) +
    scale_y_continuous(expand=c(0,0)) +
    ylab("Proportion") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))# +
  
  pop_avg_summary_sep$model <- factor(pop_avg_summary_sep$model,levels=c("raw","mock"))
  
  pop.avg.stack.sep <-  
    ggplot(pop_avg_summary_sep) + 
    geom_bar(aes(x=model,fill=Species,y=grand.mean), position="fill", stat="identity") +
    #geom_bar(aes(x=2,fill=Species,y=est_prop), position="fill", stat="identity") +
    scale_fill_jco(palette = c("default"),"Species",alpha=1) +
    scale_x_discrete("",labels=c("None","Mock")) +
    scale_y_continuous(expand=c(0,0)) +
    ylab("Proportion") +
    theme_bw() #+
    #theme(legend.position="none")
  

  
### COMBINE THE PLOTS
  
  lay <- rbind(c(1,3,4),
               c(2,2,4))
  
  pointSize  = 0.01  
  textSize   = 8
  spaceLegend =0.4  
  
  mod.leg <-  
    # guides(shape = guide_legend(override.aes = list(size = pointSize)),
    #      color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
  
  
  
  quartz(file="./plots/SRKW Diet.pdf",height=6,width=7,dpi=600,type="pdf")
  
  grid.arrange(layout_matrix= lay,
                              widths=c(1,0.35,1),
                            grobs=list(
                                  sample.stack.community.groups.sep.fin+  
                                      theme(legend.position="none",
                                            plot.margin=unit(c(0.2,0.05,-0.3,0.2),"cm")) +
                                      annotate(geom="text",x=0.55,y=0.95,label="A"),
                                  pop.avg.stack.sep +
                                    theme(legend.position="none") +
                                    annotate(geom="text",x=0.5,y=0.95,label="B"),
                                  get_legend(pop.avg.stack.sep + mod.leg),
                                  p_clr_srkw +theme(axis.text.y = element_text(size=7))+
                                    annotate(geom="text",x=-0.12,y=10.3,label="C"))
  )
  dev.off()
  
  quartz(file="./plots/SRKW Diet all sept.pdf",height=6,width=7,dpi=600,type="pdf")
  sample.stack.community.groups.sep
  dev.off()