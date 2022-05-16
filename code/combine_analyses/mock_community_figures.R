# Make figure comparing true and raw observed for three data sets
# hanfling
# Ocean fish data
# KW diet data

library(tidyverse)
library(grid)
library(gridExtra)
library(cowplot)
library(ggsci)

load("./data/summarized_data/Fish_Oceanic_and_North_mock_2_crossValidate_Skew.Rdata")
fish.dat <- Output
load("./data/summarized_data/kw_diet_output_fullmodel.Rdata")
orca.dat <- orca_poo
lake.dat <- read.csv("./data/Hanfling2016/Hanfling_12s_mock.csv")

########################################
##### These are the components of Figure 1 - version 1. highlight by species.
########################################

# Shared components 
 SIZE <- c(2,2,2,1)
 ALPHA <- c(1,1,1,0.25)
  
# Hanfling (lakes) First
  # Calculate fold-error
  lake.dat <- lake.dat %>% rename(true_prop=proportion) %>% 
                  mutate(log_obs_true = log(est_prop) - log(true_prop))
  # make a new column to enable highlighting specific species.
  lake.dat <- lake.dat %>% mutate(sp.new = "Other") %>%
                      mutate(sp.new = ifelse(Species=="T. tinca","T. tinca",sp.new)) %>%
                      mutate(sp.new = ifelse(Species=="A. alburnus","A. alburnus",sp.new)) %>%
                      mutate(sp.new = ifelse(Species=="G. cernua","G. cernua",sp.new))
  lake.dat$sp.new = factor(lake.dat$sp.new,levels=c("A. alburnus","T. tinca","G. cernua","Other") )

  lakes.p1a <- 
    ggplot(lake.dat) +
    #geom_point(aes(x=true_prop,y=est_prop,color=sp,shape=as.factor(community)),fill="white",alpha=0.5) + 
    geom_jitter(aes(x=true_prop,y=est_prop,color=Species,fill=Species),
                width=0.01, height=0, alpha=0.7,size=2) + 
    geom_abline(intercept=0,slope=1,color="red",linetype="dashed") +
    scale_color_viridis_d("Species",begin=0,end=0.9,option = "plasma") +
    scale_fill_viridis_d("Species",begin=0,end=0.9,option = "plasma") +
    #scale_shape_manual("Community",values=c(21,22)) +
    scale_x_continuous("True Proportion",limits = c(0,0.42),expand=c(0,0)) +
    scale_y_continuous("Estimated Proportion",limits = c(0,0.45),expand=c(0,0)) +
    #facet_wrap(~community) +
    theme_bw()

  lakes.p2a <-
    ggplot(lake.dat) +
    #geom_point(aes(x=true_prop,y=est_prop,color=sp,shape=as.factor(community)),fill="white",alpha=0.5) + 
    geom_jitter(aes(x=true_prop,y=log_obs_true,color=Species,fill=Species),
                width=0.005, height=0, alpha=0.7) + 
    geom_abline(intercept=0,slope=0,color="red",linetype="dashed") +
    scale_color_viridis_d("Species",begin=0,end=0.9,option = "plasma") +
    scale_fill_viridis_d("Species",begin=0,end=0.9,option = "plasma") +
    #scale_shape_manual("Community",values=c(21,22)) +
    scale_x_continuous("True Proportion",limits = c(0,0.38),expand=c(0,0)) +
    scale_y_continuous(expression("Fold error (log["*p[true]*"/"*p[obs]*"])")) + #,limits = c(0,0.45),expand=c(0,0)) +
    #facet_wrap(~community) +
    theme_bw()

  ####
  COL <- c(viridis::plasma(3,end=0.75),"black")

  lakes.p1b <- 
    ggplot(lake.dat) +
    #geom_point(aes(x=true_prop,y=est_prop,color=sp,shape=as.factor(community)),fill="white",alpha=0.5) + 
    geom_jitter(aes(x=true_prop,y=est_prop,color=sp.new,fill=sp.new,shape=sp.new,size=sp.new,alpha=sp.new),
                width=0.005, height=0) + 
    geom_abline(intercept=0,slope=1,color="black",linetype="dashed",alpha=0.5) +
    scale_color_manual("Species",values=COL) +
    scale_fill_manual("Species",values=COL) +
    scale_shape_manual("Species",values=c(24,22,23,21)) +
    scale_size_manual("Species",values=SIZE) +
    scale_alpha_manual("Species",values=ALPHA) +
    #scale_shape_manual("Community",values=c(21,22)) +
    scale_x_continuous("True Proportion",limits = c(0,0.42),expand=c(0,0)) +
    scale_y_continuous("Observed Proportion",limits = c(0,0.45),expand=c(0,0)) +
    #facet_wrap(~community) +
    theme_bw()
  
  lakes.p1b
  
  lakes.p2b <-
    ggplot(lake.dat) +
    #geom_point(aes(x=true_prop,y=est_prop,color=sp,shape=as.factor(community)),fill="white",alpha=0.5) + 
    geom_jitter(aes(x=true_prop,y=log_obs_true,color=sp.new,fill=sp.new,shape=sp.new,size=sp.new,alpha=sp.new),
                width=0.005, height=0) + 
    geom_abline(intercept=0,slope=0,color="black",linetype="dashed") +
    scale_color_manual("Species",values=COL) +
    scale_fill_manual("Species",values=COL) +
    scale_shape_manual("Species",values=c(24,22,23,21)) +
    scale_size_manual("Species",values=SIZE) +
    scale_alpha_manual("Species",values=ALPHA) +
    #scale_shape_manual("Community",values=c(21,22)) +
    scale_x_continuous("True Proportion",limits = c(0,0.38),expand=c(0,0)) +
    scale_y_continuous(expression("Fold error (log["*p[true]*"/"*p[obs]*"])")) + #,limits = c(0,0.45),expand=c(0,0)) +
    #facet_wrap(~community) +
    theme_bw()
  
  lakes.p2b
  
 
###########################################3
##### OK.  Repeat for Ocean fish data.
###########################################3

 ocean.dat <- fish.dat$ocean_dat %>% group_by(community,tech_rep,Cycles) %>% 
                  mutate(tot_conc = sum(start_conc_ng),
                         tot_Reads = sum(nReads)) %>%
                  mutate(true_prop = start_conc_ng/tot_conc,
                         est_prop  = nReads / tot_Reads,
                         log_true_est = log(true_prop) - log(est_prop)) %>%
                  rename(Species = ID_mifish) %>%
                  filter(!community == "North_Skew_2_39")
  
 ocean.dat <- ocean.dat %>% filter(Cycles==39)  
  
 ocean.dat.m <- ocean.dat %>% ungroup() %>% group_by(Species,community,true_prop) %>%
                  summarise(mean_est_prop = mean(est_prop),
                            min_est_prop = min(est_prop),
                            max_est_prop = max(est_prop)) %>%
                  mutate(mean_log_est_true = log(mean_est_prop) - log(true_prop))
                          
 # make a new column to enable highlighting specific species.
 ocean.dat.m <- ocean.dat.m %>% mutate(sp.new = "Other") %>%
                  mutate(sp.new = ifelse(Species=="Ophiodon elongatus","O. elongatus",sp.new)) %>%
                  mutate(sp.new = ifelse(Species=="Merluccius productus","M. productus",sp.new)) %>%
                  mutate(sp.new = ifelse(Species=="Engraulis mordax","E. mordax",sp.new))
 ocean.dat.m$sp.new = factor(ocean.dat.m$sp.new,levels=c("O. elongatus","M. productus","E. mordax","Other") )
 
 ## Plot
  ocean.p1b <- 
    ggplot(ocean.dat.m) +
    #geom_point(aes(x=true_prop,y=est_prop,color=sp,shape=as.factor(community)),fill="white",alpha=0.5) + 
    geom_jitter(aes(x=true_prop,y=mean_est_prop,color=sp.new,fill=sp.new,shape=sp.new,size=sp.new,alpha=sp.new),
                width=0.001, height=0) + 
    geom_abline(intercept=0,slope=1,color="black",linetype="dashed",alpha=0.5) +
    scale_color_manual("Species",values=COL) +
    scale_fill_manual("Species",values=COL) +
    scale_shape_manual("Species",values=c(24,22,23,21)) +
    scale_size_manual("Species",values=SIZE) +
    scale_alpha_manual("Species",values=ALPHA) +
    #scale_shape_manual("Community",values=c(21,22)) +
    scale_x_continuous("True Proportion",limits = c(0,0.25),expand=c(0,0)) +
    scale_y_continuous("Observed Proportion",limits = c(0,0.45),expand=c(0,0)) +
    #facet_wrap(~community) +
    theme_bw()
  
  ocean.p2b <-
    ggplot(ocean.dat.m) +
    #geom_point(aes(x=true_prop,y=est_prop,color=sp,shape=as.factor(community)),fill="white",alpha=0.5) + 
    geom_jitter(aes(x=true_prop,y=mean_log_est_true,color=sp.new,fill=sp.new,shape=sp.new,size=sp.new,alpha=sp.new),
                width=0.001, height=0) + 
    geom_abline(intercept=0,slope=0,color="black",linetype="dashed") +
    scale_color_manual("Species",values=COL) +
    scale_fill_manual("Species",values=COL) +
    scale_shape_manual("Species",values=c(24,22,23,21)) +
    scale_size_manual("Species",values=SIZE) +
    scale_alpha_manual("Species",values=ALPHA) +
    #scale_shape_manual("Community",values=c(21,22)) +
    scale_x_continuous("True Proportion",limits = c(0,0.25),expand=c(0,0)) +
    scale_y_continuous(expression("Fold error (log["*p[true]*"/"*p[obs]*"])")) + #,limits = c(0,0.45),expand=c(0,0)) +
    #facet_wrap(~community) +
    theme_bw()
  
  ocean.p2b
  
#############################################  
###### Repeat for Killer Whale Poop data.
#############################################  
  
  kw.dat <- orca.dat$dat.long.trim %>% filter(type_samp =="mock", true_prop>0) %>%
    rename(Species=sp) %>% mutate(log_obs_true = log(est_prop)-log(true_prop))
  
  kw.dat.m <- kw.dat %>% mutate(log_obs_true = log(est_prop)-log(true_prop)) %>% 
    group_by(community,Species,true_prop) %>% 
    summarise(mean_est_prop=mean(est_prop),
              min_est = min(est_prop),
              max_est=max(est_prop),
              mean_log_obs_true=mean(log_obs_true),
              min_log_obs_true = min(log_obs_true),
              max_log_obs_true = max(log_obs_true)) 
  
  kw.dat <- kw.dat %>% mutate(sp.new = "Other") %>%
    mutate(sp.new = ifelse(Species=="Oncorhynchus.mykiss","O. mykiss",sp.new)) %>%
    mutate(sp.new = ifelse(Species=="Oncorhynchus.tshawytscha","O. tshawytscha",sp.new)) %>%
    mutate(sp.new = ifelse(Species=="Clupea.pallasii","C. pallasii",sp.new))
  kw.dat$sp.new = factor(kw.dat$sp.new,levels=c("O. mykiss","O. tshawytscha","C. pallasii","Other") )
  
  kw.dat.m <- kw.dat.m %>% mutate(sp.new = "Other") %>%
    mutate(sp.new = ifelse(Species=="Oncorhynchus.mykiss","O. mykiss",sp.new)) %>%
    mutate(sp.new = ifelse(Species=="Oncorhynchus.tshawytscha","O. tshawytscha",sp.new)) %>%
    mutate(sp.new = ifelse(Species=="Clupea.pallasii","C. pallasii",sp.new))
  kw.dat.m$sp.new = factor(kw.dat.m$sp.new,levels=c("O. mykiss","O. tshawytscha","C. pallasii","Other") )
  
  
  ## Plot
  kw.p1b <- 
    ggplot() +
    #geom_point(aes(x=true_prop,y=est_prop,color=sp,shape=as.factor(community)),fill="white",alpha=0.5) + 
    geom_jitter(data= kw.dat.m, aes(x=true_prop,y=mean_est_prop,color=sp.new,fill=sp.new,shape=sp.new,size=sp.new,alpha=sp.new),
                width=0.001, height=0) + 
    geom_jitter(data= kw.dat, aes(x=true_prop,y=est_prop,color=sp.new,fill=sp.new,shape=sp.new),
                width=0.001, height=0,alpha=0.3,size=0.5) + 
    geom_abline(intercept=0,slope=1,color="black",linetype="dashed",alpha=0.5) +
    scale_color_manual("Species",values=COL) +
    scale_fill_manual("Species",values=COL) +
    scale_shape_manual("Species",values=c(24,22,23,21)) +
    scale_size_manual("Species",values=SIZE) +
    scale_alpha_manual("Species",values=ALPHA) +
    #scale_shape_manual("Community",values=c(21,22)) +
    scale_x_continuous("True Proportion",limits = c(0,0.42),expand=c(0,0)) +
    scale_y_continuous("Observed Proportion",limits = c(0,0.45),expand=c(0,0)) +
    #facet_wrap(~community) +
    theme_bw()
  
  kw.p2b <-
    ggplot(kw.dat.m) +
    #geom_point(aes(x=true_prop,y=est_prop,color=sp,shape=as.factor(community)),fill="white",alpha=0.5) + 
    geom_jitter(aes(x=true_prop,y=mean_log_obs_true,color=sp.new,fill=sp.new,shape=sp.new,size=sp.new,alpha=sp.new),
                width=0.001, height=0) +
    geom_jitter(data= kw.dat, aes(x=true_prop,y=log_obs_true,color=sp.new,fill=sp.new,shape=sp.new),
                width=0.001, height=0,alpha=0.3,size=0.75) + 
    geom_abline(intercept=0,slope=0,color="black",linetype="dashed") +
    scale_color_manual("Species",values=COL) +
    scale_fill_manual("Species",values=COL) +
    scale_shape_manual("Species",values=c(24,22,23,21)) +
    scale_size_manual("Species",values=SIZE) +
    scale_alpha_manual("Species",values=ALPHA) +
    #scale_shape_manual("Community",values=c(21,22)) +
    scale_x_continuous("True Proportion",limits = c(0,0.25),expand=c(0,0)) +
    scale_y_continuous(expression("Fold error (log["*p[true]*"/"*p[obs]*"])")) + #,limits = c(0,0.45),expand=c(0,0)) +
    #facet_wrap(~community) +
    theme_bw()
  
  kw.p2b
  
  ######## COMBINE PLOTS INTO ON 3x2 grid
  # combine into one figure for Hanfling.
  lay <- rbind(c(1,2,3),
               c(4,5,6),
               c(7,8,9))

  quartz(file="./plots/Mock_plot.pdf",type="pdf",dpi=600,height=8,width=7)  
    mock_plot <- grid.arrange(layout_matrix= lay,
               widths=c(1,1,0.4),
               lakes.p1b + theme(legend.position = "none")+ labs(subtitle="A",hjust=0,vjust=-1), 
               lakes.p2b + theme(legend.position = "none")+ labs(subtitle="B",hjust=0,vjust=-1),
               get_legend(lakes.p1b),
               ocean.p1b + theme(legend.position = "none")+ labs(subtitle="C",hjust=0,vjust=-1), 
               ocean.p2b + theme(legend.position = "none")+ labs(subtitle="D",hjust=0,vjust=-1),
               get_legend(ocean.p1b),
               kw.p1b + theme(legend.position = "none")+ labs(subtitle="E",hjust=0,vjust=-1), 
               kw.p2b + theme(legend.position = "none")+ labs(subtitle="F",hjust=0,vjust=-1),
               get_legend(kw.p1b)
                )
  
  dev.off()
  
  ########################################
  ##### Figure 1 - version 2. highlight by pulling out a single community.
  ########################################
  SIZE <- c(2,1)
  ALPHA <- c(0.75,0.25)
  COL <- c(viridis::plasma(1,end=0.75),"black")
  
  lake.comm <- "MC06"
  ocean.comm <- "Skew_Oceanic_2_39"
  ocean.comm.lab <- "Oceanic (Skew 2)"
  kw.comm <- 2
  kw.comm.lab <- "Salmonids"
  
  lake.dat <- lake.dat %>% mutate(comm.new = "Other") %>%
                  mutate(comm.new = ifelse(community ==lake.comm,lake.comm,comm.new))
  lake.dat.trim <- lake.dat %>% filter(community == lake.comm)

  ocean.dat.m <- ocean.dat.m %>% mutate(comm.new = "Other") %>%
    mutate(comm.new = ifelse(community ==ocean.comm,ocean.comm.lab,comm.new))
  ocean.dat.m$comm.new <-  factor(ocean.dat.m$comm.new,levels = c(ocean.comm.lab,"Other"))
  ocean.dat.trim <- ocean.dat.m %>% filter(community == ocean.comm)

  # Trim down the species so there are only 10 categories
  ocean.dat.trim <- ocean.dat.trim %>% mutate(sp.new = Species) %>%
                      mutate(sp.new = 
                               ifelse(Species %in% c("Stenobrachius leucopsarus", 
                                                     "Trachurus symmetricus",
                                                     "Triphoturus mexicanus"),
                                      "Other",sp.new))
  
  ocean.dat.trim <- ocean.dat.trim %>% 
    mutate(sp.new = case_when(sp.new =="Citharichthys sordidus" ~ "C. sordidus",    
                              sp.new =="Citharichthys stigmaeus" ~ "C. stigmaeus" ,
                              sp.new =="Citharichthys xanthostigma" ~ "C. xanthostigma",
                              sp.new =="Engraulis mordax" ~ "E. mordax",
                              sp.new =="Leuroglossus stilbius" ~ "L. stilbius",
                              sp.new =="Merluccius productus" ~ "M. productus",
                              sp.new =="Nannobrachium ritteri" ~ "N. ritteri",
                              sp.new =="Sardinops sagax" ~ "S. sagax",
                              sp.new =="Sebastes paucispinis" ~ "S. paucispinis",
                              TRUE ~ as.character(sp.new)))
    
  ocean.dat.trim$sp.new <- factor(ocean.dat.trim$sp.new, levels=c(ocean.dat.trim$sp.new %>% unique()))
  
  kw.dat <- kw.dat %>% mutate(comm.new = "Other") %>%
    mutate(comm.new = ifelse(community == kw.comm,kw.comm.lab,comm.new))
  kw.dat$comm.new <-  factor(kw.dat$comm.new,levels = c(kw.comm.lab,"Other"))
  
  kw.dat.m <- kw.dat.m %>% mutate(comm.new = "Other") %>%
    mutate(comm.new = ifelse(community == kw.comm,kw.comm.lab,comm.new))
  kw.dat.m$comm.new <-  factor(kw.dat.m$comm.new,levels = c(kw.comm.lab,"Other"))
  
  kw.dat.trim <- kw.dat.m %>% filter(community == kw.comm)
  
  kw.dat.trim <- kw.dat.trim %>% 
    mutate(sp.new = case_when(Species =="Clupea.pallasii" ~ "C. pallasii",    
                              Species =="Hippoglossus.stenolepis" ~ "H. stenolepis" ,
                              Species =="Oncorhynchus.tshawytscha" ~ "O. tshawytscha",
                              Species =="Ophiodon.elongatus" ~ "O. elongatus",
                              Species =="Oncorhynchus.gorbuscha" ~ "O. gorbuscha",
                              Species =="Oncorhynchus.keta" ~ "O. keta",
                              Species =="Oncorhynchus.kisutch" ~ "O. kisutch",
                              Species =="Oncorhynchus.mykiss" ~ "O. mykiss",
                              Species =="Oncorhynchus.nerka" ~ "O. nerka",
                              Species =="Salmo.salar" ~ "S. salar",
                              TRUE ~ as.character(sp.new)))
  
  
  # Make plots
  lakes.p1c <- 
    ggplot(lake.dat) +
    #geom_point(aes(x=true_prop,y=est_prop,color=sp,shape=as.factor(community)),fill="white",alpha=0.5) + 
    geom_jitter(aes(x=true_prop,y=est_prop,color=comm.new,
                    fill=comm.new,
                    shape=comm.new,
                    size=comm.new,
                    alpha=comm.new),
                width=0.005, height=0) + 
    geom_abline(intercept=0,slope=1,color="black",linetype="dashed",alpha=0.5) +
    scale_color_manual("Community",values=COL) +
    scale_fill_manual("Community",values=COL) +
    scale_shape_manual("Community",values=c(24,21)) +
    scale_size_manual("Community",values=SIZE) +
    scale_alpha_manual("Community",values=ALPHA) +
    #scale_shape_manual("Community",values=c(21,22)) +
    scale_x_continuous("True Proportion",limits = c(0,0.42),expand=c(0,0)) +
    scale_y_continuous("Observed Proportion",limits = c(0,0.45),expand=c(0,0)) +
     
    #facet_wrap(~community) +
    theme_bw() +
    theme(legend.position = c(0.01, 0.99),
          legend.justification = c("left", "top"))
  
  lakes.p1c
  
  
  ocean.p1c <- 
    ggplot(ocean.dat.m) +
    #geom_point(aes(x=true_prop,y=est_prop,color=sp,shape=as.factor(community)),fill="white",alpha=0.5) + 
    geom_jitter(aes(x=true_prop,y=mean_est_prop,
                    color=comm.new,
                    fill=comm.new,
                    shape=comm.new,
                    size=comm.new,
                    alpha=comm.new),
                width=0.001, height=0) + 
    geom_abline(intercept=0,slope=1,color="black",linetype="dashed",alpha=0.5) +
    scale_color_manual("Community",values=COL) +
    scale_fill_manual("Community",values=COL) +
    scale_shape_manual("Community",values=c(24,21)) +
    scale_size_manual("Community",values=SIZE) +
    scale_alpha_manual("Community",values=ALPHA) +
    #scale_shape_manual("Community",values=c(21,22)) +
    scale_x_continuous("True Proportion",limits = c(0,0.25),expand=c(0,0)) +
    scale_y_continuous("Observed Proportion",limits = c(0,0.35),expand=c(0,0)) +
    #facet_wrap(~community) +
    theme_bw() +
    theme(legend.position = c(0.01, 0.99),
          legend.justification = c("left", "top"))
  

  ocean.p1c
  

  kw.p1c <- 
    ggplot() +
    #geom_point(aes(x=true_prop,y=est_prop,color=sp,shape=as.factor(community)),fill="white",alpha=0.5) + 
    geom_point(data= kw.dat.m, aes(x=true_prop,y=mean_est_prop,
                                    color=comm.new,
                                    fill=comm.new,
                                    shape=comm.new,
                                    size=comm.new,
                                    alpha=comm.new),
                width=0.001, height=0) + 
    # geom_jitter(data= kw.dat, aes(x=true_prop,y=est_prop,
    #                              color=comm.new,
    #                              fill=comm.new,shape=comm.new),
    #              width=0.001, height=0,alpha=0.3,size=0.75) + 
    geom_abline(intercept=0,slope=1,color="black",linetype="dashed",alpha=0.5) +
    scale_color_manual("Community",values=COL) +
    scale_fill_manual("Community",values=COL) +
    scale_shape_manual("Community",values=c(24,21)) +
    scale_size_manual("Community",values=SIZE) +
    scale_alpha_manual("Community",values=ALPHA) +
    #scale_shape_manual("Community",values=c(21,22)) +
    scale_x_continuous("True Proportion",limits = c(0,0.42),expand=c(0,0)) +
    scale_y_continuous("Observed Proportion",limits = c(0,0.45),expand=c(0,0)) +
    #facet_wrap(~community) +
    theme_bw() +
    theme(legend.position = c(0.01, 0.99),legend.justification = c("left", "top"))
  
  kw.p1c
  
  ### Make stacked bar charts to accompany the above figures.
  LAB.stack <- c("True","Observed")
  
  lakes.p1.stack <-  
    ggplot(lake.dat.trim) + 
    geom_bar(aes(x=1,fill=Species,y=true_prop), position="fill", stat="identity") +
    geom_bar(aes(x=2,fill=Species,y=est_prop), position="fill", stat="identity") +
    scale_fill_jco(palette = c("default"),"MC06\nSpecies",alpha=1) +
    scale_x_continuous("",labels=LAB.stack,breaks=c(1,2)) +
    scale_y_continuous(expand=c(0,0)) +
    ylab("Proportion") +
    theme_bw() # +
    #theme(legend.position="none")
  
  lakes.p1.stack
  
  ocean.p1.stack <-  
    ggplot(ocean.dat.trim) + 
    geom_bar(aes(x=1,fill=sp.new,y=true_prop), position="fill", stat="identity") +
    geom_bar(aes(x=2,fill=sp.new,y=mean_est_prop), position="fill", stat="identity") +
    scale_fill_jco(palette = c("default"),"Ocean (Skew 2)\nSpecies",alpha=1) +
    scale_x_continuous("",labels=LAB.stack,breaks=c(1,2)) +
    scale_y_continuous(expand=c(0,0)) +
    ylab("Proportion") +
    theme_bw() 
  
  ocean.p1.stack
  
  kw.p1.stack <-  
    ggplot(kw.dat.trim) + 
    geom_bar(aes(x=1,fill=sp.new,y=true_prop), position="fill", stat="identity") +
    geom_bar(aes(x=2,fill=sp.new,y=mean_est_prop), position="fill", stat="identity") +
    scale_fill_jco(palette = c("default"),"Salmonid\nSpecies",alpha=1) +
    scale_x_continuous("",labels=LAB.stack,breaks=c(1,2)) +
    scale_y_continuous(expand=c(0,0)) +
    ylab("Proportion") +
    theme_bw() 
  
  kw.p1.stack
  
  ###############
  
  ######## COMBINE PLOTS INTO ON 3x2 grid
  # combine into one figure for Hanfling.
  lay <- rbind(c(1,2,3),
               c(4,5,6),
               c(7,8,9))
  
  
  
pointSize  = 0.01  
textSize   = 10
spaceLegend =0.1  

 mod.leg <-  
    # guides(shape = guide_legend(override.aes = list(size = pointSize)),
    #      color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
 
 textSize2   = 10
 spaceLegend = 0.5
  mod.leg2 <-  
   # guides(shape = guide_legend(override.aes = list(size = pointSize)),
   #      color = guide_legend(override.aes = list(size = pointSize))) +
   theme(legend.title = element_text(size = textSize2), 
         legend.text  = element_text(size = textSize2),
         legend.key.size = unit(spaceLegend, "lines"))
 
  
  quartz(file="./plots/Mock_plot_v2.pdf",type="pdf",dpi=600,height=9,width=7)  
  mock_plot_2 <- grid.arrange(layout_matrix= lay,
                            widths=c(1,1,0.4),
                            lakes.p1c + mod.leg + labs(subtitle="A",hjust=0),
                            lakes.p1.stack + theme(legend.position = "none")+ labs(subtitle="B",hjust=0,vjust=-1),
                            get_legend(lakes.p1.stack+ mod.leg2),
                            ocean.p1c + mod.leg + labs(subtitle="C",hjust=0), 
                            ocean.p1.stack + theme(legend.position = "none") + labs(subtitle="D",hjust=0),
                            get_legend(ocean.p1.stack+ mod.leg2) ,
                            kw.p1c + mod.leg + labs(subtitle="E",hjust=0),
                            kw.p1.stack + theme(legend.position = "none") + labs(subtitle="F",hjust=0),
                            get_legend(kw.p1.stack+ mod.leg2) 
  )
  dev.off()

  