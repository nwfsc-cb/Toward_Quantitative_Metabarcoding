# Compile fits for Hanfling fish data and make figures of Hanfling fits for the publication.
library(tidyverse)
library(data.table)
library(gridExtra)
library(robCompositions)
library(grid)
library(cowplot)
library(ggsci)


# Read in the posteriors for the 
load("./data/summarized_data/Hanfling_no_overdisp_no_calibration.Rdata")
raw <- Output
load("./data/summarized_data/Hanfling_no_overdisp_odd.Rdata")
mock1 <- Output # This is the cross validation that uses all of the even communities
load("./data/summarized_data/Hanfling_no_overdisp_even.Rdata")
mock2 <- Output # This is the cross validation that uses only the 39 cycle even communities

# mock1 - predict communities even using odd as mock communities.
# mock2 - predict communities odd using even as mock communities.

#########################################################
# Raw estimates from Reads
#########################################################
# summarize raw estimates from reads for each species.
raw.raw <- raw$env %>% group_by(community,Cycles,tech_rep) %>%
  mutate(sum.ng = sum(start_conc_ng),
         true.prop = start_conc_ng / sum.ng) %>%
  ungroup() %>%
  group_by(Species,community,Cycles,true.prop) %>%
  summarise(simple.Mean=mean(propReads),
            simple.N = length(tech_rep)) %>%
  replace_na(list(raw.Mean=0,raw.SD=0,raw.SE=0))

# extract predicted proportions from the posterior
COM <- data.frame(community = levels(raw$env$community %>% as.factor()))
COM$comm_idx <- 1:nrow(COM)
SP  <- raw$env %>% distinct(Species,sp_idx) %>% as.data.frame()

# These are the predicted intercepts for the posteriors
beta_posterior <- raw$stanMod_summary[["int_samp_small"]][, c(1,4:8)]
colnames(beta_posterior) <- paste0("raw.",substr(colnames(beta_posterior),1,nchar(colnames(beta_posterior))-1))
colnames(beta_posterior)[1] <- "raw.mean"
beta_posterior <- as.data.frame(beta_posterior)

raw.post <-expand.grid(comm_idx = COM$comm_idx,sp_idx =SP$sp_idx) %>% 
  arrange(comm_idx,sp_idx) %>% 
  left_join(.,COM) %>% 
  left_join(.,SP) %>% 
  bind_cols(.,beta_posterior)

# Combine the raw estimates and posterior estimates
raw.all <- full_join(raw.raw,raw.post)

#########################################################
# Mock1
#########################################################
# summarize raw estimates from reads for each species.
mock1.raw <- mock1$env %>% group_by(community,Cycles,tech_rep) %>%
      mutate(sum.ng = sum(start_conc_ng),
             true.prop = start_conc_ng / sum.ng) %>%
      ungroup() %>%
      group_by(Species,community,Cycles,true.prop) %>%
  summarise(simple.Mean=mean(propReads),
            simple.N = length(tech_rep)) %>%
  replace_na(list(raw.Mean=0,raw.SD=0,raw.SE=0))

# extract predicted proportions from the posterior
COM <- data.frame(community = levels(mock1$env$community %>% as.factor()))
COM$comm_idx <- 1:nrow(COM)
SP  <- mock1$env %>% distinct(Species,sp_idx) %>% as.data.frame()

# These are the predicted intercepts for the posteriors
beta_posterior <- mock1$stanMod_summary[["int_samp_small"]][, c(1,4:8)]
colnames(beta_posterior) <- paste0("mock.",substr(colnames(beta_posterior),1,nchar(colnames(beta_posterior))-1))
colnames(beta_posterior)[1] <- "mock.mean"
beta_posterior <- as.data.frame(beta_posterior)

mock1.post <-expand.grid(comm_idx = COM$comm_idx,sp_idx =SP$sp_idx) %>% 
    arrange(comm_idx,sp_idx) %>% 
    left_join(.,COM) %>% 
    left_join(.,SP) %>% 
    bind_cols(.,beta_posterior)

# Combine the raw estimates and posterior estimates
mock1.all <- full_join(mock1.raw,mock1.post)

#########################################################
# Mock2
#########################################################
# summarize raw estimates from reads for each species.
mock2.raw <- mock2$env %>% group_by(community,Cycles,tech_rep) %>%
  mutate(sum.ng = sum(start_conc_ng),
         true.prop = start_conc_ng / sum.ng) %>%
  ungroup() %>%
  group_by(Species,community,Cycles,true.prop) %>%
  summarise(simple.Mean=mean(propReads),
            simple.N = length(tech_rep)) %>%
  replace_na(list(raw.Mean=0,raw.SD=0,raw.SE=0))

# extract predicted proportions from the posterior
COM <- data.frame(community = levels(mock2$env$community %>% as.factor()))
COM$comm_idx <- 1:nrow(COM)
SP  <- mock2$env %>% distinct(Species,sp_idx) %>% as.data.frame()

# These are the predicted intercepts for the posteriors
beta_posterior <- mock2$stanMod_summary[["int_samp_small"]][, c(1,4:8)]
colnames(beta_posterior) <- paste0("mock.",substr(colnames(beta_posterior),1,nchar(colnames(beta_posterior))-1))
colnames(beta_posterior)[1] <- "mock.mean"
beta_posterior <- as.data.frame(beta_posterior)

mock2.post <-expand.grid(comm_idx = COM$comm_idx,sp_idx =SP$sp_idx) %>% 
  arrange(comm_idx,sp_idx) %>% 
  left_join(.,COM) %>% 
  left_join(.,SP) %>% 
  bind_cols(.,beta_posterior)

# Combine the raw estimates and posterior estimates
mock2.all <- full_join(mock2.raw,mock2.post)

# Combine the predictions into one file for convenience.
mock.all <- bind_rows(mock1.all,mock2.all)

########################################################################
########################################################################
########################################################################
########################################################################

# Combine mock and pcr_var results with raw reads.
result.dat <- left_join(mock.all,
                        raw.all %>% dplyr::select(-comm_idx,-sp_idx))

# make a distinct factor for large true proportions and for small true proportions
result.dat <- result.dat %>% 
                mutate(true.cat= case_when(true.prop < 0.1 ~ "small",
                                          true.prop >= 0.1 & true.prop < 0.2  ~ "medium",
                                          true.prop >= 0.2   ~ "large"))
                                      
# Add.offset
# pull out just ocean.skew.dat for plotting.
spread=0.15
mc01.dat <- result.dat %>% 
  filter(true.prop > 0,community=="MC01")
mc01.dat <- bind_cols(mc01.dat, 
                             data.frame(offset= seq(-spread,spread,length.out=nrow(mc01.dat))))
mc02.dat <- result.dat %>% 
  filter(true.prop > 0,community=="MC02")
mc02.dat <- bind_cols(mc02.dat, 
                      data.frame(offset= seq(-spread,spread,length.out=nrow(mc02.dat))))
mc03.dat <- result.dat %>% 
  filter(true.prop > 0,community=="MC03")
mc03.dat <- bind_cols(mc03.dat, 
                      data.frame(offset= seq(-spread,spread,length.out=nrow(mc03.dat))))
mc04.dat <- result.dat %>% 
  filter(true.prop > 0,community=="MC04")
mc04.dat <- bind_cols(mc04.dat, 
                      data.frame(offset= seq(-spread,spread,length.out=nrow(mc04.dat))))
mc05.dat <- result.dat %>% 
  filter(true.prop > 0,community=="MC05")
mc05.dat <- bind_cols(mc05.dat, 
                      data.frame(offset= seq(-spread,spread,length.out=nrow(mc05.dat))))
mc06.dat <- result.dat %>% 
  filter(true.prop > 0,community=="MC06")
mc06.dat <- bind_cols(mc06.dat, 
                      data.frame(offset= seq(-spread,spread,length.out=nrow(mc06.dat))))
mc07.dat <- result.dat %>% 
  filter(true.prop > 0,community=="MC07")
mc07.dat <- bind_cols(mc07.dat, 
                      data.frame(offset= seq(-spread,spread,length.out=nrow(mc07.dat))))
mc08.dat <- result.dat %>% 
  filter(true.prop > 0,community=="MC08")
mc08.dat <- bind_cols(mc08.dat, 
                      data.frame(offset= seq(-spread,spread,length.out=nrow(mc08.dat))))
mc09.dat <- result.dat %>% 
  filter(true.prop > 0,community=="MC09")
mc09.dat <- bind_cols(mc09.dat, 
                      data.frame(offset= seq(-spread,spread,length.out=nrow(mc09.dat))))
mc10.dat <- result.dat %>% 
  filter(true.prop > 0,community=="MC10")
mc10.dat <- bind_cols(mc10.dat, 
                      data.frame(offset= seq(-spread,spread,length.out=nrow(mc10.dat))))

# Make plots
BREAKS <- c(0.0,0.01,0.05,0.10,0.20,0.30,0.40,0.6,0.8)
x.labs <- c("None","Mock")
x.at   <- c(1,2)

skew_plot <- function(dat,
                      BREAKS=BREAKS,x.labs=x.labs,x.at=x.at){
  
  if((dat %>% ungroup() %>% distinct(true.prop) %>% nrow() == 3)){ 
    shape.val = c(21,22,23) 
    col.val = pal_jco(palette = c("default"), alpha = 1)(10)[c(1,4,8)]
    }
  if((dat %>% ungroup() %>% distinct(true.prop) %>% nrow()) == 2){ 
     shape.val = c(21,22)
     col.val = pal_jco(palette = c("default"), alpha = 1)(10)[c(1,4)]
    }
  if((dat %>% ungroup() %>% distinct(true.prop) %>% nrow()) == 1){ 
    shape.val = c(21)
    col.val = pal_jco(palette = c("default"), alpha = 1)(10)[c(1)]
    }
  skew.plot <-  ggplot(dat) +
    geom_errorbar(aes(x=1+offset,
                      ymin=raw.2.5, 
                      ymax=raw.97.5,color=true.cat),width=0,alpha=0.5)   +
    geom_point(aes(x=1+offset,y=raw.mean,shape=true.cat,fill=true.cat,color=true.cat),size=2) +
    # mock with mock communities at multiple PCR
    geom_errorbar(aes(x=2+offset,
                      ymin= mock.2.5, 
                      ymax= mock.97.5,color=true.cat),width=0,alpha=0.5)   +
    geom_point(aes(x=2+offset,mock.mean,shape=true.cat,fill=true.cat,color=true.cat),size=2) +
    
    scale_shape_manual(values = shape.val) +
    scale_fill_manual(values= col.val, "True value") +
    scale_color_manual(values= col.val,"True value") +
    scale_y_continuous("Proportion",
                       trans="sqrt",
                       # trans="log",
                       breaks = BREAKS,limits = c(0,max(dat$raw.97.5,dat$mock.97.5)*1.1),expand=c(NA,3)) +
    geom_hline(aes(yintercept = true.prop,color=true.cat),linetype="dashed") +
    # geom_point(aes(x=0.70,y=true.prop,shape=true.cat,fill=true.cat),size=3) +
    scale_x_continuous(name=NULL,breaks=x.at,labels = x.labs) +
    theme_classic() +
    theme(legend.position = "none")
  
  return(skew.plot)
}

p_mc01 <- skew_plot(dat=mc01.dat,
                    BREAKS=BREAKS,x.labs=x.labs,x.at=x.at) 
p_mc02 <- skew_plot(dat=mc02.dat,
                    BREAKS=BREAKS,x.labs=x.labs,x.at=x.at) 
p_mc03 <- skew_plot(dat=mc03.dat,
                    BREAKS=BREAKS,x.labs=x.labs,x.at=x.at) 
p_mc04 <- skew_plot(dat=mc04.dat,
                    BREAKS=BREAKS,x.labs=x.labs,x.at=x.at) 
p_mc05 <- skew_plot(dat=mc05.dat,
                    BREAKS=BREAKS,x.labs=x.labs,x.at=x.at) 
p_mc06 <- skew_plot(dat=mc06.dat,
                    BREAKS=BREAKS,x.labs=x.labs,x.at=x.at) 
p_mc07 <- skew_plot(dat=mc07.dat,
                    BREAKS=BREAKS,x.labs=x.labs,x.at=x.at) 
p_mc08 <- skew_plot(dat=mc08.dat,
                    BREAKS=BREAKS,x.labs=x.labs,x.at=x.at) 
p_mc09 <- skew_plot(dat=mc09.dat,
                    BREAKS=BREAKS,x.labs=x.labs,x.at=x.at) 
p_mc10 <- skew_plot(dat=mc10.dat,
                    BREAKS=BREAKS,x.labs=x.labs,x.at=x.at) 


grid.arrange(p_mc01 + ggtitle(NULL,subtitle="MC01"),
             p_mc02 + ggtitle(NULL,subtitle="MC02"),
             p_mc03 + ggtitle(NULL,subtitle="MC03"),
             p_mc04 + ggtitle(NULL,subtitle="MC04"),
             p_mc05 + ggtitle(NULL,subtitle="MC05"),
             p_mc06 + ggtitle(NULL,subtitle="MC06"),
             p_mc07 + ggtitle(NULL,subtitle="MC07"),
             p_mc08 + ggtitle(NULL,subtitle="MC08"),
             p_mc09 + ggtitle(NULL,subtitle="MC09"),
             p_mc10 + ggtitle(NULL,subtitle="MC10"),
              ncol=2
              )

####################################################33
# Calculate Aitchison Distance
####################################################33

# Aitchison distance
# Calculate for raw, mock, and pcr in 
raw.sp   <- data.frame(Species=raw.post %>% distinct(Species) %>% pull(Species))
raw.comm <- data.frame(community=raw.post %>% distinct(community) %>% pull(community))

mock1.sp   <- data.frame(Species=mock1.post %>% distinct(Species) %>% pull(Species))
mock1.comm <- data.frame(community=mock1.post %>% distinct(community) %>% pull(community))

mock2.sp   <- data.frame(Species=mock2.post %>% distinct(Species) %>% pull(Species))
mock2.comm <- data.frame(community=mock2.post %>% distinct(community) %>% pull(community))

# make true.matrix
true.mat <- result.dat %>% 
              pivot_wider(id_cols=c("Species","community"), values_from = true.prop,names_from = Species) %>%
              arrange(community) %>% as.data.frame()
rownames(true.mat) <- true.mat$community
true.mat <-  true.mat %>% ungroup() %>% dplyr::select(-community)

#####################################################################
# loop over raw
raw.AD <- data.frame(matrix(0,dim(raw$pars$int_samp_small)[1],nrow(true.mat)))
colnames(raw.AD) <- rownames(true.mat)
for(i in 1:dim(raw$pars$int_samp_small)[1]){
  raw.1 <- raw$pars$int_samp_small[i,,] %>% as.data.frame()
  rownames(raw.1) <- raw.comm$community 
  colnames(raw.1) <- raw.sp$Species
  raw.2 <- raw.1 %>% 
              filter(rownames(.) %in% result.dat$community)%>%
              dplyr::select(result.dat$Species) 
  
  for(j in 1:nrow(raw.2)){
    these <- which(true.mat[j,]>0)
    true <- true.mat[j,these]
    raw.3 <- raw.2[j,these] / sum(raw.2[j,these])
    
    raw.AD[i,rownames(true.mat)[j]] <-  aDist(true,raw.3)
  }
  #print(i)
}

# loop over mock1
mock1.AD <- data.frame(matrix(0,dim(mock1$pars$int_samp_small)[1],
                              sum(rownames(true.mat) %in% mock1.comm$community)))
colnames(mock1.AD) <- rownames(true.mat)[rownames(true.mat) %in% mock1.comm$community]
true.mat.trim <- true.mat[rownames(true.mat) %in% mock1.comm$community,]
for(i in 1:dim(mock1$pars$int_samp_small)[1]){
  mock.a <- mock1$pars$int_samp_small[i,,] %>% as.data.frame()
  rownames(mock.a) <- mock1.comm$community 
  colnames(mock.a) <- mock1.sp$Species
  mock.b <- mock.a %>% 
    filter(rownames(.) %in% result.dat$community)%>%
    dplyr::select(result.dat$Species) 
  
  for(j in 1:nrow(mock.b)){
    these <- which(true.mat.trim[j,]>0)
    true <- true.mat.trim[j,these]
    mock.c <- mock.b[j,these] / sum(mock.b[j,these])
    
    mock1.AD[i,rownames(true.mat.trim)[j]] <-  aDist(true,mock.c)
  }
  #print(i)
}

# loop over mock2
mock2.AD <- data.frame(matrix(0,dim(mock2$pars$int_samp_small)[1],
                              sum(rownames(true.mat) %in% mock2.comm$community)))
colnames(mock2.AD) <- rownames(true.mat)[rownames(true.mat) %in% mock2.comm$community]
true.mat.trim <- true.mat[rownames(true.mat) %in% mock2.comm$community,]
for(i in 1:dim(mock2$pars$int_samp_small)[1]){
  mock.a <- mock2$pars$int_samp_small[i,,] %>% as.data.frame()
  rownames(mock.a) <- mock2.comm$community 
  colnames(mock.a) <- mock2.sp$Species
  mock.b <- mock.a %>% 
    filter(rownames(.) %in% result.dat$community)%>%
    dplyr::select(result.dat$Species) 
  
  for(j in 1:nrow(mock.b)){
    these <- which(true.mat.trim[j,]>0)
    true <- true.mat.trim[j,these]
    mock.c <- mock.b[j,these] / sum(mock.b[j,these])
    
    mock2.AD[i,rownames(true.mat.trim)[j]] <-  aDist(true,mock.c)
  }
  #print(i)
}

raw.AD$model   <- "raw"
mock1.AD$model <- "mock"
mock2.AD$model <- "mock"

raw.AD.long <-  raw.AD %>% 
  pivot_longer(.,-model,
               names_to="community",
               values_to = "val")
mock1.AD.long <-  mock1.AD %>% 
  pivot_longer(.,-model,
               names_to="community",
               values_to = "val")
mock2.AD.long <-  mock2.AD %>% 
  pivot_longer(.,-model,
               names_to="community",
               values_to = "val")

all.AD <- bind_rows(raw.AD.long,mock1.AD.long) %>% bind_rows(.,mock2.AD.long)

all.AD.sum <- all.AD %>% group_by(model,community) %>% 
                summarise(Mean = mean(val),
                          Median = median(val),
                          SD=sd(val),
                          q.025 = quantile(val,probs=0.025),
                          q.05 = quantile(val,probs=0.05),
                          q.25 = quantile(val,probs=0.25),
                          q.75 = quantile(val,probs=0.75),
                          q.95 = quantile(val,probs=0.95),
                          q.975 = quantile(val,probs=0.975))  %>%
                ungroup() %>%
                mutate(offset.plot = 0,
                       offset.plot = ifelse(model=="raw",-0.1,offset.plot),
                       offset.plot = ifelse(model=="mock",0.1,offset.plot)) %>%
                # mutate(comm.name = case_when(community=="North_Skew_1_39"~"North\n Skew 1",
                #                              community=="North_Skew_2_39"~"North\n Skew 2",
                #                              community=="Skew_Oceanic_1_39"~"Ocean\n Skew 1",
                #                              community=="Skew_Oceanic_2_39"~"Ocean\n Skew 2",)) %>%
                mutate(Calibration = case_when(
                                      model=="raw"~"None",
                                      model=="mock"~"Mock")) %>%
                as.data.frame()

comm.off <-all.AD %>% distinct(community) %>% 
                mutate(comm.off = seq(-0.2,0.2,length.out=nrow(.)))

all.AD.sum$model <- factor(all.AD.sum$model,levels=c("raw","mock"))                
all.AD.sum <- all.AD.sum %>% left_join(.,comm.off) %>% 
                mutate(model.off = as.numeric(model))

all.AD.sum$comm.idx = as.numeric(as.factor(all.AD.sum$community))
                
xBREAKS <- all.AD.sum %>% distinct(comm.idx,community)
yBREAKS <- c(0,5,10,20,30)

p_AD_seven <- ggplot(all.AD.sum %>% filter(grepl("07",community))) +
          geom_errorbar(aes(x=+ offset.plot,ymin=q.025,ymax=q.975),width=0) +
          geom_errorbar(aes(x=comm.idx + offset.plot,ymin=q.25,ymax=q.75),size=2,width=0) +        
          geom_point(aes(x=comm.idx + offset.plot,y=Mean,shape=Calibration),fill="white",size=2)+
          scale_y_continuous("Aitchison Distance",expand=c(0,0),limits=c(0,NA),breaks=yBREAKS) +
          scale_shape_manual(values=c(21,22,24)) +
          scale_x_continuous(NULL,breaks=xBREAKS$comm.idx,labels=xBREAKS$comm.name,limits=c(NA,NA)) +
          theme_classic() +
          theme(legend.position = "none",
                plot.margin = margin(0,0,0,0.85,"lines"),
                legend.key.size=unit(0.1,'lines'),
                legend.text=element_text(size=7),
                legend.title=element_text(size=9))
p_AD_seven


xBREAKS <- all.AD.sum %>% distinct(model.off,Calibration)
yBREAKS <- c(0,1,2,3,4,5,6)
p_AD_rest <- ggplot(all.AD.sum %>% filter(!grepl("07",community))) +
  geom_errorbar(aes(x=comm.off + model.off,ymin=q.025,ymax=q.975,color=community),width=0) +
  # geom_errorbar(aes(x=comm.off + model.off,ymin=q.25,ymax=q.75,color=community),size=1.5,width=0) +        
  geom_point(aes(x=comm.off + model.off,y=Mean,color=community,fill=community),size=2)+
  scale_y_continuous("Aitchison Distance",expand=c(0,NA),limits=c(0,NA),breaks=yBREAKS) +
  scale_fill_jco(palette = c("default"),"Community",alpha=1) +
  scale_color_jco(palette = c("default"),"Community",alpha=1) +
#scale_shape_manual(values=c(21,24)) +
  scale_x_continuous(NULL,breaks=xBREAKS$model.off,labels=xBREAKS$Calibration,limits=c(0.7,2.7)) +
  theme_classic() +
  theme(legend.position = c(0.88,0.7),
        plot.margin = margin(0,0,0,1.2,"lines"),
        legend.key.size=unit(0.1,'lines'),
        legend.text=element_text(size=7),
        legend.title=element_text(size=9))

p_AD_rest

xBREAKS <- all.AD.sum %>% distinct(model.off,Calibration)
yBREAKS <- c(0,1,3,5,10,15,20,25,30,35,40)
p_AD_all <- ggplot(all.AD.sum) +
  geom_errorbar(aes(x=comm.off + model.off,ymin=q.025,ymax=q.975,color=community),width=0) +
  # geom_errorbar(aes(x=comm.off + model.off,ymin=q.25,ymax=q.75,color=community),size=1.5,width=0) +        
  geom_point(aes(x=comm.off + model.off,y=Mean,color=community,fill=community),size=2)+
  scale_y_continuous("Aitchison Distance",expand=c(0,NA),limits=c(0,NA),breaks=yBREAKS) +
  scale_fill_jco(palette = c("default"),"Community",alpha=1) +
  scale_color_jco(palette = c("default"),"Community",alpha=1) +
  #scale_shape_manual(values=c(21,24)) +
  scale_x_continuous(NULL,breaks=xBREAKS$model.off,labels=xBREAKS$Calibration) +
  theme_classic() +
  theme(legend.position = c(0.5,0.7),
        #plot.margin = margin(0,0,0,0.85,"lines"),
        legend.key.size=unit(0.1,'lines'),
        legend.text=element_text(size=7),
        legend.title=element_text(size=9))

p_AD_all



###############################################################333
#### Pull out estimates of alpha, convert to CLR
###############################################################333

p_space_mock1 <- (exp(mock1$pars$alpha) / rowSums(exp(mock1$pars$alpha))) %>% as.data.frame()
clr_alpha_list_mock1 <- cenLR(p_space_mock1)
clr_alpha <- clr_alpha_list_mock1$x.clr %>% as.data.frame()
colnames(clr_alpha) <- mock.sp$Species                       
clr_alpha_sum_mock1 <- clr_alpha %>% 
                    pivot_longer( .,
                          cols = colnames(clr_alpha),
                          names_to="species",values_to="val") %>%
                    group_by(species) %>%
                    summarise(Mean = mean(val),
                        SD=sd(val),
                        q.025 = quantile(val,probs=0.025),
                        q.05 = quantile(val,probs=0.05),
                        q.25 = quantile(val,probs=0.25),
                        q.75 = quantile(val,probs=0.75),
                        q.95 = quantile(val,probs=0.95),
                        q.975 = quantile(val,probs=0.975))    

clr_alpha_sum_mock1 <- clr_alpha_sum_mock1 %>% arrange(Mean) %>% mutate(model="Odd")

p_space_mock2 <- (exp(mock2$pars$alpha) / rowSums(exp(mock2$pars$alpha))) %>% as.data.frame()
clr_alpha_list_mock2 <- cenLR(p_space_mock2)
clr_alpha <- clr_alpha_list_mock2$x.clr %>% as.data.frame()
colnames(clr_alpha) <- mock.sp$Species                       
clr_alpha_sum_mock2 <- clr_alpha %>% 
  pivot_longer( .,
                cols = colnames(clr_alpha),
                names_to="species",values_to="val") %>%
  group_by(species) %>%
  summarise(Mean = mean(val),
            SD=sd(val),
            q.025 = quantile(val,probs=0.025),
            q.05 = quantile(val,probs=0.05),
            q.25 = quantile(val,probs=0.25),
            q.75 = quantile(val,probs=0.75),
            q.95 = quantile(val,probs=0.95),
            q.975 = quantile(val,probs=0.975))    

clr_alpha_sum_mock2 <- clr_alpha_sum_mock2 %>% arrange(Mean) %>% mutate(model="Even")

clr_alpha_sum <- bind_rows(clr_alpha_sum_mock1,clr_alpha_sum_mock2)

# get rid of reference species denotion
clr_alpha_sum <- clr_alpha_sum %>%
                      mutate(SP= ifelse(grepl("zRefSpecies_",species),
                                       substr(species,13,nchar(species)),
                                       as.character(species)))

clr_alpha_sum$SP <- factor(clr_alpha_sum$SP,
                                levels = clr_alpha_sum %>% distinct(SP) %>% pull(SP))

col.val = pal_jco(palette = c("default"), alpha = 1)(10)[c(6,7)]

p_clr_mock <- ggplot(clr_alpha_sum) +
    #geom_errorbarh(aes(xmin=q.25,xmax=q.75,y=SP,color=model),size=2,height=0) +
    geom_errorbarh(aes(xmin=q.025,xmax=q.975,y=SP,color=model),size=0.8,height=0,alpha=0.8) +
    geom_point(aes(x=Mean,y=SP,fill=model,color=model),size=2,shape=21,alpha=0.8) +
    geom_vline(xintercept=0,linetype="dashed") +
    scale_x_continuous("Amplification Efficiency (CLR)") +
    scale_y_discrete(NULL) +
    scale_color_manual(values=col.val,"Calibration\nCommunity") +
    scale_fill_manual(values=col.val,"Calibration\nCommunity") +
    theme_bw() +
    theme(legend.position = c(0.3,0.7),
          legend.key.size=unit(0.1,'lines'),
          legend.text=element_text(size=7),
          legend.title=element_text(size=9),
          axis.text.y = element_text(size=7))


###############################################################
#### Just use ALR estimates to make them directly comparable between odd and even communities.
###############################################################

alr_mock1 <- mock1$stanMod_summary$alpha %>%  as.data.frame() %>% # Odd
                mutate(model="Odd",
                       sp_idx = 1:nrow(.))%>% 
                left_join(.,mock1$mc %>% distinct(Species,sp_idx)) %>%
                            arrange(mean) %>%
                            as.data.frame()
                
alr_mock2 <- mock2$stanMod_summary$alpha %>%  as.data.frame() %>% # Odd
                mutate(model="Even",
                       sp_idx = 1:nrow(.))%>% 
                left_join(.,mock2$mc %>% distinct(Species,sp_idx)) %>%
                            arrange(mean) %>%
                            as.data.frame()

alr_all <- bind_rows(alr_mock1,alr_mock2)
colnames(alr_all)[4:8] <- c("q.025","q.25","median","q.75","q.975")

alr_all <- alr_all %>% mutate(SP= ifelse(grepl("zRefSpecies_",Species),
                            substr(Species,13,nchar(Species)),
                            as.character(Species))) %>% 
                       mutate(SP = as.factor(SP))

alr_all$SP <- factor(alr_all$SP,levels=c(alr_all %>% distinct(SP) %>% unlist() %>% as.character()))  


col.val = pal_jco(palette = c("default"), alpha = 1)(10)[c(6,7)]

p_alr_mock <- ggplot(alr_all) +
  #geom_errorbarh(aes(xmin=q.25,xmax=q.75,y=SP,color=model),size=2,height=0) +
  geom_errorbarh(aes(xmin=q.025,xmax=q.975,y=SP,color=model),size=0.8,height=0,alpha=0.8) +
  geom_point(aes(x=mean,y=SP,fill=model,color=model),size=2,shape=21,alpha=0.8) +
  geom_vline(xintercept=0,linetype="dashed") +
  scale_x_continuous("Amplification Efficiency (ALR)") +
  scale_y_discrete(NULL) +
  scale_color_manual(values=col.val,"Calibration\nCommunity") +
  scale_fill_manual(values=col.val,"Calibration\nCommunity") +
  theme_bw() +
  theme(legend.position = c(0.3,0.7),
        legend.key.size=unit(0.1,'lines'),
        legend.text=element_text(size=7),
        legend.title=element_text(size=9),
        axis.text.y = element_text(size=7))







##################
#################3



LAY <- rbind(c(1,4),c(2,4),c(3,4))

quartz(file="./plots/Hanfling_Figure_CLR.pdf",height=6,width=7,dpi=600,type="pdf")
  grid.arrange(p_mc03+
               annotate(geom="text",x=0.7,y=0.28,label="A"),
             p_mc08+
               annotate(geom="text",x=0.7,y=0.32,label="B"),
             p_AD_rest+
               annotate(geom="text",x=0.7,y=6,label="C"),
             p_clr_mock +
               annotate(geom="text",x=-0.3,y=21,label="D"),
             layout_matrix=LAY
  )

dev.off()


LAY <- rbind(c(1,4),c(2,4),c(3,4))

quartz(file="./plots/Hanfling_Figure_ALR.pdf",height=6,width=7,dpi=600,type="pdf")
grid.arrange(p_mc03+
               annotate(geom="text",x=0.7,y=0.28,label="A"),
             p_mc08+
               annotate(geom="text",x=0.7,y=0.32,label="B"),
             p_AD_rest+
               annotate(geom="text",x=0.7,y=6,label="C"),
             p_alr_mock +
               annotate(geom="text",x=-0.35,y=21,label="D"),
             layout_matrix=LAY
)

dev.off()







quartz(file="./plots/Hanfling_Figure_all_comm.pdf",height=10,width=6,dpi=600,type="pdf")

grid.arrange(p_mc01 + ggtitle(NULL,subtitle="MC01"),
             p_mc02 + ggtitle(NULL,subtitle="MC02"),
             p_mc03 + ggtitle(NULL,subtitle="MC03"),
             p_mc04 + ggtitle(NULL,subtitle="MC04"),
             p_mc05 + ggtitle(NULL,subtitle="MC05"),
             p_mc06 + ggtitle(NULL,subtitle="MC06"),
             p_mc07 + ggtitle(NULL,subtitle="MC07"),
             p_mc08 + ggtitle(NULL,subtitle="MC08"),
             p_mc09 + ggtitle(NULL,subtitle="MC09"),
             p_mc10 + ggtitle(NULL,subtitle="MC10"),
             ncol=2
)
dev.off()



quartz(file="./plots/Hanfling_Figure_AD_all.pdf",height=3,width=4,dpi=600,type="pdf")
  print(p_AD_all)
dev.off()
