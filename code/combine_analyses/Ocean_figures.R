# Compile fits for Ocean fish data and make figures of Ocean data for the publication.
library(tidyverse)
library(data.table)
library(gridExtra)
library(ggsci)
library(robCompositions)

# Read in the posteriors for the 
load("./data/summarized_data/Fish_Oceanic_and_North_Raw.Rdata")
raw <- Output
load("./data/summarized_data/Fish_Oceanic_and_North_mock_2_crossValidate_Skew.Rdata")
mock2 <- Output # This is the cross validation that uses only the 39 cycle even communities

# raw - no amplification adjustment.
# mock2 - many tau, single 39 cycle even communities

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
colnames(beta_posterior) <- paste0("mock2.",substr(colnames(beta_posterior),1,nchar(colnames(beta_posterior))-1))
colnames(beta_posterior)[1] <- "mock2.mean"
beta_posterior <- as.data.frame(beta_posterior)

mock2.post <-expand.grid(comm_idx = COM$comm_idx,sp_idx =SP$sp_idx) %>% 
  arrange(comm_idx,sp_idx) %>% 
  left_join(.,COM) %>% 
  left_join(.,SP) %>% 
  bind_cols(.,beta_posterior)

# Combine the raw estimates and posterior estimates
mock2.all <- full_join(mock2.raw,mock2.post)


########################################################################
########################################################################
########################################################################
########################################################################

# Combine mock and pcr_var results with raw reads.
result.dat <- left_join(mock2.all,raw.all %>% dplyr::select(-comm_idx,-sp_idx))                        

# make a distinct factor for large true proportions and for small true proportions
result.dat <- result.dat %>% mutate(true.cat= ifelse(true.prop>0.1,"big","small"))

# add manual color column for later.
result.dat <- result.dat %>% mutate(manual.col = "white",
                                    manual.col = ifelse(grepl("Hippo",Species),
                                                        pal_jco(palette = c("default"), alpha = 1)(10)[c(2)],
                                                          manual.col),
                                    manual.col = ifelse(grepl("Squalus",Species),
                                                        pal_jco(palette = c("default"), alpha = 1)(10)[c(9)],
                                                        manual.col),
                                    manual.col = ifelse(grepl("Leurog",Species),
                                                        pal_jco(palette = c("default"), alpha = 1)(10)[c(8)],
                                                        manual.col)
                                    )

# pull out just ocean.skew.dat for plotting.
spread=0.15
ocean.skew1.dat <- result.dat %>% 
  filter(true.prop > 0,community=="Skew_Oceanic_1_39")
ocean.skew1.dat <- bind_cols(ocean.skew1.dat, 
                             data.frame(offset= seq(-spread,spread,length.out=nrow(ocean.skew1.dat))))

ocean.skew2.dat <- result.dat %>% 
  filter(true.prop > 0,community=="Skew_Oceanic_2_39")
ocean.skew2.dat <- bind_cols(ocean.skew2.dat, 
                            data.frame(offset= seq(-spread,spread,length.out=nrow(ocean.skew2.dat))))

# pull out just north.skew.dat for plotting.
north.skew1.dat <- result.dat %>% 
  filter(true.prop > 0,community=="North_Skew_1_39")
north.skew1.dat <- bind_cols(north.skew1.dat, 
                             data.frame(offset= seq(-spread,spread,length.out=nrow(north.skew1.dat))))

north.skew2.dat <- result.dat %>% 
  filter(true.prop > 0,community=="North_Skew_2_39")
north.skew2.dat <- bind_cols(north.skew2.dat, 
                             data.frame(offset= seq(-spread,spread,length.out=nrow(north.skew1.dat))))

# Make plots
BREAKS <- c(0.0,0.01,0.05,0.10,0.20,0.30,0.40,0.6,0.8)
x.labs <- c("None","Mock","Mock B")
x.at   <- c(1,2,3)

skew_plot <- function(dat,
                      BREAKS=BREAKS,x.labs=x.labs,x.at=x.at){
  
  shape.val = c(21,22)
  col.val = pal_jco(palette = c("default"), alpha = 1)(10)[c(1,4)]
  
  skew.plot <-  ggplot(dat) +
    geom_errorbar(aes(x=1+offset,
                      ymin=raw.2.5, 
                      ymax=raw.97.5,color=true.cat),width=0,alpha=0.75)   +
    geom_point(aes(x=1+offset,y=raw.mean,shape=true.cat,fill=true.cat,color=true.cat),size=2) +
    # mock with mock communities at multiple PCR
    # geom_errorbar(aes(x=3+offset,
    #                   ymin= mock1.2.5, 
    #                   ymax= mock1.97.5),width=0,alpha=0.5)   +
    # geom_point(aes(x=3+offset,mock1.mean,shape=true.cat,fill=true.cat,color=true.cat)) +
    # mock with mock communities at single PCR
    geom_errorbar(aes(x=2+offset,
                      ymin= mock2.2.5, 
                      ymax= mock2.97.5,color=true.cat),width=0,alpha=0.75)   +
    geom_point(aes(x=2+offset,mock2.mean,shape=true.cat,fill=true.cat,color=true.cat),size=2) +
    # geom_point(aes(x=4+offset,pcr.mean,shape=true.cat,fill=true.cat)) +
    # geom_errorbar(aes(x=4+offset,
    #                   ymin= pcr.2.5, 
    #                   ymax= pcr.97.5),width=0)   +
    scale_shape_manual(values =c(21,22)) +
    scale_fill_manual(values= col.val, "True value") +
    scale_color_manual(values= col.val,"True value") +
    scale_y_continuous("Proportion",
                       trans="sqrt",
                       # trans="log",
                       breaks = BREAKS,expand=c(0,NA),limits = c(0,NA)) +
    geom_hline(aes(yintercept = true.prop,color=true.cat),linetype="dashed") +
    #geom_point(aes(x=0.70,y=true.prop,shape=true.cat,fill=true.cat),size=3) +
    scale_x_continuous(name=NULL,breaks=x.at,labels = x.labs) +
    theme_classic() +
    theme(legend.position = "none")
  
  return(skew.plot)
}


ocean.skew.1 <- skew_plot(dat=ocean.skew1.dat,
                          BREAKS=BREAKS,x.labs=x.labs,x.at=x.at) 
                

ocean.skew.2 <- skew_plot(dat=ocean.skew2.dat,
                          BREAKS=BREAKS,x.labs=x.labs,x.at=x.at) 

north.skew.1 <- skew_plot(dat=north.skew1.dat,
                          BREAKS=BREAKS,x.labs=x.labs,x.at=x.at) 

north.skew.2 <- skew_plot(dat=north.skew2.dat,
                          BREAKS=BREAKS,x.labs=x.labs,x.at=x.at) 


grid.arrange(ocean.skew.1 + ggtitle(NULL,subtitle="Ocean Skew 1"),
              ocean.skew.2 +ggtitle(NULL,subtitle="Ocean Skew 2"),
              north.skew.1 + ggtitle(NULL,subtitle="North Skew 1"),
              north.skew.2 + ggtitle(NULL,subtitle="North Skew 2"),
              ncol=2,nrow=2
              )

####################################################33
# Calculate Aitchison Distance
####################################################33


# Aitchison distance
# Calculate for raw, mock, and pcr in 

# Make containers for each type of output
make.cont <- function(dat){
  out <-  matrix(0,max(dat$comm_idx),max(dat$sp_idx))
  rownames(out) <- dat %>% distinct(community) %>% pull(community)
  colnames(out) <- dat %>% distinct(Species) %>% pull(Species)
  return(out)
}

raw.out <- make.cont(raw.post)
mock1.out <- make.cont(mock1.post)
mock2.out <- make.cont(mock2.post)
pcr.out <- make.cont(pcr.post)

raw.sp   <- data.frame(Species=raw.post %>% distinct(Species) %>% pull(Species))
raw.comm <- data.frame(community=raw.post %>% distinct(community) %>% pull(community))

mock.sp   <- data.frame(Species=mock1.post %>% distinct(Species) %>% pull(Species))
mock.comm <- data.frame(community=mock1.post %>% distinct(community) %>% pull(community))

pcr.sp   <- data.frame(Species=pcr.post %>% distinct(Species) %>% pull(Species))
pcr.comm <- data.frame(community=pcr.post %>% distinct(community) %>% pull(community))

# make true.matrix
true.mat <- result.dat %>% 
              pivot_wider(id_cols=c("Species","community"), values_from = true.prop,names_from = Species) %>%
              as.data.frame()
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
mock1.AD <- data.frame(matrix(0,dim(mock1$pars$int_samp_small)[1],nrow(true.mat)))
colnames(mock1.AD) <- rownames(true.mat)
for(i in 1:dim(mock1$pars$int_samp_small)[1]){
  mock.a <- mock1$pars$int_samp_small[i,,] %>% as.data.frame()
  rownames(mock.a) <- mock.comm$community 
  colnames(mock.a) <- mock.sp$Species
  mock.b <- mock.a %>% 
    filter(rownames(.) %in% result.dat$community)%>%
    dplyr::select(result.dat$Species) 
  
  for(j in 1:nrow(mock.b)){
    these <- which(true.mat[j,]>0)
    true <- true.mat[j,these]
    mock.c <- mock.b[j,these] / sum(mock.b[j,these])
    
    mock1.AD[i,rownames(true.mat)[j]] <-  aDist(true,mock.c)
  }
  #print(i)
}

# loop over mock1
mock2.AD <- data.frame(matrix(0,dim(mock2$pars$int_samp_small)[1],nrow(true.mat)))
colnames(mock2.AD) <- rownames(true.mat)
for(i in 1:dim(mock2$pars$int_samp_small)[1]){
  mock.a <- mock2$pars$int_samp_small[i,,] %>% as.data.frame()
  rownames(mock.a) <- mock.comm$community 
  colnames(mock.a) <- mock.sp$Species
  mock.b <- mock.a %>% 
    filter(rownames(.) %in% result.dat$community)%>%
    dplyr::select(result.dat$Species) 
  
  for(j in 1:nrow(mock.b)){
    these <- which(true.mat[j,]>0)
    true <- true.mat[j,these]
    mock.c <- mock.b[j,these] / sum(mock.b[j,these])
    
    mock2.AD[i,rownames(true.mat)[j]] <-  aDist(true,mock.c)
  }
  #print(i)
}

# loop over pcr
pcr.AD <- data.frame(matrix(0,dim(var_pcr$pars$int_samp_small)[1],nrow(true.mat)))
colnames(pcr.AD) <- rownames(true.mat)
for(i in 1:dim(var_pcr$pars$int_samp_small)[1]){
  pcr.a <- var_pcr$pars$int_samp_small[i,,] %>% as.data.frame()
  rownames(pcr.a) <- pcr.comm$community 
  colnames(pcr.a) <- pcr.sp$Species
  pcr.b <- pcr.a %>% 
    filter(rownames(.) %in% result.dat$community)%>%
    dplyr::select(result.dat$Species) 
  
  for(j in 1:nrow(pcr.b)){
    these <- which(true.mat[j,]>0)
    true <- true.mat[j,these]
    pcr.c <- pcr.b[j,these] / sum(pcr.b[j,these])
    
    pcr.AD[i,rownames(true.mat)[j]] <-  aDist(true,pcr.c)
  }
  print(i)
}

raw.AD$model <- "raw"
mock1.AD$model <- "mock1"
mock2.AD$model <- "mock2"
pcr.AD$model <- "pcr"

all.AD <- bind_rows(raw.AD,mock1.AD) %>% bind_rows(.,mock2.AD) %>% bind_rows(.,pcr.AD)

all.AD <- all.AD %>% 
              pivot_longer(.,-model,
                           names_to="community",
                           values_to = "val")
all.AD.sum <- all.AD %>% group_by(model,community) %>% 
                summarise(Mean = mean(val),
                          SD=sd(val),
                          q.025 = quantile(val,probs=0.025),
                          q.05 = quantile(val,probs=0.05),
                          q.25 = quantile(val,probs=0.25),
                          q.75 = quantile(val,probs=0.75),
                          q.95 = quantile(val,probs=0.95),
                          q.975 = quantile(val,probs=0.975))  %>%
                ungroup() %>%
                mutate(offset.plot = 0,
                       offset.plot = ifelse(model=="raw",-0.05,offset.plot),
                       offset.plot = ifelse(model=="mock2",0.05,offset.plot), 
                       offset.plot = ifelse(model=="pcr",0.2,offset.plot)) %>%
                mutate(comm.name = case_when(community=="North_Skew_1_39"~"North\n Skew 1",
                                             community=="North_Skew_2_39"~"North\n Skew 2",
                                             community=="Skew_Oceanic_1_39"~"Ocean\n Skew 1",
                                             community=="Skew_Oceanic_2_39"~"Ocean\n Skew 2",)) %>%
                mutate(Calibration = case_when(model=="raw"~"None",
                             model=="mock1"~"Mock B",
                             model=="mock2"~"Mock",
                             model=="pcr"~"Variable PCR",))
                as.data.frame()

all.AD.sum$comm.idx = as.numeric(as.factor(all.AD.sum$community))
                


xBREAKS <- all.AD.sum %>% distinct(comm.idx,comm.name)
yBREAKS <- c(0,5,10,20,30)

col.val <- pal_jco(palette = c("default"), alpha = 1)(10)[c(6,7)]

p_AD_north <- ggplot(all.AD.sum %>% filter(!model %in% c("pcr","mock1"),grepl("North",community))) +
          geom_errorbar(aes(x=comm.idx + offset.plot,ymin=q.025,ymax=q.975,color=model),width=0) +
          geom_errorbar(aes(x=comm.idx + offset.plot,ymin=q.25,ymax=q.75,color=model),size=2,width=0) +        
          geom_point(aes(x=comm.idx + offset.plot,y=Mean,color=model ),shape =21,fill="white",size=2) +
          scale_y_continuous(NULL,expand=c(0,NA),limits=c(0,NA),breaks=yBREAKS) +
          #scale_shape_manual(values=c(21,22,24)) +
          #scale_fill_manual(values=col.val,"Community") +
          scale_color_manual(values=col.val,"Calibration") +
          scale_x_continuous(NULL,breaks=xBREAKS$comm.idx,labels=xBREAKS$comm.name,limits=c(NA,2.3)) +
          theme_classic() +
          theme(legend.position = "none",
                plot.margin = margin(0,0,0,0.85,"lines"),
                legend.key.size=unit(0.1,'lines'),
                legend.text=element_text(size=7),
                legend.title=element_text(size=9))
          
p_AD_north


yBREAKS <- c(0,1,2,3,4,5)
p_AD_ocean <- ggplot(all.AD.sum %>% filter(!model %in% c("pcr","mock1"),grepl("Ocean",community))) +
  geom_errorbar(aes(x=comm.idx + offset.plot,ymin=q.025,ymax=q.975,color=Calibration),width=0) +
  geom_errorbar(aes(x=comm.idx + offset.plot,ymin=q.25,ymax=q.75,color=Calibration),size=2,width=0) +        
  geom_point(aes(x=comm.idx + offset.plot,y=Mean,color=Calibration),shape=21,fill="white",size=2)+
  scale_y_continuous("Aitchison Distance",expand=c(0,0),limits=c(0,NA),breaks=yBREAKS) +
  #scale_shape_manual(values=c(21,22,24)) +
  scale_color_manual(values=col.val,"Calibration") +
  scale_x_continuous(NULL,breaks=xBREAKS$comm.idx,labels=xBREAKS$comm.name,limits=c(NA,4.3)) +
  theme_classic() +
  theme(legend.position = c(0.4755,0.85),
        plot.margin = margin(0,0,0,0.85,"lines"),
        legend.key.size=unit(0.1,'lines'),
        legend.text=element_text(size=7),
        legend.title=element_text(size=9))

p_AD_ocean


###############################################################333
#### Pull out estimates of alpha, convert to CLR
###############################################################333

p_space_mock1 <- (exp(mock1$pars$alpha) / rowSums(exp(mock1$pars$alpha))) %>% as.data.frame()
clr_alpha_list_mock1 <- cenLR(p_space_mock1)
clr_alpha <- clr_alpha_list_mock1$x.clr %>% as.data.frame()
colnames(clr_alpha) <- mock.sp$Species                       
clr_alpha_sum <- clr_alpha %>% 
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

clr_alpha_sum <- clr_alpha_sum %>% arrange(Mean)

# get rid of reference species denotion
clr_alpha_sum <- clr_alpha_sum %>%
                      mutate(SP= ifelse(grepl("zRefSpecies_",species),
                                       substr(species,13,nchar(species)),
                                       as.character(species)))
clr_alpha_sum$SP <- factor(clr_alpha_sum$SP,
                                levels = clr_alpha_sum$SP)

# Add manual colors
clr_alpha_sum <- clr_alpha_sum %>% mutate(manual.col = "white",
           manual.col = ifelse(grepl("Hippo",species),
                               pal_jco(palette = c("default"), alpha = 1)(10)[c(2)],
                               manual.col),
           manual.col = ifelse(grepl("Squalus",species),
                               "#00A087FF",
                               manual.col),
           manual.col = ifelse(grepl("Leurog",species),
                               pal_jco(palette = c("default"), alpha = 1)(10)[c(8)],
                               manual.col)
)



p_clr_mock1 <- ggplot(clr_alpha_sum) +
    geom_errorbarh(aes(xmin=q.25,xmax=q.75,y=SP),size=2,height=0) +
    geom_errorbarh(aes(xmin=q.025,xmax=q.975,y=SP),size=0.8,height=0) +
    geom_point(aes(x=Mean,y=SP,fill=SP,),size=3,shape=21) +
    geom_vline(xintercept=0,linetype="dashed") +
    scale_fill_manual(values=clr_alpha_sum$manual.col %>% as.character()) +
      scale_x_continuous("Amplification Efficiency (CLR)") +
    scale_y_discrete(NULL) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.y = element_text(size=7))


#clr_alpha_sum$SP2 <- clr_alpha_sum$SP


##################
#################3

LAY <- rbind(c(1,2,3),c(4,5,3))

grid.arrange(ocean.skew.1,
             ocean.skew.2,
             p_clr_mock1,
             north.skew.1,
             north.skew.2,
             layout_matrix=LAY
)



LAY <- rbind(c(1,2),c(3,4))

quartz(file="./plots/Ocean_Fish_Figure.pdf",height=6,width=7,dpi=600,type="pdf")
  grid.arrange(ocean.skew.1 +
               annotate(geom="text",x=0.7,y=0.38,label="A"),
             ocean.skew.2+
               annotate(geom="text",x=0.7,y=0.36,label="B"),
             north.skew.1 +
               annotate(geom="text",x=0.7,y=0.60,label="C"),
             north.skew.2+
               annotate(geom="text",x=0.7,y=0.80,label="D"),
             layout_matrix=LAY
  )

dev.off()


LAY <- rbind(c(1,4),c(2,4),c(3,4))
quartz(file="./plots/Ocean_Fish_Figure2.pdf",height=6,width=8,dpi=600,type="pdf")
grid.arrange(ocean.skew.1 +
               annotate(geom="text",x=0.7,y=0.38,label="A"),
             north.skew.2+
               annotate(geom="text",x=0.7,y=0.85,label="B"),
             p_AD +
               annotate(geom="text",x=0.7,y=35,label="C"),
             p_clr_mock1 +
               annotate(geom="text",x=-0.14,y=34,label="D"),
             layout_matrix=LAY
)

dev.off()


quartz(file="./plots/Ocean_Fish_North_Ocean.pdf",height=6,width=8,dpi=600,type="pdf")
grid.arrange(ocean.skew.1 + ggtitle(NULL,subtitle="Ocean Skew 1"),
             ocean.skew.2 +ggtitle(NULL,subtitle="Ocean Skew 2"),
             north.skew.1 + ggtitle(NULL,subtitle="North Skew 1"),
             north.skew.2 + ggtitle(NULL,subtitle="North Skew 2"),
             ncol=2,nrow=2
)
dev.off()


LAY <- rbind(c(1,1,5),c(2,2,5),c(3,4,5))
quartz(file="./plots/Ocean_Fish_Figure3.pdf",height=6,width=8,dpi=600,type="pdf")
grid.arrange(ocean.skew.1 +
               geom_point(data=ocean.skew1.dat %>% filter(grepl("Leuro",Species)),
                            aes(x=1+offset,y=raw.mean,shape=true.cat),
                          fill="#3B3B3BFF",color="#3B3B3BFF") +
               geom_point(data=ocean.skew1.dat %>% filter(grepl("Leuro",Species)),
                          aes(x=2+offset,y=mock2.mean,shape=true.cat),
                          fill="#3B3B3BFF",color="#3B3B3BFF") +
               annotate(geom="text",x=0.7,y=0.38,label="A"),
             north.skew.2+
               geom_point(data=north.skew2.dat %>% filter(grepl("Hippo",Species) | grepl("Squalus",Species)),
                          aes(x=1+offset,y=raw.mean,shape=true.cat),
                          fill=c("#EFC000FF","#00A087FF"),color=c("#EFC000FF","#00A087FF")) +
               geom_point(data=north.skew2.dat %>% filter(grepl("Hippo",Species) | grepl("Squalus",Species)),
                          aes(x=2+offset,y=mock2.mean,shape=true.cat),
                         fill=c("#EFC000FF","#00A087FF"),color=c("#EFC000FF","#00A087FF")) +
               annotate(geom="text",x=0.7,y=0.85,label="B"),
             p_AD_ocean +
               annotate(geom="text",x=2.7,y=3.15,label="C"),
             p_AD_north +
               annotate(geom="text",x=0.7,y=35,label="D"),
             p_clr_mock1 +
               annotate(geom="text",x=-0.14,y=34,label="E"),
             widths=c(0.25,0.25,0.5),
             layout_matrix=LAY
)

dev.off()

##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
## Repeat some of the above analyses but remove the shark, Squalus 
## to illustrate the effect of one terrible amplifier on
## the estimated proportions.
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################


# Here we go get the posteriors for two model fits (raw, mock2)
# and then reproject for all  species except for Squalus.

# go get the posterior estimates for the beta parameters and 
# relevant design matrices
raw.beta <- raw$pars$beta
mock.beta <- mock2$pars$beta

design.raw <- raw$stan_data$model_matrix_b_samp_small
design.mock <- mock2$stan_data$model_matrix_b_samp_small

# write a function to calculate the predictions 
predict.comp <- function(beta.mat, design.mat, sp.list, sp.exclude){
  
      these.sp <- sp.list %>% 
                    mutate(id = 1:nrow(sp.list)) %>% filter(!Species %in% sp.exclude)
      N.sp <- nrow(these.sp)
      N.obs.small <- nrow(design.mat)
      
      logit_val_samp_ppd <- matrix(0,N.obs.small,N.sp) 
      prob_samp_ppd    <- matrix(0,N.obs.small,nrow(these.sp)) 
      
      int_samp_small <- NULL
      for(i in 1:(dim(beta.mat)[1])){ # loop over MCMC iterations.
      
        for(n in 1:N.sp){
           logit_val_samp_ppd[,n] = design.mat %*% beta.mat[i,these.sp$id[n],] ; 
        }
        for(m in 1:N.obs.small){
            prob_samp_ppd[m,] = exp(logit_val_samp_ppd[m,]) / sum(exp(logit_val_samp_ppd[m,]))
        }
        prob_samp_mod <- prob_samp_ppd %>% as.data.frame()
        prob_samp_mod$mcmc = i
        prob_samp_mod$comm <- 1:ncol(design.mat)
        int_samp_small =  bind_rows(int_samp_small,prob_samp_mod)
      } # end loop over MCMC
        
      colnames(int_samp_small)[1:N.sp] <- these.sp$Species
      return(int_samp_small)
      
} # end function   
        

# Evaluate for both raw and mock model fits.
raw.int_samp_small <-  predict.comp(raw.beta, 
                                    design.raw, 
                                    sp.list = data.frame(Species=raw$Species),
                                    sp.exclude = c("Squalus acanthias"))       
      
mock.int_samp_small <-  predict.comp(mock.beta, 
                                    design.mock, 
                                    sp.list = data.frame(Species=mock2$Species),
                                    sp.exclude = c("Squalus acanthias"))       

# pivot longer and summarize

raw.long.int_samp_small <- pivot_longer(raw.int_samp_small,
                                        cols = !c("mcmc","comm"),
                                        names_to = "Species",
                                        values_to = "prop") %>% 
                                  group_by(Species,comm) %>%
                                  summarize(Mean = mean(prop),
                                            q.0.025 = quantile(prop,probs=c(0.025)),
                                            q.0.975 = quantile(prop,probs=c(0.975)))
  
mock.long.int_samp_small <- pivot_longer(mock.int_samp_small,
                                        cols = !c("mcmc","comm"),
                                        names_to = "Species",
                                        values_to = "prop") %>%
                              group_by(Species,comm) %>%
                              summarize(Mean = mean(prop),
                                q.0.025 = quantile(prop,probs=c(0.025)),
                                q.0.975 = quantile(prop,probs=c(0.975)))

# Merge in community names
COM <- data.frame(community = levels(raw$env$community %>% as.factor()))
raw.comm <- data.frame(comm = 1:nrow(COM),COM)

COM <- data.frame(community = levels(mock2$env$community %>% as.factor()))
mock.comm <- data.frame(comm = 1:nrow(COM),COM)

raw.long.int_samp_small <- left_join(raw.long.int_samp_small,raw.comm)
mock.long.int_samp_small <- left_join(mock.long.int_samp_small,mock.comm)

# fix community names and drop even communities

raw.long.int_samp_small <- raw.long.int_samp_small %>% 
                              filter(grepl("Skew",community))

mock.long.int_samp_small <- mock.long.int_samp_small %>% 
                              filter(grepl("Skew",community))

#### OK.  Do calculations based upon base data, excluding c("Squalus acanthias")
raw.raw <- raw$env %>% filter(!Species==c("Squalus acanthias")) %>%
  group_by(community,Cycles,tech_rep) %>%
  mutate(sum.ng = sum(start_conc_ng),
         true.prop = start_conc_ng / sum.ng) %>%
  ungroup() %>%
  group_by(Species,community,Cycles,true.prop) %>%
  summarise(simple.Mean=mean(propReads),
            simple.N = length(tech_rep)) %>%
  replace_na(list(raw.Mean=0,raw.SD=0,raw.SE=0))


mock2.raw <- mock2$env %>% filter(!Species==c("Squalus acanthias")) %>%
  group_by(community,Cycles,tech_rep) %>%
  mutate(sum.ng = sum(start_conc_ng),
         true.prop = start_conc_ng / sum.ng) %>%
  ungroup() %>%
  group_by(Species,community,Cycles,true.prop) %>%
  summarise(simple.Mean=mean(propReads),
            simple.N = length(tech_rep)) %>%
  replace_na(list(raw.Mean=0,raw.SD=0,raw.SE=0))

# filter to include only species that have non-zero true proportions, 
# merge with predictions based on projections without Squalus
# Trim to only include the North Skew (only ones where Squalus was used)

raw.trim <- raw.raw %>% filter(true.prop > 0,grepl("North_Skew",community) ) %>% 
            left_join(.,raw.long.int_samp_small) %>%
            mutate(type="raw")

mock2.trim <- mock2.raw %>% filter(true.prop > 0,grepl("North_Skew",community) ) %>% 
            left_join(.,mock.long.int_samp_small) %>%
            mutate(type="mock")

# Only 


# Combine mock and raw
result.dat <- bind_rows(raw.trim,mock2.trim)

# make a distinct factor for large true proportions and for small true proportions
result.dat <- result.dat %>% mutate(true.cat= ifelse(true.prop>0.1,"big","small"))

# add manual color column for later.
result.dat <- result.dat %>% mutate(manual.col = "white",
                                    manual.col = ifelse(grepl("Hippo",Species),
                                                        pal_jco(palette = c("default"), alpha = 1)(10)[c(2)],
                                                        manual.col),
                                    manual.col = ifelse(grepl("Squalus",Species),
                                                        pal_jco(palette = c("default"), alpha = 1)(10)[c(9)],
                                                        manual.col),
                                    manual.col = ifelse(grepl("Leurog",Species),
                                                        pal_jco(palette = c("default"), alpha = 1)(10)[c(8)],
                                                        manual.col)
)



# pull out just ocean.skew.dat for plotting.
spread=0.15

# pull out just north.skew.dat for plotting.
north.skew1.dat <- result.dat %>% 
  filter(true.prop > 0,community=="North_Skew_1_39")
north.skew1.dat <- bind_cols(north.skew1.dat, 
                             data.frame(offset= seq(-spread,spread,length.out=nrow(north.skew1.dat))))

north.skew2.dat <- result.dat %>% 
  filter(true.prop > 0,community=="North_Skew_2_39")
north.skew2.dat <- bind_cols(north.skew2.dat, 
                             data.frame(offset= seq(-spread,spread,length.out=nrow(north.skew1.dat))))

# Make plots
BREAKS <- c(0.0,0.01,0.05,0.10,0.20,0.30,0.40,0.6,0.8)
x.labs <- c("None","Mock","Mock B")
x.at   <- c(1,2,3)


  shape.val = c(21,22)
  col.val = pal_jco(palette = c("default"), alpha = 1)(10)[c(1,4)]
  
north.skew.1.trim <-  ggplot() +
    geom_errorbar(data=north.skew1.dat %>% filter(type=="raw") ,aes(x=1+offset,
                      ymin=q.0.025, 
                      ymax=q.0.975,color=true.cat),width=0,alpha=0.75)   +
    geom_point(data=north.skew1.dat %>% filter(type=="raw"),
               aes(x=1+offset,y=Mean,shape=true.cat,fill=true.cat,color=true.cat),size=2) +
    geom_errorbar(data=north.skew1.dat %>% filter(type=="mock"),
                  aes(x=2+offset,
                      ymin= q.0.025, 
                      ymax= q.0.975,color=true.cat),width=0,alpha=0.75)   +
    geom_point(data=north.skew1.dat %>% filter(type=="mock") ,aes(x=2+offset,Mean,shape=true.cat,fill=true.cat,color=true.cat),size=2) +
    scale_shape_manual(values =shape.val) +
    scale_fill_manual(values= col.val, "True value") +
    scale_color_manual(values= col.val,"True value") +
    scale_y_continuous("Proportion",
                       trans="sqrt",
                       # trans="log",
                       breaks = BREAKS,expand=c(0,NA),limits = c(0,NA)) +
    geom_hline(data=north.skew1.dat,aes(yintercept = true.prop,color=true.cat),linetype="dashed") +
    #geom_point(aes(x=0.70,y=true.prop,shape=true.cat,fill=true.cat),size=3) +
    scale_x_continuous(name=NULL,breaks=x.at,labels = x.labs) +
    theme_classic() +
    theme(legend.position = "none")
  
north.skew.2.trim <-  ggplot() +
  geom_errorbar(data=north.skew2.dat %>% filter(type=="raw") ,aes(x=1+offset,
                                                                  ymin=q.0.025, 
                                                                  ymax=q.0.975,color=true.cat),width=0,alpha=0.75)   +
  geom_point(data=north.skew2.dat %>% filter(type=="raw"),
             aes(x=1+offset,y=Mean,shape=true.cat,fill=true.cat,color=true.cat),size=2) +
  geom_errorbar(data=north.skew2.dat %>% filter(type=="mock"),
                aes(x=2+offset,
                    ymin= q.0.025, 
                    ymax= q.0.975,color=true.cat),width=0,alpha=0.75)   +
  geom_point(data=north.skew2.dat %>% filter(type=="mock") ,aes(x=2+offset,Mean,shape=true.cat,fill=true.cat,color=true.cat),size=2) +
  scale_shape_manual(values =shape.val) +
  scale_fill_manual(values= col.val, "True value") +
  scale_color_manual(values= col.val,"True value") +
  scale_y_continuous("Proportion",
                     trans="sqrt",
                     # trans="log",
                     breaks = BREAKS,expand=c(0,NA),limits = c(0,NA)) +
  geom_hline(data=north.skew2.dat,aes(yintercept = true.prop,color=true.cat),linetype="dashed") +
  #geom_point(aes(x=0.70,y=true.prop,shape=true.cat,fill=true.cat),size=3) +
  scale_x_continuous(name=NULL,breaks=x.at,labels = x.labs) +
  theme_classic() +
  theme(legend.position = "none")



quartz(file="./plots/Ocean_Fish_North_trim_squalus.pdf",height=3,width=8,dpi=600,type="pdf")
  grid.arrange(north.skew.1.trim + ggtitle(NULL,subtitle="North Skew 1"),
             north.skew.2.trim + ggtitle(NULL,subtitle="North Skew 2"),
             ncol=2
  )
dev.off()


# Caluclate  Aitichson Distance based on posterior mean predictions.


A <- north.skew1.dat %>% filter(type=="raw")
B <- north.skew1.dat %>% filter(type=="mock")
aDist(A$true.prop,A$Mean)
#[1] 12.06688
aDist(B$true.prop,B$Mean)
#[1] 10.6452

C <- north.skew2.dat %>% filter(type=="raw")
D <- north.skew2.dat %>% filter(type=="mock")
aDist(C$true.prop,C$Mean)
#[1] 5.266042
aDist(D$true.prop,D$Mean)
#[1] 3.555512



