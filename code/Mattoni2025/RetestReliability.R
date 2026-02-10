library(tidyverse)
library(psych)
library(performance)
library(lme4)
library(lmerTest)

subj <- c('sub-101','sub-103','sub-104','sub-105')
l1stats <- read.delim('C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/derivatives/extractions/extractions_L1stats-revised_smoothed.tsv')
l1_mid <- l1stats[l1stats$task=='mid' & l1stats$label=='anticipation_reward>neutral',]
l1_sr <- l1stats[l1stats$task=='sharedreward' & l1stats$label=='S-C_rew>pun',]

subs <- substr(subj,5,8)
sessions <- seq(1:12)
runs <- seq(1:2)
spaces <- c("mni")
echoes <- c("single", "multiecho")
confounds <- c("fmriprep", "tedana")


#read motion
mriqc <- read.csv('C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/derivatives/data-outputs/mriqc/mriqc_metrics.csv')
mriqc$sub <- as.numeric(sub("sub-", "", mriqc$sub))
mriqc$ses <- as.numeric(sub("ses-", "", mriqc$ses))
mriqc$run <- as.numeric(sub("run-", "", mriqc$run))
mriqc$task <- sub("^task-", "", mriqc$task)
mriqc <- mriqc[mriqc$echo=='avg',]

mid_icc <- sr_icc <- mid_icc_unadj <- sr_icc_unadj <- list()
for (ispace in spaces) {
  for (iecho in echoes) {
    for (iconfound in confounds) {
      for (irun in 1:3){
        # subset data
        if (irun==1 | irun ==2){
          subdat <- subset(l1_mid, l1_mid$space == ispace & l1_mid$run ==irun &
                           l1_mid$acq == iecho & l1_mid$confounds == iconfound) 
          }
          else {
            subdat <- subset(l1_mid, l1_mid$space == ispace &
                               l1_mid$acq == iecho & l1_mid$confounds == iconfound)
        }
        colname <- paste0('run-',irun,'_',ispace, "_",
                          iecho, "_", iconfound)
        
        # check if subset is empty
        if (nrow(subdat) == 0) {
          mid_icc[[colname]] <- NA
          next  # skip to next iteration
        }
        subdat <- merge(subdat,mriqc,all.x=T,by=c('sub','ses','run','task'))
        if (irun==1 | irun ==2){
          mlm_vs <- lmer(NAcc_zstat_mean ~ ses + mean_fd + (1 | sub),subdat)
          mlm_brs <- lmer(BRS_corr ~ ses + mean_fd + (1 | sub),subdat) 
          }
        else {
          mlm_vs <- lmer(NAcc_zstat_mean ~ ses + run + mean_fd + (1 | sub),subdat)
          mlm_brs <- lmer(BRS_corr ~ ses + run + mean_fd + (1 | sub),subdat)
        }
        
      mid_icc[[colname]]['VS'] <- icc(mlm_vs)[1]
      mid_icc[[colname]]['BRS'] <- icc(mlm_brs)[1]
      
      mid_icc_unadj[[colname]]['VS'] <- icc(mlm_vs)[2]
      mid_icc_unadj[[colname]]['BRS'] <- icc(mlm_brs)[2]
      
}}}}


mid_icc_df <- bind_rows(lapply(names(mid_icc), function(nm) {
  vals <- mid_icc[[nm]]
  data.frame(name = nm, t(vals))
}), .id = NULL)
mid_icc_df <- mid_icc_df %>%
  separate(name, into = c("run", "space", "sequence", "confound"), sep = "_", remove = TRUE)

mid_icc_avg <- mid_icc_df %>%
  mutate(
    run = case_when(
      run %in% c("run-1", "run-2") ~ "avg_run12",
      run == "run-3" ~ "ses",
      TRUE ~ run),
    VS = as.numeric(VS),
    BRS = as.numeric(BRS)) %>%
  group_by(run, space, sequence, confound) %>%
  summarise(
    VS = mean(VS, na.rm = TRUE),
    BRS = mean(BRS, na.rm = TRUE),
    .groups = "drop")



#shared reward
for (ispace in spaces) {
  for (iecho in echoes) {
    for (iconfound in confounds) {
      for (irun in 1:3){
        # subset data
        if (irun==1 | irun ==2){
          subdat <- subset(l1_sr, l1_sr$space == ispace & l1_sr$run ==irun &
                             l1_sr$acq == iecho & l1_sr$confounds == iconfound) 
        }
        else {
          subdat <- subset(l1_sr, l1_sr$space == ispace &
                             l1_sr$acq == iecho & l1_sr$confounds == iconfound)
        }
        colname <- paste0('run-',irun,'_',ispace, "_",
                          iecho, "_", iconfound)
        
        # check if subset is empty
        if (nrow(subdat) == 0) {
          sr_icc[[colname]] <- NA
          next  # skip to next iteration
        }
        subdat <- merge(subdat,mriqc,all.x=T,by=c('sub','ses','run','task'))
        if (irun==1 | irun ==2){
          mlm_vs <- lmer(NAcc_zstat_mean ~ ses + mean_fd + (1 | sub),subdat)
          mlm_brs <- lmer(BRS_corr ~ ses + mean_fd +  (1 | sub),subdat) 
        }
        else {
          mlm_vs <- lmer(NAcc_zstat_mean ~ ses + run + mean_fd + (1 | sub),subdat)
          mlm_brs <- lmer(BRS_corr ~ ses + run + mean_fd + (1 | sub),subdat)
        }
        
        sr_icc[[colname]]['VS'] <- icc(mlm_vs)[1]
        sr_icc[[colname]]['BRS'] <- icc(mlm_brs)[1]
        
        sr_icc_unadj[[colname]]['VS'] <- icc(mlm_vs)[2]
        sr_icc_unadj[[colname]]['BRS'] <- icc(mlm_brs)[2]
        
      }}}}


sr_icc_df <- bind_rows(lapply(names(sr_icc), function(nm) {
  vals <- sr_icc[[nm]]
  data.frame(name = nm, t(vals))
}), .id = NULL)
sr_icc_df <- sr_icc_df %>%
  separate(name, into = c("run", "space", "sequence", "confound"), sep = "_", remove = TRUE)

sr_icc_avg <- sr_icc_df %>%
  mutate(
    run = case_when(
      run %in% c("run-1", "run-2") ~ "avg_run12",
      run == "run-3" ~ "ses",
      TRUE ~ run),
    VS = as.numeric(VS),
    BRS = as.numeric(BRS)) %>%
  group_by(run, space, sequence, confound) %>%
  summarise(
    VS = mean(VS, na.rm = TRUE),
    BRS = mean(BRS, na.rm = TRUE),
    .groups = "drop")


mid_icc_avg$acq_cnfd <- paste0(mid_icc_avg$sequence,'_',mid_icc_avg$confound)
mid_icc_avg <- mid_icc_avg[mid_icc_avg$acq_cnfd!="single_tedana",]
mid_icc_avg <- mid_icc_avg %>%
  pivot_longer(cols = c(VS, BRS),
               names_to = "type",
               values_to = "estimate")

mid_icc_avg$acq_cnfd <- factor(mid_icc_avg$acq_cnfd)
mid_icc_avg$pattern_type <- ifelse(mid_icc_avg$type == "BRS", "stripe", "none")

mid_icc_avg$acq_cnfd <- factor(
  mid_icc_avg$acq_cnfd,
  levels = c("single_fmriprep", "multiecho_fmriprep", "multiecho_tedana"),
  labels = c("A", "B", "C"))

mid_icc_avg$estimate_signed <- ifelse(mid_icc_avg$type == "BRS",
                                      -mid_icc_avg$estimate,
                                      mid_icc_avg$estimate)
mid_icc_avg$type <- factor(mid_icc_avg$type, levels = c("VS", "BRS"))
mid_icc_avg$run[mid_icc_avg$run == "ses"] <- "Session-Level"
mid_icc_avg$run[mid_icc_avg$run == "avg_run12"] <- "Run-Level Average"



mean(mid_icc_avg$estimate[mid_icc_avg$run=='Run-Level Average' & mid_icc_avg$acq_cnfd=='A' & mid_icc_avg$type=='VS'],na.rm=T)
mean(mid_icc_avg$estimate[mid_icc_avg$run=='Run-Level Average' & mid_icc_avg$acq_cnfd=='B' & mid_icc_avg$type=='VS'],na.rm=T)
mean(mid_icc_avg$estimate[mid_icc_avg$run=='Run-Level Average' & mid_icc_avg$acq_cnfd=='C' & mid_icc_avg$type=='VS'],na.rm=T)

mean(mid_icc_avg$estimate[mid_icc_avg$run=='Run-Level Average' & mid_icc_avg$acq_cnfd=='A' & mid_icc_avg$type=='BRS'],na.rm=T)
mean(mid_icc_avg$estimate[mid_icc_avg$run=='Run-Level Average' & mid_icc_avg$acq_cnfd=='B' & mid_icc_avg$type=='BRS'],na.rm=T)
mean(mid_icc_avg$estimate[mid_icc_avg$run=='Run-Level Average' & mid_icc_avg$acq_cnfd=='C' & mid_icc_avg$type=='BRS'],na.rm=T)

mean(mid_icc_avg$estimate[mid_icc_avg$run=='Session-Level' & mid_icc_avg$acq_cnfd=='A' & mid_icc_avg$type=='VS'],na.rm=T)
mean(mid_icc_avg$estimate[mid_icc_avg$run=='Session-Level' & mid_icc_avg$acq_cnfd=='B' & mid_icc_avg$type=='VS'],na.rm=T)
mean(mid_icc_avg$estimate[mid_icc_avg$run=='Session-Level' & mid_icc_avg$acq_cnfd=='C' & mid_icc_avg$type=='VS'],na.rm=T)

mean(mid_icc_avg$estimate[mid_icc_avg$run=='Session-Level' & mid_icc_avg$acq_cnfd=='A' & mid_icc_avg$type=='BRS'],na.rm=T)
mean(mid_icc_avg$estimate[mid_icc_avg$run=='Session-Level' & mid_icc_avg$acq_cnfd=='B' & mid_icc_avg$type=='BRS'],na.rm=T)
mean(mid_icc_avg$estimate[mid_icc_avg$run=='Session-Level' & mid_icc_avg$acq_cnfd=='C' & mid_icc_avg$type=='BRS'],na.rm=T)



sr_icc_avg$acq_cnfd <- paste0(sr_icc_avg$sequence,'_',sr_icc_avg$confound)
sr_icc_avg <- sr_icc_avg[sr_icc_avg$acq_cnfd!="single_tedana",]
sr_icc_avg <- sr_icc_avg %>%
  pivot_longer(cols = c(VS, BRS),
               names_to = "type",
               values_to = "estimate")

sr_icc_avg$acq_cnfd <- factor(sr_icc_avg$acq_cnfd)
sr_icc_avg$pattern_type <- ifelse(sr_icc_avg$type == "BRS", "stripe", "none")

sr_icc_avg$acq_cnfd <- factor(
  sr_icc_avg$acq_cnfd,
  levels = c("single_fmriprep", "multiecho_fmriprep", "multiecho_tedana"),
  labels = c("A", "B", "C"))

sr_icc_avg$estimate_signed <- ifelse(sr_icc_avg$type == "BRS",
                                     -sr_icc_avg$estimate,
                                     sr_icc_avg$estimate)
sr_icc_avg$type <- factor(sr_icc_avg$type, levels = c("VS", "BRS"))
sr_icc_avg$run[sr_icc_avg$run == "ses"] <- "Session-Level"
sr_icc_avg$run[sr_icc_avg$run == "avg_run12"] <- "Run-Level Average"


mean(sr_icc_avg$estimate[sr_icc_avg$run=='Run-Level Average' & sr_icc_avg$acq_cnfd=='A' & sr_icc_avg$type=='VS'],na.rm=T)
mean(sr_icc_avg$estimate[sr_icc_avg$run=='Run-Level Average' & sr_icc_avg$acq_cnfd=='B' & sr_icc_avg$type=='VS'],na.rm=T)
mean(sr_icc_avg$estimate[sr_icc_avg$run=='Run-Level Average' & sr_icc_avg$acq_cnfd=='C' & sr_icc_avg$type=='VS'],na.rm=T)

mean(sr_icc_avg$estimate[sr_icc_avg$run=='Run-Level Average' & sr_icc_avg$acq_cnfd=='A' & sr_icc_avg$type=='BRS'],na.rm=T)
mean(sr_icc_avg$estimate[sr_icc_avg$run=='Run-Level Average' & sr_icc_avg$acq_cnfd=='B' & sr_icc_avg$type=='BRS'],na.rm=T)
mean(sr_icc_avg$estimate[sr_icc_avg$run=='Run-Level Average' & sr_icc_avg$acq_cnfd=='C' & sr_icc_avg$type=='BRS'],na.rm=T)

mean(sr_icc_avg$estimate[sr_icc_avg$run=='Session-Level' & sr_icc_avg$acq_cnfd=='A' & sr_icc_avg$type=='VS'],na.rm=T)
mean(sr_icc_avg$estimate[sr_icc_avg$run=='Session-Level' & sr_icc_avg$acq_cnfd=='B' & sr_icc_avg$type=='VS'],na.rm=T)
mean(sr_icc_avg$estimate[sr_icc_avg$run=='Session-Level' & sr_icc_avg$acq_cnfd=='C' & sr_icc_avg$type=='VS'],na.rm=T)

mean(sr_icc_avg$estimate[sr_icc_avg$run=='Session-Level' & sr_icc_avg$acq_cnfd=='A' & sr_icc_avg$type=='BRS'],na.rm=T)
mean(sr_icc_avg$estimate[sr_icc_avg$run=='Session-Level' & sr_icc_avg$acq_cnfd=='B' & sr_icc_avg$type=='BRS'],na.rm=T)
mean(sr_icc_avg$estimate[sr_icc_avg$run=='Session-Level' & sr_icc_avg$acq_cnfd=='C' & sr_icc_avg$type=='BRS'],na.rm=T)

