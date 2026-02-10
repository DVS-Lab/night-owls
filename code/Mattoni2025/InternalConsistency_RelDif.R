library(tidyverse)
library(psych)
library(ggplot2)
library(ggpattern)
library(patchwork)
library(viridis)

subj <- c('sub-101','sub-103','sub-104','sub-105')
lss <- read_delim("C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/derivatives/extractions/extractions_LSS_smoothed.tsv", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
lss <- lss[do.call(order, lss[, 1:8]), ]

lss$label<-NA
lss$ses <- as.numeric(lss$ses)
lss_mid <- lss[lss$task=='mid',]
lss_sr <- lss[lss$task=='sharedreward',]


subs <- substr(subj,5,8)
sessions <- seq(1:12)
runs <- seq(1:2)
#spaces <- c("MNI152NLin6Asym", "T1w")
spaces <- c("MNI152NLin6Asym")

echoes <- c("single", "multiecho")
confounds <- c("base", "tedana")


#Assign labels
##mid
for (isub in subs) {
for (ises in sessions) {
for (irun in runs) {
  
  filepath <- paste0('C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/stimuli/mid/data/sub-',
                     isub,'/sub-',isub,'_task-mid_ses-',ises,'_run-',irun,'.csv')
  if (!file.exists(filepath)) next
  log <- read.csv(filepath)
  log <- log[-nrow(log),]
  
  for (ispace in spaces) {
  for (iecho in echoes) {
  for (iconfound in confounds) {
  
  idx <- which(lss_mid$sub==isub & lss_mid$ses==ises & lss_mid$run==irun &
                 lss_mid$space==ispace & lss_mid$acq==iecho & lss_mid$confounds==iconfound)[1]
  if (is.na(idx)) next
  
  length <- length(which(lss_mid$sub==isub & lss_mid$ses==ises & lss_mid$run==irun &
                            lss_mid$space==ispace & lss_mid$acq==iecho & lss_mid$confounds==iconfound))
  lss_mid$label[idx:(idx+length-1)] <- log$cue.color
  }}}
}}}
###subset to rew trials
lss_mid_rew <- lss_mid[lss_mid$label=='Green',]

##sharedreward
for (isub in subs) {
  for (ises in sessions) {
    for (irun in runs) {
      ises<-as.numeric(ises)
      filepath <- paste0('C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/stimuli/sharedreward/logs/sub-',
                         isub,'/sub-',isub,'_task-sharedreward_ses-',ises,'_run-',irun,'_raw.csv')
      if (!file.exists(filepath)) next
      log <- read.csv(filepath)

      for (ispace in spaces) {
        for (iecho in echoes) {
          for (iconfound in confounds) {
            
            idx <- which(lss_sr$sub==isub & lss_sr$ses==ises & lss_sr$run==irun &
                           lss_sr$space==ispace & lss_sr$acq==iecho & lss_sr$confounds==iconfound)[1]
            if (is.na(idx)) next
            
            length <- length(which(lss_sr$sub==isub & lss_sr$ses==ises & lss_sr$run==irun &
                                     lss_sr$space==ispace & lss_sr$acq==iecho & lss_sr$confounds==iconfound))
            
            length_replacement <- length(log$Feedback)
            
            if (length != length_replacement) {
              message("Length mismatch at iteration ", 
                      ": target=", ", replacement=", length_replacement)}
              
              
            lss_sr$label[idx:(idx+length-1)] <- log$Feedback
          }}}
    }}}
###Remove trial with missed response
lss_sr <- lss_sr[-which(lss_sr$sub==105 & lss_sr$ses==2 & lss_sr$run==1 & lss_sr$trial==22),]
###subset to rew trials
lss_sr_rew <- lss_sr[lss_sr$label==3,]

#Num trials
#Internal consistency across runs
mid_rel_ntrials_vs <- mid_rel_ntrials_brs <- sr_rel_ntrials_vs <- sr_rel_ntrials_brs <- list()

# mid
ntrials <- seq(7, 336, 7)

for (isub in subs) {
  
  # subset data
  subdat <- subset(lss_mid_rew, lss_mid_rew$sub == isub & lss_mid_rew$space == "MNI152NLin6Asym" &
                     lss_mid_rew$acq == "multiecho" & lss_mid_rew$confounds == "tedana")  
  subdat_vs <- subdat[,c('sub','ses','NAcc_zstat_mean','run')]
  subdat_vs <- subdat_vs %>%
    arrange(ses, run) %>%
    mutate(
      rownum = seq(1:nrow(subdat_vs)),
      scan = (ses - 1) * 2 + run)
  
  subdat_brs <- subdat[,c('sub','ses','BRS_corr','run')]
  subdat_brs <- subdat_brs %>%
    arrange(ses, run) %>%
    mutate(
      rownum = seq(1:nrow(subdat_brs)),
      scan = (ses - 1) * 2 + run)
  
  for (itrial in ntrials) {
    
    results_vs <- results_brs <- numeric(1000)

    for (i in 1:1000) {
      
      # Holdout: 6 full runs
      holdout <- sample(subdat_vs$rownum, 168)
      
      holdout_mean_vs <- subdat_vs %>%
        filter(rownum %in% holdout) %>%
        pull(NAcc_zstat_mean) %>%
        mean()
      
      holdout_sd_vs <- subdat_vs %>%
        filter(rownum %in% holdout) %>%
        pull(NAcc_zstat_mean) %>%
        sd()
      
      holdout_mean_brs <- subdat_brs %>%
        filter(rownum %in% holdout) %>%
        pull(BRS_corr) %>%
        mean()
      
      holdout_sd_brs <- subdat_brs %>%
        filter(rownum %in% holdout) %>%
        pull(BRS_corr) %>%
        sd()
      
      remaining <- setdiff(subdat_vs$rownum, holdout)
      comparison <- sample(remaining, itrial)
      
      comparison_trials_vs <- subdat_vs %>%
        filter(rownum %in% comparison) %>%
        pull(NAcc_zstat_mean)
      
      comparison_mean_vs <- mean(comparison_trials_vs)
      
      results_vs[i] <- abs((comparison_mean_vs - holdout_mean_vs) / holdout_sd_vs)
      
      comparison_trials_brs <- subdat_brs %>%
        filter(rownum %in% comparison) %>%
        pull(BRS_corr)
      
      comparison_mean_brs <- mean(comparison_trials_brs)
      
      results_brs[i] <- abs((comparison_mean_brs - holdout_mean_brs) / holdout_sd_brs)
    }
    
    # Store results
    colname <- paste(isub, itrial, sep = "_")
    mid_rel_ntrials_vs[[colname]] <- mean(results_vs)
    mid_rel_ntrials_brs[[colname]] <- mean(results_brs)
    
  }}



# sr
ntrials <- seq(5, 264, 7)

for (isub in subs) {
  
  # subset data
  subdat <- subset(lss_sr_rew, lss_sr_rew$sub == isub & lss_sr_rew$space == "MNI152NLin6Asym" &
                     lss_sr_rew$acq == "multiecho" & lss_sr_rew$confounds == "tedana")  
  subdat_vs <- subdat[,c('sub','ses','NAcc_zstat_mean','run')]
  subdat_vs <- subdat_vs %>%
    arrange(ses, run) %>%
    mutate(
      rownum = seq(1:nrow(subdat_vs)),
      scan = (ses - 1) * 2 + run)
  
  subdat_brs <- subdat[,c('sub','ses','BRS_corr','run')]
  subdat_brs <- subdat_brs %>%
    arrange(ses, run) %>%
    mutate(
      rownum = seq(1:nrow(subdat_brs)),
      scan = (ses - 1) * 2 + run)
  
  scans <- unique(subdat_vs$scan)
  
  for (itrial in ntrials) {
    
    results_vs <- results_brs <- numeric(1000)

    for (i in 1:1000) {
      
      # Holdout: 6 full runs
      holdout <- sample(subdat_vs$rownum, 132)
      
      holdout_mean_vs <- subdat_vs %>%
        filter(rownum %in% holdout) %>%
        pull(NAcc_zstat_mean) %>%
        mean()
      
      holdout_sd_vs <- subdat_vs %>%
        filter(rownum %in% holdout) %>%
        pull(NAcc_zstat_mean) %>%
        sd()
      
      holdout_mean_brs <- subdat_brs %>%
        filter(rownum %in% holdout) %>%
        pull(BRS_corr) %>%
        mean()
      
      holdout_sd_brs <- subdat_brs %>%
        filter(rownum %in% holdout) %>%
        pull(BRS_corr) %>%
        sd()
      
      remaining <- setdiff(subdat_vs$rownum, holdout)
      comparison <- sample(remaining, itrial)
      
      comparison_trials_vs <- subdat_vs %>%
        filter(rownum %in% comparison) %>%
        pull(NAcc_zstat_mean)
      
      comparison_mean_vs <- mean(comparison_trials_vs)
      
      results_vs[i] <- abs((comparison_mean_vs - holdout_mean_vs) / holdout_sd_vs)
      
      comparison_trials_brs <- subdat_brs %>%
        filter(rownum %in% comparison) %>%
        pull(BRS_corr)
      
      comparison_mean_brs <- mean(comparison_trials_brs)
      
      results_brs[i] <- abs((comparison_mean_brs - holdout_mean_brs) / holdout_sd_brs)
    }
    
    # Store results
    colname <- paste(isub, itrial, sep = "_")
    sr_rel_ntrials_vs[[colname]] <- mean(results_vs)
    sr_rel_ntrials_brs[[colname]] <- mean(results_brs)
    
  }}


# plot
mid_vs_df <- data.frame(
  name = names(mid_rel_ntrials_vs),
  rel = unlist(mid_rel_ntrials_vs, use.names = FALSE)) %>%
  separate(name, into = c("sub", "trials"), sep = "_", remove = FALSE) %>%
  mutate(
    sub = as.integer(sub),
    trials = as.integer(trials)) %>%
  select(sub, trials, rel)

subject_colors_num <- c(
  "101" = "#fde725",
  "103" = "#35b779",
  "104" = "#31688e",
  "105" = "#440154",
  "Group" = "black")

mid_vs_p <- ggplot(mid_vs_df, aes(x = trials, y = rel, color = factor(sub), group = sub)) +
  geom_line(size=1.5,alpha=.7) +
  geom_vline(xintercept = 28, linetype = "dashed", color = "black") +
  annotate("text", x = 28, y = max(mid_vs_df$rel), label = "Run", 
           hjust = -0.1, vjust = 1, size=4) +
  geom_vline(xintercept = 56, linetype = "dashed", color = "black") +
  annotate("text", x = 57, y = max(mid_vs_df$rel-.05), label = "Session", 
           hjust = -0.1, vjust = 1, size=4) +
  scale_color_manual(values = subject_colors_num) +
  scale_x_continuous(breaks = seq(0, 300, 50)) +
  ylim(0,.40) +
  xlim(0,250) +
  labs(title="NAcc",
    x = "Trials",
    y = "MID Split-Half Difference\n(Holdout SDs)",
    color = "Subject") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(face = "bold", size = 14),
        plot.title = element_text(hjust = 0.5))


mid_brs_df <- data.frame(
  name = names(mid_rel_ntrials_brs),
  rel = unlist(mid_rel_ntrials_brs, use.names = FALSE)) %>%
  separate(name, into = c("sub", "trials"), sep = "_", remove = FALSE) %>%
  mutate(
    sub = as.integer(sub),
    trials = as.integer(trials)) %>%
  select(sub, trials, rel)

mid_brs_p <- ggplot(mid_brs_df, aes(x = trials, y = rel, color = factor(sub), group = sub)) +
  geom_line(size=1.5,alpha=.7) +
  geom_vline(xintercept = 28, linetype = "dashed", color = "black") +
  annotate("text", x = 28, y = max(mid_brs_df$rel), label = "Run", 
           hjust = -0.1, vjust = 1, size=4) +
  geom_vline(xintercept = 56, linetype = "dashed", color = "black") +
  annotate("text", x = 57, y = max(mid_brs_df$rel-.05), label = "Session", 
           hjust = -0.1, vjust = 1, size=4) +
  scale_color_manual(values = subject_colors_num) +
  scale_x_continuous(breaks = seq(0, 300, 50)) +
  ylim(0,.40) +
  xlim(0,250) +
  labs(title="BRS",
    x = "Trials",
    y = NULL,
    color = "Subject") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(face = "bold", size = 14),
        plot.title = element_text(hjust = 0.5))
mid_p <- mid_vs_p + mid_brs_p
mid_p
ggsave("C:/Users/mmatt/Desktop/Projects/NightOwls/NOSC-Analysis/Revised/mid_reldif.png",
       height = 6, width = 8, units = "in")



# plot
sr_vs_df <- data.frame(
  name = names(sr_rel_ntrials_vs),
  rel = unlist(sr_rel_ntrials_vs, use.names = FALSE)) %>%
  separate(name, into = c("sub", "trials"), sep = "_", remove = FALSE) %>%
  mutate(
    sub = as.integer(sub),
    trials = as.integer(trials)) %>%
  select(sub, trials, rel)

sr_vs_p <- ggplot(sr_vs_df, aes(x = trials, y = rel, color = factor(sub), group = sub)) +
  geom_line(size=1.5,alpha=.7) +
  geom_vline(xintercept = 22, linetype = "dashed", color = "black") +
  annotate("text", x = 22, y = max(sr_vs_df$rel), label = "Run", 
           hjust = -0.1, vjust = 1, size=4) +
  geom_vline(xintercept = 44, linetype = "dashed", color = "black") +
  annotate("text", x = 45, y = max(sr_vs_df$rel-.05), label = "Session", 
           hjust = -0.1, vjust = 1, size=4) +
  scale_color_manual(values = subject_colors_num) +
  scale_x_continuous(breaks = seq(0, 300, 50)) +
  ylim(0,.40) +
  xlim(0,250) +
  labs(title="NAcc",
       x = "Trials",
       y = "SR Split-Half Difference\n(Holdout SDs)",
       color = "Subject") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(face = "bold", size = 14),
        plot.title = element_text(hjust = 0.5))


sr_brs_df <- data.frame(
  name = names(sr_rel_ntrials_brs),
  rel = unlist(sr_rel_ntrials_brs, use.names = FALSE)) %>%
  separate(name, into = c("sub", "trials"), sep = "_", remove = FALSE) %>%
  mutate(
    sub = as.integer(sub),
    trials = as.integer(trials)) %>%
  select(sub, trials, rel)

sr_brs_p <- ggplot(sr_brs_df, aes(x = trials, y = rel, color = factor(sub), group = sub)) +
  geom_line(size=1.5,alpha=.7) +
  geom_vline(xintercept = 22, linetype = "dashed", color = "black") +
  annotate("text", x = 22, y = max(sr_brs_df$rel), label = "Run", 
           hjust = -0.1, vjust = 1, size=4) +
  geom_vline(xintercept = 44, linetype = "dashed", color = "black") +
  annotate("text", x = 45, y = max(sr_brs_df$rel-.05), label = "Session", 
           hjust = -0.1, vjust = 1, size=4) +
  scale_color_manual(values = subject_colors_num) +
  scale_x_continuous(breaks = seq(0, 300, 50)) +
  ylim(0,.40) +
  xlim(0,250) +
  labs(title="BRS",
       x = "Trials",
       y = NULL,
       color = "Subject") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(face = "bold", size = 14),
        plot.title = element_text(hjust = 0.5))
sr_p <- sr_vs_p + sr_brs_p
sr_p
ggsave("C:/Users/mmatt/Desktop/Projects/NightOwls/NOSC-Analysis/Revised/sr_reldif.png",
       height = 6, width = 8, units = "in")

