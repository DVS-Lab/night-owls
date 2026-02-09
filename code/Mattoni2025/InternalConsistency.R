library(tidyverse)
library(psych)

subj <- c('sub-101','sub-103','sub-104','sub-105')
lss <- read_delim('C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/derivatives/extractions/extractions_LSS_smoothed.tsv',
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



#Internal consistency across runs

mid_rel <- sr_rel <- list()

#mid
for (isub in subs) {
  for (irun in runs)
  for (ispace in spaces) {
    for (iecho in echoes) {
      for (iconfound in confounds) {
        
        # subset data
        subdat <- subset(lss_mid_rew, lss_mid_rew$sub == isub & lss_mid_rew$run==irun & lss_mid_rew$space == ispace &
                           lss_mid_rew$acq == iecho & lss_mid_rew$confounds == iconfound)  
        subdat_vs <- subdat[,c('sub','ses','run','NAcc_zstat_mean')]
        subdat_brs <- subdat[,c('sub','ses','run','BRS_corr')]
        
        colname <- paste0("sub-", isub,'_run-',irun, "_", ispace, "_",
                          iecho, "_", iconfound)
        
        # check if subset is empty
        if (nrow(subdat) == 0) {
          mid_rel[[colname]] <- NA
          next  
        }
        
        # reshape to wide 
        widedat_vs <- subdat_vs %>%
          group_by(sub, ses, run) %>%
          mutate(trial = row_number()) %>%
          ungroup() %>%
          pivot_wider(
            id_cols = c(sub, ses, run),
            names_from = trial,
            values_from = NAcc_zstat_mean,
            names_prefix = "NAcc_zstat_mean")
        widedat_brs <- subdat_brs %>%
          group_by(sub, ses, run) %>%
          mutate(trial = row_number()) %>%
          ungroup() %>%
          pivot_wider(
            id_cols = c(sub, ses, run),
            names_from = trial,
            values_from = BRS_corr,
            names_prefix = "BRS_corr")
        
        # keep only trial columns
        trials_vs <- widedat_vs[ , grepl("NAcc_zstat_mean", names(widedat_vs)), drop = FALSE]
        trials_brs <- widedat_brs[ , grepl("BRS_corr", names(widedat_brs)), drop = FALSE]
        
        
        # estimate split-half reliability across 1000 iterations
        split_corr_vs <- numeric(1000)
        split_corr_brs <- numeric(1000)
        
        for (i in 1:1000) {
          idx1 <- sample(1:ncol(trials_vs), ncol(trials_vs)/2)
          idx2 <- setdiff(1:ncol(trials_vs), idx1)
          
          # correlation between halves
          r_vs <- cor(rowSums(trials_vs[,idx1]), rowSums(trials_vs[,idx2]),use="pairwise.complete.obs")
          r_brs <- cor(rowSums(trials_brs[,idx1]), rowSums(trials_brs[,idx2]),use="pairwise.complete.obs")
          
          split_corr_vs[i] <- r_vs
          split_corr_brs[i] <- r_brs
          
        }
        
        mid_rel[[colname]]['VS'] <- mean(split_corr_vs)
        mid_rel[[colname]]['BRS'] <- mean(split_corr_brs)
        
  }}}
}


mid_rel_df <- bind_rows(lapply(names(mid_rel), function(nm) {
  vals <- mid_rel[[nm]]
  data.frame(name = nm, t(vals))
}), .id = NULL)
mid_rel_df <- mid_rel_df %>%
  separate(
    col = name,
    into = c("sub", "run", "space", "acq", "confounds"),
    sep = "_",remove = TRUE  )

#Average across runs
mid_rel_avg <- mid_rel_df %>%
  group_by(sub, space, acq, confounds) %>%
  summarise(VS = mean(VS, na.rm = TRUE),
    BRS = mean(BRS, na.rm = TRUE),.groups = "drop")
# Add group average
mid_rel_grp <- mid_rel_avg %>%
  group_by(space, acq, confounds) %>%
  summarise(VS=mean(VS, na.rm=TRUE),BRS=mean(BRS,na.rm=T)) %>%
  mutate(sub="Group") %>%
  ungroup()


#Shared Reward
for (isub in subs) {
for (irun in runs)
  for (ispace in spaces) {
    for (iecho in echoes) {
      for (iconfound in confounds) {
        
        # subset data
        subdat <- subset(lss_sr_rew, lss_sr_rew$sub == isub & lss_sr_rew$run == irun & lss_sr_rew$space == ispace &
                           lss_sr_rew$acq == iecho & lss_sr_rew$confounds == iconfound)  
        subdat_vs <- subdat[,c('sub','ses','run','NAcc_zstat_mean')]
        subdat_brs <- subdat[,c('sub','ses','run','BRS_corr')]
        
        colname <- paste0("sub-", isub, "_run-",irun, "_", ispace, "_",
                          iecho, "_", iconfound)
        
        # check if subset is empty
        if (nrow(subdat) == 0) {
          sr_rel[[colname]] <- NA
          next  
        }
        
        # reshape to wide 
        widedat_vs <- subdat_vs %>%
          group_by(sub, ses, run) %>%
          mutate(trial = row_number()) %>%  
          ungroup() %>%
          pivot_wider(
            id_cols = c(sub, ses, run),
            names_from = trial,
            values_from = NAcc_zstat_mean,
            names_prefix = "NAcc_zstat_mean")
        widedat_brs <- subdat_brs %>%
          group_by(sub, ses, run) %>%
          mutate(trial = row_number()) %>%  
          ungroup() %>%
          pivot_wider(
            id_cols = c(sub, ses, run),
            names_from = trial,
            values_from = BRS_corr,
            names_prefix = "BRS_corr")
        
        # keep only trial columns
        trials_vs <- widedat_vs[ , grepl("NAcc_zstat_mean", names(widedat_vs)), drop = FALSE]
        trials_brs <- widedat_brs[ , grepl("BRS_corr", names(widedat_brs)), drop = FALSE]
        
        
        # estimate split-half reliability across 1000 iterations
        split_corr_vs <- numeric(1000)
        split_corr_brs <- numeric(1000)
        
        for (i in 1:1000) {
          idx1 <- sample(1:ncol(trials_vs), ncol(trials_vs)/2)
          idx2 <- setdiff(1:ncol(trials_vs), idx1)
          
          # correlation between halves
          r_vs <- cor(rowSums(trials_vs[,idx1]), rowSums(trials_vs[,idx2]),use="pairwise.complete.obs")
          r_brs <- cor(rowSums(trials_brs[,idx1]), rowSums(trials_brs[,idx2]),use="pairwise.complete.obs")
          
          split_corr_vs[i] <- r_vs
          split_corr_brs[i] <- r_brs
          
        }
        
        sr_rel[[colname]]['VS'] <- mean(split_corr_vs)
        sr_rel[[colname]]['BRS'] <- mean(split_corr_brs)
        
  }}}
}


sr_rel_df <- bind_rows(lapply(names(sr_rel), function(nm) {
  vals <- sr_rel[[nm]]
  data.frame(name = nm, t(vals))
}), .id = NULL)
sr_rel_df <- sr_rel_df %>%
  separate(
    col = name,
    into = c("sub", "run", "space", "acq", "confounds"),
    sep = "_",remove = TRUE  )

#Average across runs
sr_rel_avg <- sr_rel_df %>%
  group_by(sub, space, acq, confounds) %>%
  summarise(VS = mean(VS, na.rm = TRUE),
            BRS = mean(BRS, na.rm = TRUE),.groups = "drop")
# Add group average
sr_rel_grp <- sr_rel_avg %>%
  group_by(space, acq, confounds) %>%
  summarise(VS=mean(VS, na.rm=TRUE),BRS=mean(BRS,na.rm=T)) %>%
  mutate(sub="Group") %>%
  ungroup()





#Number of Trials Affecting IC
mid_rel_ntrial <- sr_rel_ntrial <- list()

#mid
ntrials_mid <- seq(8,56,2)
for (isub in subs) {
    for (ispace in spaces) {
      for (iecho in echoes) {
        for (iconfound in confounds) {
          for (itrial in ntrials_mid){
          
          # subset data
          subdat <- subset(lss_mid_rew, lss_mid_rew$sub == isub & lss_mid_rew$space == ispace &
                             lss_mid_rew$acq == iecho & lss_mid_rew$confounds == iconfound)  
          subdat_vs <- subdat[,c('sub','ses','NAcc_zstat_mean','run')]
          subdat_vs$NAcc_zstat_mean_resid <- resid(lm(NAcc_zstat_mean ~ run,subdat_vs))
          subdat_brs <- subdat[,c('sub','ses','BRS_corr','run')]
          subdat_brs$BRS_corr_resid <- resid(lm(BRS_corr ~ run,subdat_brs))
          
          
          colname <- paste0("sub-", isub, "_", ispace, "_",
                            iecho, "_", iconfound,"_",itrial)
          
          # check if subset is empty
          if (nrow(subdat) == 0) {
            mid_rel_ntrial[[colname]] <- NA
            next  # skip to next iteration
          }
          
          # reshape to wide 
          widedat_vs <- subdat_vs %>%
            group_by(sub, ses) %>%
            mutate(trial = row_number()) %>%   # create a trial index
            ungroup() %>%
            pivot_wider(
              id_cols = c(sub, ses),
              names_from = trial,
              values_from = NAcc_zstat_mean_resid,
              names_prefix = "NAcc_zstat_mean")
          widedat_brs <- subdat_brs %>%
            group_by(sub, ses) %>%
            mutate(trial = row_number()) %>%   # create a trial index
            ungroup() %>%
            pivot_wider(
              id_cols = c(sub, ses),
              names_from = trial,
              values_from = BRS_corr_resid,
              names_prefix = "BRS_corr")
          
          # keep only trial columns
          trials_vs <- widedat_vs[ , grepl("NAcc_zstat_mean", names(widedat_vs)), drop = FALSE]
          trials_brs <- widedat_brs[ , grepl("BRS_corr", names(widedat_brs)), drop = FALSE]
          
          
          # estimate split-half reliability across 1000 iterations
          split_corr_vs <- numeric(1000)
          split_corr_brs <- numeric(1000)
          
          for (i in 1:1000) {
            cols <-  sample(ncol(trials_vs), itrial)
            trials_vs_cut <- trials_vs[,cols]
            trials_brs_cut <- trials_brs[,cols]
            idx1 <- sample(1:ncol(trials_vs_cut), ncol(trials_vs_cut)/2)
            idx2 <- setdiff(1:ncol(trials_vs_cut), idx1)

            
            # correlation between halves
            r_vs <- cor(rowSums(trials_vs[,idx1]), rowSums(trials_vs[,idx2]),use="pairwise.complete.obs")
            r_brs <- cor(rowSums(trials_brs[,idx1]), rowSums(trials_brs[,idx2]),use="pairwise.complete.obs")
            
            split_corr_vs[i] <- r_vs
            split_corr_brs[i] <- r_brs
            
          }
          
          mid_rel_ntrial[[colname]]['VS'] <- mean(split_corr_vs)
          mid_rel_ntrial[[colname]]['BRS'] <- mean(split_corr_brs)
        }}}}
}

mid_rel_ntrial_df <- bind_rows(lapply(names(mid_rel_ntrial), function(nm) {
  vals <- mid_rel_ntrial[[nm]]
  data.frame(name = nm, t(vals))
}), .id = NULL)
mid_rel_ntrial_df <- mid_rel_ntrial_df %>%
  separate(
    col = name,
    into = c("sub", "space", "acq", "confounds","ntrial"),
    sep = "_",remove = TRUE  )
#Focus ME tedana
mid_rel_ntrial_df <- mid_rel_ntrial_df[mid_rel_ntrial_df$space!="T1w" & mid_rel_ntrial_df$acq == "multiecho" &
                                         mid_rel_ntrial_df$confounds=="tedana",]
mid_rel_ntrial_df$ntrial <- as.numeric(mid_rel_ntrial_df$ntrial)
mid_rel_ntrial_df <- mid_rel_ntrial_df[,c('sub','ntrial','VS','BRS')]
# Add group average
mid_rel_ntrial_grp <- mid_rel_ntrial_df %>%
  group_by(ntrial) %>%
  summarise(VS=mean(VS, na.rm=TRUE),BRS=mean(BRS,na.rm=T))
mid_rel_ntrial_grp$sub <- 'Group'



#Shared Reward
ntrials_sr <- seq(8,44,2)
for (isub in subs) {
    for (ispace in spaces) {
      for (iecho in echoes) {
        for (iconfound in confounds) {
          for (itrial in ntrials_sr){
          
          # subset data
          subdat <- subset(lss_sr_rew, lss_sr_rew$sub == isub & lss_sr_rew$space == ispace &
                             lss_sr_rew$acq == iecho & lss_sr_rew$confounds == iconfound)  
          subdat_vs <- subdat[,c('sub','ses','NAcc_zstat_mean','run')]
          subdat_vs$NAcc_zstat_mean_resid <- resid(lm(NAcc_zstat_mean ~ run,subdat_vs))
          subdat_brs <- subdat[,c('sub','ses','BRS_corr','run')]
          subdat_brs$BRS_corr_resid <- resid(lm(BRS_corr ~ run,subdat_brs))
          
          colname <- paste0("sub-", isub, "_", ispace, "_",
                            iecho, "_", iconfound,'_',itrial)
          
          # check if subset is empty
          if (nrow(subdat) == 0) {
            sr_rel_ntrial[[colname]] <- NA
            next  # skip to next iteration
          }
          
          # reshape to wide 
          widedat_vs <- subdat_vs %>%
            group_by(sub, ses) %>%
            mutate(trial = row_number()) %>%   # create a trial index
            ungroup() %>%
            pivot_wider(
              id_cols = c(sub, ses),
              names_from = trial,
              values_from = NAcc_zstat_mean_resid,
              names_prefix = "NAcc_zstat_mean")
          widedat_brs <- subdat_brs %>%
            group_by(sub, ses) %>%
            mutate(trial = row_number()) %>%   # create a trial index
            ungroup() %>%
            pivot_wider(
              id_cols = c(sub, ses),
              names_from = trial,
              values_from = BRS_corr_resid,
              names_prefix = "BRS_corr")
          
          # keep only trial columns
          trials_vs <- widedat_vs[ , grepl("NAcc_zstat_mean", names(widedat_vs)), drop = FALSE]
          trials_brs <- widedat_brs[ , grepl("BRS_corr", names(widedat_brs)), drop = FALSE]
          
          
          # estimate split-half reliability across 1000 iterations
          split_corr_vs <- numeric(1000)
          split_corr_brs <- numeric(1000)
          
          for (i in 1:1000) {
            cols <-  sample(ncol(trials_vs), itrial)
            trials_vs_cut <- trials_vs[,cols]
            trials_brs_cut <- trials_brs[,cols]
            idx1 <- sample(1:ncol(trials_vs_cut), ncol(trials_vs_cut)/2)
            idx2 <- setdiff(1:ncol(trials_vs_cut), idx1)
            
            # correlation between halves
            r_vs <- cor(rowSums(trials_vs[,idx1]), rowSums(trials_vs[,idx2]),use="pairwise.complete.obs")
            r_brs <- cor(rowSums(trials_brs[,idx1]), rowSums(trials_brs[,idx2]),use="pairwise.complete.obs")
            
            split_corr_vs[i] <- r_vs
            split_corr_brs[i] <- r_brs
            
          }
          
          sr_rel_ntrial[[colname]]['VS'] <- mean(split_corr_vs)
          sr_rel_ntrial[[colname]]['BRS'] <- mean(split_corr_brs)
          
        }}}}
}

sr_rel_ntrial_df <- bind_rows(lapply(names(sr_rel_ntrial), function(nm) {
  vals <- sr_rel_ntrial[[nm]]
  data.frame(name = nm, t(vals))
}), .id = NULL)
sr_rel_ntrial_df <- sr_rel_ntrial_df %>%
  separate(
    col = name,
    into = c("sub", "space", "acq", "confounds","ntrial"),
    sep = "_",remove = TRUE  )
#Focus ME tedana
sr_rel_ntrial_df <- sr_rel_ntrial_df[sr_rel_ntrial_df$space!="T1w" & sr_rel_ntrial_df$acq == "multiecho" &
                                         sr_rel_ntrial_df$confounds=="tedana",]
sr_rel_ntrial_df$ntrial <- as.numeric(sr_rel_ntrial_df$ntrial)
# Add group average
sr_rel_ntrial_grp <- sr_rel_ntrial_df %>%
  group_by(ntrial) %>%
  summarise(VS=mean(VS, na.rm=TRUE),BRS=mean(BRS,na.rm=T))
sr_rel_ntrial_grp$sub <- 'Group'


