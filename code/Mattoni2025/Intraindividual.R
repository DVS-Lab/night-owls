library(tidyverse)
library(psych)
library(performance)
library(lme4)
library(lmerTest)
library(MuMIn)
library(partR2)


subj <- c('sub-101','sub-103','sub-104','sub-105')
l1stats <- read.delim('C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/derivatives/extractions/extractions_L1stats-revised_smoothed.tsv')
l1_mid <- l1stats[l1stats$task=='mid' & l1stats$label=='anticipation_reward>neutral',]
l1_sr <- l1stats[l1stats$task=='sharedreward' & l1stats$label=='S-C_rew>pun',]

l1_mid_prec <- l1_mid[l1_mid$space=="mni" & l1_mid$acq=="multiecho" & l1_mid$confounds=="tedana",]
l1_sr_prec <- l1_sr[l1_sr$space=="mni" & l1_sr$acq=="multiecho" & l1_sr$confounds=="tedana",]


subs <- substr(subj,5,8)
sessions <- seq(1:12)
runs <- seq(1:2)
spaces <- c("mni")
echoes <- c("single", "multiecho")
confounds <- c("fmriprep", "tedana")

behav_long <- read.csv('C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/derivatives/data-outputs/behavioral_data.csv')


#read motion
mriqc <- read.csv('C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/derivatives/data-outputs/mriqc/mriqc_metrics.csv')
mriqc$sub <- as.numeric(sub("sub-", "", mriqc$sub))
mriqc$ses <- as.numeric(sub("ses-", "", mriqc$ses))
mriqc$run <- as.numeric(sub("run-", "", mriqc$run))
mriqc$task <- sub("^task-", "", mriqc$task)
mriqc <- mriqc[mriqc$task!="rest",]
mriqc <- mriqc[mriqc$echo=='avg',]

l1_mid_prec <- merge(l1_mid_prec,mriqc[mriqc$task=="mid",],by=c("sub","ses","run","task"),all.x=T)
l1_sr_prec <- merge(l1_sr_prec,mriqc[mriqc$task=="sharedreward",],by=c("sub","ses","run","task"),all.x=T)

scan_notes <- read.csv('C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/derivatives/data-outputs/scan_sessions.csv')
behav_long$sub <- as.numeric(substr(behav_long$subj,5,8))
behav_long_cb <- merge(behav_long,scan_notes,by=c('sub','ses'))

l1_mid_prec_behav <- merge(l1_mid_prec,behav_long_cb,by=c('sub','ses'))
l1_mid_intra <- l1_mid_prec_behav[l1_mid_prec_behav$counterbalance=='A' & l1_mid_prec_behav$obs==1 & l1_mid_prec_behav$run==1,]
l1_mid_intra <- rbind(l1_mid_intra,l1_mid_prec_behav[l1_mid_prec_behav$counterbalance=='A' & l1_mid_prec_behav$obs==4 & l1_mid_prec_behav$run==2,])
l1_mid_intra <- rbind(l1_mid_intra,l1_mid_prec_behav[l1_mid_prec_behav$counterbalance=='B' & l1_mid_prec_behav$obs==2 & l1_mid_prec_behav$run==1,])
l1_mid_intra <- rbind(l1_mid_intra,l1_mid_prec_behav[l1_mid_prec_behav$counterbalance=='B' & l1_mid_prec_behav$obs==5 & l1_mid_prec_behav$run==2,])
l1_mid_intra <- l1_mid_intra[order(l1_mid_intra$sub, l1_mid_intra$ses, l1_mid_intra$run),]
l1_mid_intra$ses_f <- factor(l1_mid_intra$ses)
l1_mid_intra$sub <- factor(l1_mid_intra$sub)
l1_mid_intra$kss.r <- 9-l1_mid_intra$kss
l1_mid_intra$run_f <- factor(l1_mid_intra$run)



l1_sr_prec_behav <- merge(l1_sr_prec,behav_long_cb,by=c('sub','ses'))
l1_sr_intra <- l1_sr_prec_behav[l1_sr_prec_behav$counterbalance=='A' & l1_sr_prec_behav$obs==2 & l1_sr_prec_behav$run==1,]
l1_sr_intra <- rbind(l1_sr_intra,l1_sr_prec_behav[l1_sr_prec_behav$counterbalance=='A' & l1_sr_prec_behav$obs==5 & l1_sr_prec_behav$run==2,])
l1_sr_intra <- rbind(l1_sr_intra,l1_sr_prec_behav[l1_sr_prec_behav$counterbalance=='B' & l1_sr_prec_behav$obs==1 & l1_sr_prec_behav$run==1,])
l1_sr_intra <- rbind(l1_sr_intra,l1_sr_prec_behav[l1_sr_prec_behav$counterbalance=='B' & l1_sr_prec_behav$obs==4 & l1_sr_prec_behav$run==2,])
l1_sr_intra <- l1_sr_intra[order(l1_sr_intra$sub, l1_sr_intra$ses, l1_sr_intra$run),]
l1_sr_intra$ses_f <- factor(l1_sr_intra$ses)
l1_sr_intra$sub <- factor(l1_sr_intra$sub)
l1_sr_intra$kss.r <- 9-l1_sr_intra$kss
l1_sr_intra$run_f <- factor(l1_sr_intra$run)




#Effect of session
##WP Standardize
l1_mid_intra <- l1_mid_intra %>%
  group_by(sub) %>%
  mutate(VS.s = as.numeric(scale(VS_mean)),
         BRS.s = as.numeric(scale(BRS_corr)),
         ses.s = as.numeric(scale(ses)),
         run.s = as.numeric(scale(run)),
         kss.s = as.numeric(scale(kss.r)),
         mean_fd.s = as.numeric(scale(mean_fd)),
         mood.s = as.numeric(scale(mood)),
         panas_pa.s = as.numeric(scale(panas_pa))) %>%
  ungroup()

l1_sr_intra <- l1_sr_intra %>%
  group_by(sub) %>%
  mutate(VS.s = as.numeric(scale(VS_mean)),
         BRS.s = as.numeric(scale(BRS_corr)),
         ses.s = as.numeric(scale(ses)),
         run.s = as.numeric(scale(run)),
         kss.s = as.numeric(scale(kss.r)),
         mean_fd.s = as.numeric(scale(mean_fd)),
         mood.s = as.numeric(scale(mood)),
         panas_pa.s = as.numeric(scale(panas_pa))) %>%
  ungroup()


##MLMs
mid_changes_vs <- lmer(VS.s ~ ses.s + run.s + mean_fd.s + (1 | sub) + (1 | ses.s), l1_mid_intra)
summary(mid_changes_vs)
mid_changes_brs <- lmer(BRS.s ~ ses.s + run.s + mean_fd.s + (1 | sub) + (1 | ses.s), l1_mid_intra)
summary(mid_changes_brs)

sr_changes_vs <- lmer(VS.s ~ ses.s + run.s + mean_fd.s + (1 | sub) + (1 | ses.s), l1_sr_intra)
summary(sr_changes_vs)
sr_changes_brs <- lmer(BRS.s ~ ses.s + run.s + mean_fd.s + (1 | sub) + (1 | ses.s), l1_sr_intra)
summary(sr_changes_brs)


##Ind subject
mid_intra_sub <- split(l1_mid_intra, l1_mid_intra$sub)
mid_intra_sub_mlm <- map(mid_intra_sub, ~ {lmer(VS.s ~ ses.s + run.s + mean_fd.s + (1 | ses.s),data = .x)})
summary(mid_intra_sub_mlm[[1]])
summary(mid_intra_sub_mlm[[2]])
summary(mid_intra_sub_mlm[[3]])
summary(mid_intra_sub_mlm[[4]])

mid_intra_sub_mlm <- map(mid_intra_sub, ~ {lmer(BRS.s ~ ses.s + run.s + mean_fd.s + (1 | ses.s),data = .x)})
summary(mid_intra_sub_mlm[[1]])
summary(mid_intra_sub_mlm[[2]])
summary(mid_intra_sub_mlm[[3]])
summary(mid_intra_sub_mlm[[4]])

sr_intra_sub <- split(l1_sr_intra, l1_sr_intra$sub)
sr_intra_sub_mlm <- map(sr_intra_sub, ~ {lmer(VS.s ~ run.s + ses.s + mean_fd.s + (1 | ses.s),data = .x)})
summary(sr_intra_sub_mlm[[1]])
summary(sr_intra_sub_mlm[[2]])
summary(sr_intra_sub_mlm[[3]])
summary(sr_intra_sub_mlm[[4]])

sr_intra_sub_mlm <- map(sr_intra_sub, ~ {lmer(BRS.s ~ run.s + ses.s + mean_fd.s + (1 | ses.s),data = .x)})
summary(sr_intra_sub_mlm[[1]])
summary(sr_intra_sub_mlm[[2]])
summary(sr_intra_sub_mlm[[3]])
summary(sr_intra_sub_mlm[[4]])


#Mood associations

##MID
mid_intra_sub_mlm <- map(mid_intra_sub, ~ {lmer(VS.s ~ mood.s + kss.s + panas_pa.s + (1 | ses.s),data = .x)})

summary(mid_intra_sub_mlm[[1]])
summary(mid_intra_sub_mlm[[2]])
summary(mid_intra_sub_mlm[[3]])
summary(mid_intra_sub_mlm[[4]])
partR2(mid_intra_sub_mlm[[1]],data=mid_intra_sub[[1]],partvars <- c("mood.s", "kss.s", "panas_pa.s"))
partR2(mid_intra_sub_mlm[[2]],data=mid_intra_sub[[2]],partvars <- c("mood.s", "kss.s", "panas_pa.s"))
partR2(mid_intra_sub_mlm[[3]],data=mid_intra_sub[[3]],partvars <- c("mood.s", "kss.s", "panas_pa.s"))
partR2(mid_intra_sub_mlm[[4]],data=mid_intra_sub[[4]],partvars <- c("mood.s", "kss.s", "panas_pa.s"))


mid_intra_sub_mlm <- map(mid_intra_sub, ~ {lmer(BRS.s ~ mood.s + kss.s + panas_pa.s + (1 | ses.s),data = .x)})

summary(mid_intra_sub_mlm[[1]])
summary(mid_intra_sub_mlm[[2]])
summary(mid_intra_sub_mlm[[3]])
summary(mid_intra_sub_mlm[[4]])
partR2(mid_intra_sub_mlm[[1]],data=mid_intra_sub[[1]],partvars <- c("mood.s", "kss.s", "panas_pa.s"))
partR2(mid_intra_sub_mlm[[2]],data=mid_intra_sub[[2]],partvars <- c("mood.s", "kss.s", "panas_pa.s"))
partR2(mid_intra_sub_mlm[[3]],data=mid_intra_sub[[3]],partvars <- c("mood.s", "kss.s", "panas_pa.s"))
partR2(mid_intra_sub_mlm[[4]],data=mid_intra_sub[[4]],partvars <- c("mood.s", "kss.s", "panas_pa.s"))


##sr
sr_intra_sub_mlm <- map(sr_intra_sub, ~ {lmer(VS.s ~ mood.s + kss.s + panas_pa.s + (1 | ses.s),data = .x)})

summary(sr_intra_sub_mlm[[1]])
summary(sr_intra_sub_mlm[[2]])
summary(sr_intra_sub_mlm[[3]])
summary(sr_intra_sub_mlm[[4]])
partR2(sr_intra_sub_mlm[[1]],data=mid_intra_sub[[1]],partvars <- c("mood.s", "panas_pa.s", "kss.s"))
partR2(sr_intra_sub_mlm[[2]],data=mid_intra_sub[[2]],partvars <- c("mood.s", "panas_pa.s", "kss.s"))
partR2(sr_intra_sub_mlm[[3]],data=mid_intra_sub[[3]],partvars <- c("mood.s", "panas_pa.s", "kss.s"))
partR2(sr_intra_sub_mlm[[4]],data=mid_intra_sub[[4]],partvars <- c("mood.s", "panas_pa.s", "kss.s"))


sr_intra_sub_mlm <- map(sr_intra_sub, ~ {lmer(BRS.s ~ mood.s + kss.s + panas_pa.s + (1 | ses.s),data = .x)})

summary(sr_intra_sub_mlm[[1]])
summary(sr_intra_sub_mlm[[2]])
summary(sr_intra_sub_mlm[[3]])
summary(sr_intra_sub_mlm[[4]])
partR2(sr_intra_sub_mlm[[1]],data=mid_intra_sub[[1]],partvars <- c("mood.s", "panas_pa.s", "kss.s"))
partR2(sr_intra_sub_mlm[[2]],data=mid_intra_sub[[2]],partvars <- c("mood.s", "panas_pa.s", "kss.s"))
partR2(sr_intra_sub_mlm[[3]],data=mid_intra_sub[[3]],partvars <- c("mood.s", "panas_pa.s", "kss.s"))
partR2(sr_intra_sub_mlm[[4]],data=mid_intra_sub[[4]],partvars <- c("mood.s", "panas_pa.s", "kss.s"))
