library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(patchwork)

behav_long <- read.csv('C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/derivatives/data-outputs/behavioral_data.csv')

#Change across all obs
##All subj
pos_mood_mlm <- lmer(mood ~ obs + (1 | subj/ses), data = behav_long)
summary(pos_mood_mlm)
##Ind Subj
mood_subj <- split(behav_long, behav_long$subj)
mood_subj_mlm <- map(mood_subj, ~ {lmer(mood ~ obs + (1 | ses), data = .x)})
summary(mood_subj_mlm[[1]])
summary(mood_subj_mlm[[2]])
summary(mood_subj_mlm[[3]])
summary(mood_subj_mlm[[4]])
###No systematic change across obs


#Induction - Obs 3 vs 4 
mood_induct <- behav_long[behav_long$obs==3 | behav_long$obs==4,]
mood_induct$obs <- factor(mood_induct$obs, levels = c(3, 4), labels = c("pre", "post"))
##All subj
induct_mlm <- lmer(mood ~ obs + (1 | subj/ses), data = mood_induct)
summary(induct_mlm)
##Ind Subj
# split that by subject
induct_subj <- split(mood_induct, mood_induct$subj)
induct_subj_mlm <- map(induct_subj, ~ {lmer(mood ~ obs + (1 | ses), data = .x)})
summary(induct_subj_mlm[[1]])
summary(induct_subj_mlm[[2]])
summary(induct_subj_mlm[[3]])
summary(induct_subj_mlm[[4]])


#Induction - Average Pre vs post
##All subj
mood_avg <- behav_long %>%
  mutate(obs = ifelse(obs <= 3, "pre", "post")) %>%
  group_by(subj, ses, obs) %>%
  summarise(mood = mean(mood, na.rm = TRUE), .groups = "drop") %>%
  mutate(obs = factor(obs, levels = c("pre", "post"))) %>% 
  arrange(subj, ses, obs)


induct_avg_mlm <- lmer(mood ~ obs + (1 | subj/ses), data = mood_avg)
summary(induct_avg_mlm)
##Ind Subj
# split that by subject
induct_avg_subj <- split(mood_avg, mood_avg$subj)
induct_avg_subj_mlm <- map(induct_avg_subj, ~ {lmer(mood ~ obs + (1 | ses), data = .x)})
summary(induct_avg_subj_mlm[[1]])
summary(induct_avg_subj_mlm[[2]])
summary(induct_avg_subj_mlm[[3]])
summary(induct_avg_subj_mlm[[4]])


rm(list=setdiff(ls(), c('behav_long','mood_induct','mood_avg')))


#PANAS
## M and SD
behav_long %>%
  group_by(subj, ses) %>% 
  slice(1) %>%
  group_by(subj) %>%
  summarise(panas_pa_sd = mean(panas_pa, na.rm = TRUE))
behav_long %>%
  group_by(subj, ses) %>% 
  slice(1) %>%
  group_by(subj) %>%
  summarise(panas_pa_sd = sd(panas_pa, na.rm = TRUE))


#Alertness differences
alertness_test <- behav_long[behav_long$obs==1 | behav_long$obs==4,]
alertness_test$alert <- 9 - alertness_test$kss

alertness_test$obs <- factor(alertness_test$obs, levels = c(1, 4), labels = c("pre", "post"))
##All subj
alert_mlm <- lmer(alert ~ obs + (1 | subj/ses), data = alertness_test)
summary(alert_mlm)
##Ind Subj
# split that by subject
alert_subj <- split(alertness_test, alertness_test$subj)
alert_subj_mlm <- map(alert_subj, ~ {lmer(alert ~ obs + (1 | ses), data = .x)})
summary(alert_subj_mlm[[1]])
summary(alert_subj_mlm[[2]])
summary(alert_subj_mlm[[3]])
summary(alert_subj_mlm[[4]])




