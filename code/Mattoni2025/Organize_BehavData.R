library(tidyverse)


subj <- c('sub-101','sub-103','sub-104','sub-105')


#Read mood rating data
mood.csv <- sapply(subj,function(x) NULL)
for (sub in 1:length(subj)){
  dir <- paste0('C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/stimuli/mood/logs/'
                ,subj[sub])
  #sessions <- list.dirs(dir)
  #for (ses in 1:length(sessions)){
  mood.files <- list.files(dir,recursive = T,full.names = T)
  mood.names <- sub("\\.csv$", "",list.files(dir,recursive = T))
  ses_num <- as.integer(sub("^ses-([0-9]+).*", "\\1", mood.names))
  mood.files <- mood.files[order(ses_num)]
  mood.names <- mood.names[order(ses_num)]
  mood.names <- sub("^ses-[^/]+/", "", mood.names)  
  mood.csv[[sub]] <- setNames(lapply(mood.files, read.csv), mood.names)
}

#Create wide df
mood_wide <- as.data.frame(matrix(NA,length(subj),6*12+1))

obs_seq <- expand.grid(obs = 1:6,ses = 1:12)
obs_seq <- sprintf("ses-%02d_obs-%02d", obs_seq$ses, obs_seq$obs)

mood_wide[1] <- subj
colnames(mood_wide) <- c('subj',obs_seq)

#Fill df
for (sub in seq_along(subj)){
  for (ses in 1:12){
    for (obs in 1:6){
      
      ses_str <- sprintf("ses-%02d", ses)
      obs_str <- sprintf("obs-%02d", obs)
      list_name <- paste0(subj[sub], "_ses-", ses, "_obs-", obs, "_mood")
      col_name <- paste0(ses_str, "_", obs_str)
      
      
      if (list_name %in% names(mood.csv[[sub]])) {
        mood_wide[sub, col_name] <- mood.csv[[sub]][[list_name]][["Response"]][1]
      } else {
        mood_wide[sub, col_name] <- NA}
    }}}


#Convert to long
mood_long <- mood_wide %>%
  pivot_longer(
    cols = -subj,
    names_to = c("ses", "obs"),
    names_pattern = "ses-(\\d+)_obs-(\\d+)",
    values_to = "mood") %>%
  mutate(
    ses = as.integer(ses),
    obs = as.integer(obs))


#Read alertness data
alert <- read.csv('C:/Users/mmatt/Desktop/Projects/NightOwls/NightOwls-AlertnessData_DATA_2025-08-22_1159.csv')
alert$ses <- as.numeric(substr(alert$redcap_event_name,4,5))
alert$subj <- paste0('sub-',alert$night_owl_id)
alert <- alert[1:47,-c(1,2)]   

#merge
behav_long <- mood_long %>%
  left_join(alert, by = c("subj", "ses")) %>%
  mutate(
    kss = case_when(
      obs %in% 1:3 ~ alertness,
      obs %in% 4:6 ~ fmri_alertness
    )
  ) %>%
  select(subj, ses, obs, mood, kss, caffeine)


#Read PANAS data
panas <- read.csv('C:/Users/mmatt/Desktop/Projects/NightOwls/NightOwls-PANASData_DATA_2025-08-22_1153.csv')
panas$ses <- as.numeric(substr(panas$redcap_event_name,4,5))
panas$subj <- paste0('sub-',panas$night_owl_id)
panas <- panas[,c(26,25,3,4)]
colnames(panas) <- c('subj','ses','panas_pa','panas_na')

#merge
behav_long <- merge(behav_long,panas,by=c("subj","ses"),all.x=T)

#Remove 101 ses12
behav_long <- behav_long[-c(67:72),]

              

#Save
write.csv(behav_long,
          'C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/derivatives/data-outputs/behavioral_data.csv',
          row.names = F)
rm(list=setdiff(ls(), "behav_long"))