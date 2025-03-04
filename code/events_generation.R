library(readr)

events_gen <- function(sub,ses){
  sub <- as.character(sub)
  ses <- as.character(ses)
  for (r in 1:2){
    
    #Reset dataframes
    mid.events <- sr.events <- as.data.frame(matrix(NA,0,4))
    
    
    #Generate folders, if don't exist
    out.path <- paste0('C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/bids/sub-',
                       sub,'/ses-0',ses,'/func/')
    if (!dir.exists(out.path)) {
      dir.create(out.path, recursive = TRUE)
    }
    
    
    #MID
    
    raw.mid <- read.csv(paste0('C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/stimuli/mid/data/sub-',sub,'/sub-',sub,
                               '_task-mid_ses-',ses,'_run-',r,'.csv'))
    ##Trial Timing
    mid.cue <- cbind(raw.mid$cue_onset_time[-57],raw.mid$cue_duration[-57],rep(NA,56),'n/a')
    mid.isi <- cbind(raw.mid$isi_onset_time[-57],raw.mid$isi_duration_act[-57],rep(NA,56),'n/a')
    mid.target <- cbind(raw.mid$target_onset_time[-57],raw.mid$target_duration[-57],rep('target',56),'n/a')
    mid.target[,4] <- round(raw.mid$Target_Resp.rt[-57],4)
    mid.target[,4][is.na(mid.target[,4])] <- 'n/a'
    mid.feedback <- cbind(raw.mid$feedback_onset_time[-57],raw.mid$feedback_duration[-57],rep(NA,56),'n/a')
    
    ##Trial Types
    for (t in 1:56){
      if (raw.mid$cue.color[t] == 'Green'){
        mid.cue[t,3] <- 'cue_reward'
        mid.isi[t,3] <- 'isi_reward'
        if (raw.mid$.response[t] == 1){
          mid.feedback[t,3] <- 'feedback_positive_reward'}
        else {
          mid.feedback[t,3] <- 'feedback_negative_reward'}
      }
      else {
        mid.cue[t,3] <- 'cue_neutral'
        mid.isi[t,3] <- 'isi_neutral'
        if (raw.mid$.response[t] == 1){
          mid.feedback[t,3] <- 'feedback_positive_neutral'}
        else {
          mid.feedback[t,3] <- 'feedback_negative_neutral'}
      }
    }
    
    ##Merge base events
    mid.events <- rbind(mid.events,mid.cue,mid.isi,mid.target,mid.feedback)
    
    ##Create combined cue regressor
    mid.antic <- mid.cue
    mid.antic[,2] <- as.numeric(mid.cue[,2]) + as.numeric(mid.isi[,2])
    mid.antic[,3][mid.antic[,3]=='cue_reward'] <- 'anticipation_reward'
    mid.antic[,3][mid.antic[,3]=='cue_neutral'] <- 'anticipation_neutral'
    
    ##Merge, Round, Order
    mid.events <- rbind(mid.events,mid.antic)
    colnames(mid.events) <- c('onset','duration','trial_type','response_time')
    mid.events$onset <- round(as.numeric(mid.events$onset),4)
    mid.events$duration <- round(as.numeric(mid.events$duration),4)
    mid.events <- mid.events[order(mid.events$onset),]
    
    ##Export
    mid.out <- paste0('sub-',sub,'_ses-0',ses,'_task-mid_run-',r,'_events.tsv')
    write_delim(mid.events,paste0(out.path,mid.out),na='n/a',delim = "\t")     
    
  } 
}

#Ses-01 COMPLETED
events_gen(101,1)

#Ses-02 Completed
events_gen(101,02)

#Ses-03 Completed
events_gen(101,03)
