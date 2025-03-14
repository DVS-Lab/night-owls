set.seed(1125) # Set a seed for reproducibility


#Practice
ntrials = 12

colors <- c(rep('Green', ntrials/2), rep('Blue', ntrials/2))
csequence <- sample(colors) 
print(csequence)

itis <- c(rep(1, ntrials/2), rep(2, ntrials/4), rep(3, ntrials/6), rep(4, ntrials/12))
isequence <- sample(itis)
print(isequence)

mid_trials <- as.data.frame(cbind(csequence,isequence))
colnames(mid_trials) <- c('CueColor','itiTime')
write.csv(mid_trials,'C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/mid/timing/prac_MID_trials.csv',row.names = F)


sum(itis) + sum(isis)

set.seed(125)
#Full
ntrials = 56
colors <- c(rep('Green', ntrials/2), rep('Blue', ntrials/2))
isis <- c(rep(1.5, ntrials/2), rep(2, ntrials/4),  rep(2.5, ntrials/8), rep(3, ntrials/8))
itis <- c(rep(2, ntrials/2), rep(3, ntrials/4), rep(5, ntrials/8), rep(7, ntrials/8))

for (ses in 1:12){
  for (run in 1:2){
    csequence <- sample(colors) 
    isisequence <- sample(isis)
    itisequence <- sample(itis)
    mid_trials <- as.data.frame(cbind(csequence,isisequence,itisequence))
    mid_trials[4] <- as.numeric(mid_trials[,3]) - (1.5-as.numeric(mid_trials[,2]))
    mid_trials <- mid_trials[-3]
    colnames(mid_trials) <- c('CueColor','isiTime','itiTime')
    print(sum(as.numeric(mid_trials[,2]))+sum(mid_trials[3]))
    write.csv(mid_trials,paste0('C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/mid/timing/ses-',ses,'_run-',run,'_MID_trials.csv'),row.names = F)
  }
}

