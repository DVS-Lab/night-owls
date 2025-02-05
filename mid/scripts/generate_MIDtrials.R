set.seed(1125) # Set a seed for reproducibility


#Practice
ntrials = 12

colors <- c('Green','Blue')
csequence <- sample(colors, size = ntrials, replace = TRUE)
print(csequence)
table(csequence) / length(csequence)

# Define the numbers and their probabilities
numbers <- c(1, 2, 4, 7)
#numbers <- c(1,1,1,1)
probabilities <- c(0.50, 0.25, 0.15, 0.10)

# Generate the sequence of numbers
sequence <- sample(numbers, size = ntrials, replace = TRUE, prob = probabilities)

# Display the sequence
print(sequence)

# Check the proportion of each number
table(sequence) / length(sequence)


mid_trials <- as.data.frame(cbind(csequence,sequence))
colnames(mid_trials) <- c('CueColor','itiTime')
write.csv(mid_trials,'C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/mid/timing/prac_MID_trials.csv',row.names = F)




set.seed(125)
#Full
ntrials = 75
for (ses in 1:12){
  for (run in 1:2){
    csequence <- sample(colors, size = ntrials, replace = TRUE)
    print(csequence)
    table(csequence) / length(csequence)
    
    
    sequence <- sample(numbers, size = ntrials, replace = TRUE, prob = probabilities)
    print(sequence)
    table(sequence) / length(sequence)
    
    
    mid_trials <- as.data.frame(cbind(csequence,sequence))
    colnames(mid_trials) <- c('CueColor','itiTime')
    write.csv(mid_trials,paste0('C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/mid/timing/ses-',ses,'_run-',run,'_MID_trials.csv'),row.names = F)
  }
}


