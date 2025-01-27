set.seed(125) # Set a seed for reproducibility


colors <- c('Green','Blue')
csequence <- sample(colors, size = 45, replace = TRUE)
print(csequence)
table(csequence) / length(csequence)

# Define the numbers and their probabilities
#numbers <- c(1, 2, 4, 7)
numbers <- c(1,1,1,1)
probabilities <- c(0.50, 0.25, 0.15, 0.10)

# Generate the sequence of 45 numbers
sequence <- sample(numbers, size = 45, replace = TRUE, prob = probabilities)

# Display the sequence
print(sequence)

# Check the proportion of each number
table(sequence) / length(sequence)


mid_trials <- as.data.frame(cbind(csequence,sequence))
colnames(mid_trials) <- c('CueColor','itiTime')
write.csv(mid_trials,'C:/Users/mmatt/Desktop/Projects/NightOwls/night-owls/mid/MID_trials.csv',row.names = F)
