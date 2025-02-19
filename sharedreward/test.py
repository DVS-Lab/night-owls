from psychopy import core, data
import pandas as pd
import random

# Load the CSV file
df = pd.read_csv("logs/sub-999/sub-999_task-sharedreward_ses-3_run-2_raw.csv")

# Randomly select 10 trials
selected_trials = df.sample(n=20, random_state=None)

total_earnings = 0

# Process each trial
for index, row in selected_trials.iterrows():
    resp = row['resp']
    outcome_val = row['outcome_val']
    
    if outcome_val == 5:
        change = 0
    elif resp == 2 and outcome_val > 5:
        change = 5
    elif resp == 2 and outcome_val < 5:
        change = -2.5
    elif resp == 3 and outcome_val > 5:
        change = -2.5
    elif resp == 3 and outcome_val < 5:
        change = 5
    else:
        change = 0
    
    total_earnings += change


# Adjust final sum based on conditions
if total_earnings < 3:
    total_earnings = 4
elif total_earnings > 20:
    total_earnings = 20
    
    
# Print the sum of the 10 trials
print(f"Total earnings from 10 trials: ${total_earnings}")

# End the experiment
core.quit()
