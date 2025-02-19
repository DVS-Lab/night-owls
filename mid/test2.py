import pandas as pd

outfile = f"data/sub-999/mid_sub-999_ses-2_run-1.csv"


# Load the CSV file

df = pd.read_csv(outfile)

# Randomly select 10 trials
selected_trials = df.sample(n=10, random_state=None)

total_earnings = 0

# Process each trial
for index, row in selected_trials.iterrows():
    resp = row['.response']
    rcue = row['cue.color']
    
    if rcue == 'Green' and resp == 1:
        change = 3
    else:
        change = 0     
    total_earnings += change
    
print(f"{total_earnings}")