import pandas as pd

#MID
total_earnings_across_sessions = 0

# Iterate over sessions
for ses in range(6, 13):
    total_earnings_for_session = 0
    
    # Process MID runs
    for run in range(1, 3):
        outfile = f'C:\\Users\\mmatt\\Desktop\\Projects\\NightOwls\\night-owls\\stimuli\\mid\\data\\sub-105\\sub-105_task-mid_ses-{ses}_run-{run}.csv'
        #outfile = f'C:\\Users\\Public\\LAB PROJECTS\\Smith-Lab\\GitHub\\night-owls\\stimuli\\mid\\data\\sub-104\\sub-104_task-mid_ses-{ses}_run-{run}.csv'
        df = pd.read_csv(outfile)

        # Randomly select 10 trials
        selected_trials = df.sample(n=10, random_state=125)

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

        if total_earnings > 20:
            total_earnings = 20
        
        total_earnings_for_session += total_earnings
        print(f'Earnings for MID run {run} in session {ses:02d}: {total_earnings}')  # Print earnings for each run

    # Process Shared Reward runs
    for run in range(1, 3):
        outfile = f'C:\\Users\\mmatt\\Desktop\\Projects\\NightOwls\\night-owls\\stimuli\\sharedreward\\logs\\sub-105\\sub-105_task-sharedreward_ses-{ses}_run-{run}_raw.csv'  # Correct file path for Shared Reward
        #outfile = f'C:\\Users\\Public\\LAB PROJECTS\\Smith-Lab\\GitHub\\night-owls\\stimuli\\sharedreward\\logs\\sub-104\\sub-104_task-sharedreward_ses-{ses}_run-{run}_raw.csv'  # Correct file path for Shared Reward

        df = pd.read_csv(outfile)

        # Randomly select 10 trials
        selected_trials = df.sample(n=10, random_state=125)

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

        total_earnings_for_session += total_earnings
        print(f'Earnings for Shared Reward run {run} in session {ses:02d}: {total_earnings}')  # Print earnings for each run

    # Check if total earnings for session exceeds 50
    if total_earnings_for_session > 50:
        print(f'Sum for session {ses:02d}: {total_earnings_for_session} (set to 50 for sum)')  # Print sum within session
        total_earnings_for_session = 50  # Set to 50 if it exceeds
    else:
        print(f'Sum for session {ses:02d}: {total_earnings_for_session}')  # Print sum within session

    total_earnings_across_sessions += total_earnings_for_session

print(f'Total earnings across all sessions: {total_earnings_across_sessions}')  # Print total sum

#Sub-101 ses 1-4: $167
#Sub-101 ses 5-9: 

#Sub-105 ses 1-5: 161.5
#Sub-105 ses 6-12: 276


#Sub-104 ses- 1-4: $170.50
#5-10: 225

#Sub-105: Ses 1-5: 
#Ses 6-..