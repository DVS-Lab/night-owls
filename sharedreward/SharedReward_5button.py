###RF-1 sequence pilot
###shared reward, ER design
###Dominic Fareri

from psychopy import visual, core, event, gui, data, sound, logging
import csv
import datetime
import random
import numpy
import os
import pyglet
import pandas as pd

#parameters
useFullScreen = True
useDualScreen=2
DEBUG = False

frame_rate=1
instruct_dur=3
initial_fixation_dur = 4
#final_fixation_dur = 8
decision_dur=2.5
outcome_dur=1

responseKeys=('3','2','z')

#get subjID
subjDlg=gui.Dlg(title="Shared Reward Task")
subjDlg.addField('Participant:')
#subjDlg.addField('Enter Friend Name: ') #1
#subjDlg.addField('Enter Partner Name: ')#NOTE: PARTNER IS THE CONFEDERATE/STRANGER #2
subjDlg.addField('Session:', choices=['1', '2','3','4','5','6','7','8','9','10','11','12'])
subjDlg.addField('Run:', choices=['1', '2'])
#subjDlg.addField('MB:', choices=['2', '3'])
#subjDlg.addField('ME:', choices=['1', '3', '4'])
#subjDlg.addField('IP:', choices=['0', '2'])
subjDlg.show()

if subjDlg.show():  # This displays the dialog
    subj_id = subjDlg.data[0]  # Extract the Subject ID
    ses = subjDlg.data[1]  # Extract the selected session
    run = subjDlg.data[2]  # Extract the selected Run
    #mb = subjDlg.data[3]  # Extract the MB value
    #me = subjDlg.data[4]  # Extract the ME value
else:
    core.quit()  # Gracefully exit if "Cancel" is pressed

run_data = {
    'Participant ID': subj_id,
    'Date': str(datetime.datetime.now()),
    'Description': 'NOSC - SharedReward Task',
    'Session': ses,
    'Run': run,
    #'MB': mb,
    #'ME': me
    }

#window setup
display = pyglet.canvas.get_display()
screens = display.get_screens()
win = visual.Window([800,600], monitor="testMonitor", units="deg", fullscr=useFullScreen, allowGUI=False, screen=len(screens)-2)

#checkpoint
print("got to check 1")

#define stimulus
fixation = visual.TextStim(win, text="+", height=5,color="black")

#create ellipsis stimulus
dots = visual.TextStim(win, 
    pos=[0, 0], 
    text='...', 
    height=5, 
    color='black')
DotsClock = core.Clock()

#waiting for trigger
ready_screen = visual.TextStim(win, text="Please wait for the block of trials to begin. \n\nRemember to make your choice when the question mark is on the screen and keep your head still!\n You do NOT need to be as quick as possible, but don't be too slow.", height=2.5,wrapWidth=40)

#decision screen
nameStim = visual.TextStim(win=win,font='Arial',pos=(0, 4.5), height=2, color='white', colorSpace='rgb', opacity=1,depth=-1.0);
cardStim = visual.Rect(win=win, name='polygon', width=8, height=14, ori=0, pos=(0, -4),lineWidth=5, lineColor=[1,1,1], lineColorSpace='rgb',fillColor=[0,0,0], fillColorSpace='rgb',opacity=1, depth=0.0, interpolate=True)
question = visual.TextStim(win=win, name='text',text='?',font='Arial',pos=(0, -4), height=4, wrapWidth=None, ori=0, color='white', colorSpace='rgb', opacity=1,depth=-1.0);
pictureStim =  visual.ImageStim(win, pos=(0,10.5), size=(10,10))

#outcome screen
outcome_cardStim = visual.Rect(win=win, name='polygon', width=8, height=14, ori=0, pos=(0, -4),lineWidth=5, lineColor=[1,1,1], lineColorSpace='rgb',fillColor=[0,0,0], fillColorSpace='rgb',opacity=1, depth=0.0, interpolate=True)
outcome_text = visual.TextStim(win=win, name='text',text='',font='Arial',pos=(0, -4), height=4, wrapWidth=None, ori=0, color='white', colorSpace='rgb', opacity=1,depth=-1.0);
outcome_money = visual.TextStim(win=win, name='text',text='',font='Wingdings 3',pos=(0, 1), height=4.5, wrapWidth=None, ori=0, colorSpace='rgb', opacity=1,depth=-1.0);

#instructions
instruct_screen = visual.TextStim(win, text='In this game you will have to guess the numerical value of a card for a chance to win some money.\nIf you think the value of the card will be lower than 5, press with your middle finger.\n\nIf you think the value of the card will be higher than 5, press with your index finger.', pos = (0,1), wrapWidth=40, height = 2.5)
instruct_screen2 = visual.TextStim(win, text='Remember, you will be sharing monetary outcomes on each trial with the partner displayed at the top of the screen––either the computer or a previous participant.\n\nIf you guess correctly, you and your partner earn $10 ($5 each).\n If you guess incorrectly, you and your partner lose $5 ($2.50 each).', pos = (0,1), wrapWidth=40, height = 2.5)

#exit
exit_screen = visual.TextStim(win, text='Thanks for playing! Please wait for instructions from the experimenter.', pos = (0,1), wrapWidth=20, height = 2.5)

#logging
expdir = os.getcwd()
subjdir = '%s/logs/sub-%s' % (expdir, subj_id)
if not os.path.exists(subjdir):
    os.makedirs(subjdir)
log_file = os.path.join(f'sub-{subj_id}_task-sharedreward_ses-{ses}_run-{run}_raw.csv')

globalClock = core.Clock()
logging.setDefaultClock(globalClock)

timer = core.Clock()

#trial handler
trial_data = [r for r in csv.DictReader(open('%s/event-related/params/sub-' % (os.getcwd()) + subj_id + '/sub-'
    + subj_id + '_ses-' + ses + '_run-' + run + '_design.csv','rU'))]
outfile = f"logs/sub-{subj_id}/sub-{subj_id}_task-sharedreward_ses-{ses}_run-{run}_raw.csv"

#trial_data = [r for r in csv.DictReader(open('SharedReward_design.csv','rU'))]
#trials = data.TrialHandler(trial_data[:], 1, method="sequential") #change to [] for full run

trials_run = data.TrialHandler(trial_data[:], 1, method="sequential") #change to [] for full run

#set partner names
# 3 = friend, 2 = confederate, 1 = computer
# change names accordingly here

stim_map = {
  '2': 'Jack',
  '1': 'Computer',
  }

image_map = {
  '2': 'stranger',
  '1': 'computer',
}

outcome_map = {
  '3': 'reward',
  '2': 'neutral',
  '1': 'punish',
  }

#checkpoint
print("got to check 2")

'''
runs=[]
for run in range(1):
    run_data = []
    for t in range(8):
        sample = random.sample(range(len(trial_data)),1)[0]
        run_data.append(trial_data.pop(sample))
    runs.append(run_data)
'''

# Instructions
instruct_screen.draw()
win.flip()
event.waitKeys(keyList=('2'))

instruct_screen2.draw()
win.flip()
event.waitKeys(keyList=('2'))

# main task loop
def do_run(run, trials):
    resp=[]
    fileName=log_file.format(subj_id,ses,run)



    #wait for trigger
    ready_screen.draw()
    win.flip()
    event.waitKeys(keyList=('equal'))
    globalClock.reset()
    studyStart = globalClock.getTime()
    trials.addData('studyStart',studyStart)

    #Initial Fixation screen
    dots.draw()
    win.flip()
    initial_fixation_Onset = globalClock.getTime()
    trials.addData('InitFixOnset',initial_fixation_Onset)
    core.wait(initial_fixation_dur)
    initial_fixation_offset = globalClock.getTime()
    trials.addData('InitFixOffset',initial_fixation_offset)

    #globalClock.reset()

    for trial in trials:
        condition_label = stim_map[trial['Partner']]
        image_label = image_map[trial['Partner']]
        imagepath = os.path.join(expdir,'Images')
        image = os.path.join(imagepath, "%s.png") % image_label
        nameStim.setText(condition_label)
        pictureStim.setImage(image)

        #decision phase
        timer.reset()
        event.clearEvents()

        resp=[]
        resp_val=None
        resp_onset=None

        decision_onset = globalClock.getTime()
        trials.addData('decision_onset',decision_onset)

        while timer.getTime() < decision_dur:
            cardStim.draw()
            question.draw()
            pictureStim.draw()
            nameStim.draw()
            #dec_screen_pre = globalClock.getTime()
            win.flip()
            #dec_screen_post = globalClock.getTime()

            resp = event.getKeys(keyList = responseKeys)
            #trials.addData('dec_screen_pre',dec_screen_pre)
            #trials.addData('dec_screen_post',dec_screen_post)

            if len(resp)>0:
                if resp[0] == 'z':
                #trials.saveAsText(fileName=log_file.format(subj_id),delim=',',dataOut='all_raw')
                    os.chdir(subjdir)
                    trials.saveAsWideText(fileName)
                    os.chdir(expdir)
                    win.close()
                    core.quit()
                resp_val = int(resp[0])
                if resp_val==2:
                    #resp_onset = globalClock.getTime()
                    question.setColor('darkorange')
                    #rt = resp_onset - decision_onset
                    #core.wait(decision_dur - rt)
                if resp_val==3:
                    #resp_onset = globalClock.getTime()
                    question.setColor('darkorange')
                    #rt = resp_onset - decision_onset
                    #core.wait(decision_dur -rt)
                cardStim.draw()
                question.draw()
                pictureStim.draw()
                nameStim.draw()
                #orangeText_pre = globalClock.getTime()
                win.flip()
                #orangeText_post = globalClock.getTime()
                resp_onset = globalClock.getTime()
                rt = resp_onset - decision_onset
                core.wait(0.1)
                #core.wait(decision_dur - rt)
                break
            else:
                resp_val = 0
                #resp_onset = 999
                resp_onset = globalClock.getTime()
                #rt = 0
                rt = resp_onset - decision_onset
                core.wait(0.1)

        trials.addData('resp', int(resp_val))
        trials.addData('resp_onset', resp_onset)
        trials.addData('rt', rt)
        #trials.addData('orangeText_pre',orangeText_pre)
        #trials.addData('orangeText_post',orangeText_post)

        ###reset question mark color
        question.setColor('white')

        #ISI

        timer.reset()

        given_ISI = float(trial['ISI'])
        isi_for_trial = float(given_ISI+(2.5-(rt+0.1)))
        ISI_onset = globalClock.getTime()
        trials.addData('ISI_onset', ISI_onset)
        trials.addData('isi_for_trial', isi_for_trial)

        fixation.draw()
        #ISI_pre_onset = globalClock.getTime()
        win.flip()
        #ISI_post_onset = globalClock.getTime()
        core.wait(isi_for_trial)

        ISI_offset = globalClock.getTime()
        trials.addData('ISI_offset', ISI_offset)
        #trials.addData('ISI_pre_onset',ISI_pre_onset)
        #trials.addData('ISI_post_onset',ISI_post_onset)
        #outcome phase
        timer.reset()
        #win.flip()
        outcome_onset = globalClock.getTime()

        while timer.getTime() < outcome_dur:
            outcome_cardStim.draw()
            #pictureStim.draw()
            #nameStim.draw()
            #win.flip()

            if trial['Feedback'] == '3' and resp_val == 3:
                outcome_txt = int(random.randint(1,4))
                outcome_moneyTxt= 'h'
                outcome_color='lime'
                trials.addData('outcome_val', int(outcome_txt))
            elif trial['Feedback'] == '3' and resp_val == 2:
                outcome_txt = int(random.randint(6,9))
                outcome_moneyTxt= 'h'
                outcome_color='lime'
                trials.addData('outcome_val', int(outcome_txt))
            elif trial['Feedback'] == '2' and resp_val == 3:
                outcome_txt = int(5)
                outcome_moneyTxt= 'n'
                outcome_color='white'
                trials.addData('outcome_val', int(outcome_txt))
            elif trial['Feedback'] == '2' and resp_val == 2:
                outcome_txt = int(5)
                outcome_moneyTxt= 'n'
                outcome_color='white'
                trials.addData('outcome_val', int(outcome_txt))
            elif trial['Feedback'] == '1' and resp_val == 3:
                outcome_txt = int(random.randint(6,9))
                outcome_moneyTxt= 'i'
                outcome_color='darkred'
                trials.addData('outcome_val', int(outcome_txt))
            elif trial['Feedback'] == '1' and resp_val == 2:
                outcome_txt = int (random.randint(1,4))
                outcome_moneyTxt= 'i'
                outcome_color='darkred'
                trials.addData('outcome_val', int(outcome_txt))
            elif resp_val == 0:
                outcome_txt='#'
                outcome_moneyTxt = ''
                outcome_color='white'
                outcome_value='999'
                trials.addData('outcome_val',int(outcome_value))


            outcome_text.setText(outcome_txt)
            outcome_money.setText(outcome_moneyTxt)
            outcome_money.setColor(outcome_color)
            outcome_text.draw()
            outcome_money.draw()
            #outcome_pre_onset = globalClock.getTime()
            win.flip()
            #outcome_post_onset = globalClock.getTime()
            core.wait(outcome_dur)
            #trials.addData('outcome_val', outcome_txt)
            trials.addData('outcome_onset', outcome_onset)
            #trials.addData('outcome_pre_onset', outcome_pre_onset)
            #trials.addData('outcome_post_onset', outcome_post_onset)

            outcome_offset = globalClock.getTime()
            trials.addData('outcome_offset', outcome_offset)

            duration = outcome_offset - decision_onset
            trials.addData('trialDuration', duration)

            event.clearEvents()
        print("got to check 3")


        #ITI
        logging.log(level=logging.DATA, msg='ITI') #send fixation log event
        timer.reset()
        ITI_onset = globalClock.getTime()
        iti_for_trial = float(trial['ITI'])
        #while timer.getTime() < iti_for_trial:
        dots.draw()
        #ITI_pre_onset = globalClock.getTime()
        win.flip()
        #ITI_post_onset = globalClock.getTime()
        core.wait(iti_for_trial)
        ITI_offset = globalClock.getTime()
        trials.addData('ITIonset', ITI_onset)
        trials.addData('ITIoffset', ITI_offset)
        #trials.addData('ITI_pre_onset',ITI_pre_onset)
        #trials.addData('ITI_post_onset',ITI_post_onset)

    #Final Fixation screen after trials completed
    dots.draw()
    win.flip()
    expected_dur = 428
    #buffer_dur = 4
    #total_dur = expected_dur + buffer_dur
    #if globalClock.getTime() < expected_dur:
    endTime = (expected_dur - globalClock.getTime())
   # else:
    #    endTime = buffer_dur
    core.wait(endTime)
    final_fixation_offset = globalClock.getTime()
    trials.addData('final_fix_offset', final_fixation_offset)
    
    os.chdir(subjdir)
    trials.saveAsWideText(fileName)
    os.chdir(expdir)    


for run, trials in enumerate([trials_run]):
    do_run(run, trials)

# Exit
# Load the CSV file

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

exit_screen = visual.TextStim(win, text=f"Thanks for playing, you earned ${total_earnings} for this run!\n\nPlease wait for instructions from the experimenter.", pos = (0,1), wrapWidth=40, height = 2.5)
exit_screen.draw()
win.flip()
event.waitKeys()

# End the experiment
core.quit()
