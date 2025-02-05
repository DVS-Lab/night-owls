# -*- coding: utf-8 -*-
"""
MID.py

Monetary incentive delay task. Participants attend a circular cue, and are required to respond to a subsequently presented white rectangle while it is presented.
Rectanhle presentation times vary based on a stepwise procedure calibrated to reach 63% performance. 
Circular cues can be of one of two types predicting two possible outcomes - reward (green cues) or no reward (blue cues). 
Responding in time for a reward cues yields a monetary gain; responding in time to a no-reward cue does nothing.

Current version: 45 trials, 30 potentially yielding $1 reward per trial, 15 no reward (no punishment trials). 
WRITTEN FOR fMRI - affects the beginning and the end of the task.

Originally written for PsychoPy v 1.84.2
(Peirce, JW (2007) PsychoPy - Psychophysics software in Python. J Neurosci Methods, 162(1-2):8-13)
Oct 13 2020 - Updated for compatibility with Python3 (PsychoPy v 2020.2.4

External dependencies:ses
    a folder (inst_dir parameter) with a txt file containing the instructions (inst_file parameter)
    a csv file (trials_file parameter) in the root directory (the directory running the code) with N+1 lines, where N is the number of trials. First line in the CSV file should be "CueColor".

User-defined variables: see # setting up some user-defined variables

Last updated: October 13 2020

@author: nivreggev (reggevn@bgu.ac.il)
"""
from __future__ import division
from psychopy import gui, visual, core, data, event, logging, monitors
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)
from psychopy.hardware.emulator import launchScan
import numpy as np  # whole numpy lib is available, prepend 'np.'
import pandas as pd
from numpy.random import random, shuffle
import random
import pyglet
import csv
import os  # handy system and path functions
import sys  # to get file system encoding

## setting up some user-defined variables

expName = "MID"
version = "1.0" #
data_dir = "data" # location of outputs to be generated; includes data for participants as well as trial selection and trial presentation sequence
inst_dir = "text" # location of instructions directory
inst_file = ["instructions_MID.txt"] # name of instructions files (needs to be .txt)
study_times = [0.5, 2.25, 0.5, 1, 0.5] # component duration (s): cue, delay, target, feedback, lastfixation
initial_fix_dur = 8 # added time to make sure homogenicity of magnetic field is reached
closing_fix_dur = 10 # added time to make sure haemodynamic responses of the last trials are properly modeled 
min_target_dur = 0.13 # sets the minimum presentation time for target (in seconds)

# settings for fMRI emulation:
MR_settings = {
    'TR': 2.000,     # duration (sec) per whole-brain volume
    'volumes': 110,    # number of whole-brain 3D volumes per scanning run
    'sync': 'equal', # character to use as the sync timing event; assumed to come at start of a volume
    'skip': 2,       # number of volumes lacking a sync pulse at start of scan (for T1 stabilization)
    'sound': True    # in test mode: play a tone as a reminder of scanner noise
    }

## defining some initialization functions

def initialization(expName,version):
    """Present initial dialog; initialize some parameters"""
    # Store info about the experiment session
    expName = expName + version  
    expInfo = {u'participant': '', u'session':['1', '2','3','4','5','6','7','8','9','10','11','12'],u'run':['1','2']}
    dlg = gui.DlgFromDict(dictionary=expInfo, title=expName)
    if dlg.OK == False:
        core.quit()  # user pressed cancel
    expInfo['date'] = data.getDateStr()  # add a simple timestamp
    expInfo['expName'] = expName
    sn = int(expInfo['participant'])
    ses = expInfo['session']
    run = expInfo['run']
    
      
    return(expInfo,expName,sn,ses,run)


def make_screen():
    """Generates screen variables"""
    platform = pyglet.canvas.get_display()
    display = pyglet.canvas.get_display()
    screens = display.get_screens()
    win_res = [screens[-1].width, screens[-1].height]
    exp_mon = monitors.Monitor('exp_mon')
    exp_mon.setSizePix(win_res)
    win = visual.Window(size=win_res, screen=len(screens)-1, allowGUI=True,
                        fullscr=True, monitor=exp_mon, units='height',
                        color=(0.2, 0.2, 0.2))
    return(win_res, win)

def start_datafiles(_thisDir, expName, expInfo, data_dir, sn, ses, run):

        # Create the subject folder path
    subject_folder = os.path.join(_thisDir, data_dir, f"sub-{sn}")
    
    # Check if the folder exists, if not, create it
    if not os.path.exists(subject_folder):
        os.makedirs(subject_folder)

    fname = f"mid_sub-{sn}_ses-{ses}_run-{run}"
    curdirlist = os.listdir(_thisDir + os.sep + data_dir)
    for i in curdirlist:
            if i == fname + '.csv':
                warndlg = gui.Dlg(title='Warning!')
                warndlg.addText('A data file with this number already exists.')
                warndlg.addField('Overwrite?\t\t', initial="no")
                warndlg.addField('If no, new SN:\t', initial='0')
                warndlg.show()
                if gui.OK:
                    over = warndlg.data[0].lower() == 'no'
                else:
                    core.quit()
                if over:
                    sn = int(warndlg.data[1])
                    pad = 4-len(str(sn))
                    snstr = '0'*pad + str(sn)
                    fname=f"mid_sub-{sn}_ses-{ses}"
    filename = os.path.join(subject_folder, fname)
    return(filename)
    
def display_inst(instr_part,task,forwardKey,backKey,startKey,instructFinish):
    """ display instructions for a specific experimental task; input includes: 
    instr_part: instructions extracted from text
    task: task serial number (in actual serial order, starting at 1; convetred to Python's representation, where 1 is 0, in the function"""
    endOfInstructions = False
    instructLine = 0
    inst = instr_part[task-1]
    while not endOfInstructions:
        instructPrompt.setText(inst[instructLine])
        instructPrompt.draw()
        if instructLine == 0:
            instructFirst.draw()
            win.flip()
            instructRep = event.waitKeys(keyList=[forwardKey])
        else:
            instructMove.draw()
            win.flip()
            instructRep = event.waitKeys(keyList=[forwardKey, backKey])
        if event.getKeys(keyList=[endKey]):
            core.quit()       
        if instructRep[0] == backKey:
            instructLine -= 1
        elif instructRep[0] == forwardKey:
            instructLine += 1
        if inst[instructLine] == "end":
            endOfInstructions = True
        # check for quit:
    instructFinish.draw()
    win.flip()
    event.waitKeys(keyList=[startKey])
    
### START SET UP OF STUDY 

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))#.decode(sys.getfilesystemencoding())
os.chdir(_thisDir)

# present initialization dialog
[expInfo,expName,sn,ses,run] = initialization(expName,version) 
trials_file = f"/timing/ses-{ses}_run-{run}_MID_trials.csv"


# Data file name creation; later add .psyexp, .csv, .log, etc
filename = start_datafiles(_thisDir, expName, expInfo, data_dir, sn, ses, run)

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version=version, extraInfo=expInfo, runtimeInfo=None,
    originPath=None, savePickle=True, saveWideText=True, dataFileName=filename)

# save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp

# Setup the window and presentation constants
[win_res, win] = make_screen()
yScr = 1.
xScr = float(win_res[0])/win_res[1]
fontH = yScr/25
wrapW = xScr/1.5
textCol = 'black'
# store frame rate of monitor if we can measure it
expInfo['frameRate'] = win.getActualFrameRate()
if expInfo['frameRate'] != None and expInfo['frameRate'] <200: # inserted a manual byapss here to aviod illogical refresh rates
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess

# set random seed - session dependent
random.seed(ses)

# determine accepted inputs 
forwardKey = "2"
backKey = "1"
startKey = "2"
expKeys = ["1","2","3","4","5"] # including all response button keys to catch misaligned fingers/responses
endKey = "5"
# Initialize components for Routine "instructions"
instructFirst = visual.TextStim(win, text="Use your index finger to continue", height=fontH, color=textCol, pos=[0, -yScr/4])
instructMove = visual.TextStim(win, text="Use your index finger to continue, or middle to go back.", height=fontH, color=textCol, pos=[0, -yScr/4])

#import instructions
instr_part = [[],[],[]]
inst_file = inst_file
for inst in range (0,len(inst_file)):
    inname = _thisDir + os.sep + inst_dir + os.sep + inst_file[inst]
    infile = open(inname, 'r')
    for line in infile:
        instr_part[inst].append(line.rstrip())
    instr_part[inst].append("end")
    infile.close()

## START component code to be run before the window creation

#create fixation stimulus
fix = visual.TextStim(win, pos=[0, 0], text='+', height=fontH*2, color=textCol)
FixClock = core.Clock()

#create ellipsis stimulus
dots = visual.TextStim(win, 
    pos=[0, 0], 
    text='...', 
    height=fontH*2, 
    color='black')
DotsClock = core.Clock()

# Initialize components for Routine "instructions"
instructPrompt = visual.TextStim(win=win, font='Arial', pos=(0, yScr/10), height=fontH, wrapWidth=wrapW, color=textCol);
instructFinish = visual.TextStim(win, text="You have reached the end of the instructions. When you are ready to begin the task, place your fingers on the keys notify the experimenter.",
                                     height=fontH, color=textCol, pos=[0, 0], wrapWidth=wrapW)    

# Initialize components for task transitions
wait = visual.TextStim(win, pos=[0, 0], text="The task will begin momentarily. Get ready...", height=fontH, color=textCol)
wait_str = "The task will begin momentarily. Get ready..."
endf = visual.TextStim(win, pos=[0, 0], text="Thank you. This part of the experiment is now complete.",wrapWidth=wrapW, height=fontH, color=textCol)                                     

# Initialize components for Routine "cue" 

Cue = visual.Circle(win,radius=0.2, edges=128,lineWidth=0, pos=(0, 0))
CueClock = core.Clock()
Cue_trials_template = _thisDir + os.sep + trials_file

# Initialize components for Routine "Target"
TargetClock = core.Clock()
Target = visual.Rect(win,width=0.5, height=0.5, fillColor = "white", lineWidth=0, pos=(0,0)) 

# Initialize components for Routine "Feedback"
FeedbackClock = core.Clock()
Trial_FB = visual.TextStim(win=win, name='Trial_FB', text='✓', font='Arial', pos=(0, 0), height=fontH*4, wrapWidth=None, ori=0, 
    color='lime', colorSpace='rgb', opacity=1, bold=True)
Blank_FB_Rectangle = visual.ImageStim(win=win, name='Blank_FB', mask=None,ori=0, pos=(0, 0), size=(xScr/8, xScr/8),texRes=128, interpolate=True)

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

# create content to be displayed
stimuli = pd.read_csv(Cue_trials_template) # read template stimuli
stimuli = stimuli.reindex(np.random.permutation(stimuli.index)) # shuffle order of stimuli

## Displaying Instructions

# keyboard checking is just starting
event.clearEvents(eventType='keyboard')    
event.Mouse(visible=False)  # hide mouse cursor
display_inst(instr_part,1,forwardKey,backKey,startKey,instructFinish)

# reset the non-slip timer for next routine
routineTimer.reset()
event.clearEvents(eventType='keyboard')    
 
### START EXPERIMENTAL LOOP

# wait for TR signal if in scanner
wait.draw()
win.flip()
event.waitKeys(keyList=['equal'])
task_start_time = globalClock.getTime()  # Record the start time of the task
nominalTime = 0 # set up virtual time keeper to align actual with a-priori time allocation

    
# set up counters for trials (to determine cue color and for total earnings
trial_counter = 0
Tot_Earn = 0

# create the staircase handler to adjust for individual threshold (stairs defined in units of screen frames; actual minimum presentation duration is determined by the min_target_dur parameter, the staircase procedure can only add frame rates to that minimum value)
trials = data.StairHandler(startVal=10.0,
    stepType='lin',
    stepSizes=[6, 3, 3, 2, 2, 1, 1],  # reduce step size every two reversals
    minVal=0, maxVal=15,
    nUp=1, nDown=2,  # will home in on the 65% threshold (nUp=1, nDown=3 homes in on 80%)
    nTrials=64,
    extraInfo=expInfo)
print(f"Number of stimuli rows: {len(stimuli)}")

    
thisExp.addLoop(trials)  # add the loop to the experiment
#nominalTime = 0 # set up virtual time keeper to align actual with a-priori time allocation
#globalClock.reset() # to align actual time with virtual time keeper
    
# present initial fixation
t_start = globalClock.getTime()
t = t_start
while t < t_start + initial_fix_dur:
    t = globalClock.getTime()
    fix.draw()
    win.flip()              
nominalTime = t # set up virtual time keeper to align actual with a-priori time allocation

for thisTrial in trials:
    trials.addOtherData('time.onset', globalClock.getTime()) # add trial onset time to the data file
    currentLoop = trials
        
    # update component parameters for each repeat
    Choice_Resp = event.BuilderKeyResponse()
    print('trial counter:',trial_counter)
    CueColor = stimuli.iloc[trial_counter][0] # get cue color from the externally imported stimuli list, based on trial_counter
    trials.addOtherData('cue.color', CueColor) # add trial onset time to the data file

    Cue.fillColor = CueColor
    trial_counter += 1
    
    # ------Prepare to start Routine "Cue"-------
    t = 0
    CueClock.reset()  # clock
    # reset the non-slip timer for next routine
    routineTimer.reset()                
    continueRoutine = True
    routineTimer.add(study_times[0]) # set time limit for current phase
    nominalTime += study_times[0] # update nominal time keeper
    
    # keep track of which components have finished
    CueComponents = [Cue]
    for thisComponent in CueComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    
    # -------Start Routine "Cue"-------
    while continueRoutine and routineTimer.getTime() > 0:
        # get current time
        t = CueClock.getTime()
        
        # first screen updates
        if t >= 0.0 and Cue.status == NOT_STARTED:
            # keep track of start time/frame for later
            Cue.tStart = t
            Cue.setAutoDraw(True)
        frameRemains = 0.0 + study_times[0] - win.monitorFramePeriod * 0.75  # most of one frame period left
        if Cue.status == STARTED and t >= frameRemains:
            Cue.setAutoDraw(False)

        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in CueComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # check for quit (the Esc key)
        if endExpNow or event.getKeys(keyList=[endKey]):
            core.quit()
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "Cue"-------
    for thisComponent in CueComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)

    # ------Prepare to start Routine "ISI"-------
    t = 0
    FixClock.reset()  # clock    
    # reset the non-slip timer for next routine
    routineTimer.reset()   
    continueRoutine = True
    
    # set ISI duration 
    isi_time = stimuli.iloc[trial_counter-1]['isiTime']  # Get ISI from current trial row
    nominalTime += isi_time
    routineTimer.add(isi_time)

    trials.addOtherData('isi_duration', isi_time)   # store total ISI duration


    # keep track of which components have finished
    isiComponents = [fix]
    for thisComponent in isiComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    
    # -------Start Routine "ISI"-------
    while continueRoutine and routineTimer.getTime() > 0:
        t = FixClock.getTime()
        if t >= 0 and fix.status == NOT_STARTED:
            fix.tStart = t
            fix.setAutoDraw(True)
        frameRemains = 0.0 + isi_time - win.monitorFramePeriod * 0.75  
        if fix.status == STARTED and t >= frameRemains:
            fix.setAutoDraw(False)
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in isiComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # check for quit (the Esc key)
        if endExpNow or event.getKeys(keyList=[endKey]):
            core.quit()
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "ISI"-------
    for thisComponent in isiComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)    
    
    # ------Prepare to start Routine "Target"-------
    t = 0
    TargetClock.reset()  # clock
    
    # reset the non-slip timer for next routine
    routineTimer.reset()                   
    continueRoutine = True
    routineTimer.add(study_times[2])
    nominalTime += study_times[2]

    # update component parameters for each repeat
    Target_Resp = event.BuilderKeyResponse()
    
    # keep track of which components have finished
    TargetComponents = [Target, Target_Resp]
    for thisComponent in TargetComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    
    # -------Start Routine "Target"-------
    while continueRoutine and routineTimer.getTime() > 0:
        # get current time
        t = TargetClock.getTime()
        
        # selection screen updates
        if t >= 0.0 and Target.status == NOT_STARTED:
            # keep track of start time/frame for later
            Target.tStart = t
            # display target 
            Target.setAutoDraw(True)
            # open response options
            Target_Resp.tStart = t
            Target_Resp.status = STARTED
            # keyboard checking is just starting
            win.callOnFlip(Target_Resp.clock.reset)  # t=0 on next screen flip
            event.clearEvents(eventType='keyboard')  
            theseKeys = []

        frameRemainsResp = min_target_dur + frameDur*thisTrial #- win.monitorFramePeriod * 0.75  # most of one frame period left;  range: min_target_dur (160ms) to one frame (~17ms, depends on refresh rate) * 6
        if Target.status == STARTED and t >= frameRemainsResp:
            print('thisTrial:',thisTrial) # print for QA purpose
            print('frameDur:',frameDur) # print for QA purpose
            print('frameRemainsResp:',frameRemainsResp) # print for QA purpose
            
            Target.setAutoDraw(False)
            theseKeys = event.getKeys(keyList=expKeys)
            ThisResp = 0 # set response to no response - change only if response was given in the allowed time frame
            
            if len(theseKeys) > 0:  # at least one key was pressed
                ThisResp = 1
                Target_Resp.rt = Target_Resp.clock.getTime()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in TargetComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # check for quit (the Esc key)
        if endExpNow or event.getKeys(keyList=[endKey]):
            core.quit()
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "Target"-------
    for thisComponent in TargetComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    
    # add the data to the staircase so it can be used to calculate the next level
    trials.addResponse(ThisResp)

    # check responses to add RT
    if ThisResp:  # we had a response
        trials.addOtherData('Target_Resp.rt', Target_Resp.rt)
    
    # ------Prepare to start Routine "Feedback"-------
    t = 0
    FeedbackClock.reset()  # clock
    # reset the non-slip timer for next routine
    routineTimer.reset()                   
    continueRoutine = True
    routineTimer.add(study_times[3])
    nominalTime += study_times[3]
    
    # update trial components
    if ThisResp and CueColor=='Green': # if it was a rewarded trial and a response was given
        Tot_Earn += 1
        Trial_FB.setText('✓')  # Show checkmark
        Trial_FB.setColor('lime')  # Bright green color
        Trial_FB.setHeight(fontH*4)  # Make check larger
        trials.addOtherData('Trial.rewardType', '1')
    else:    
        Trial_FB.setText('−')  # Show horizontal bar
        Trial_FB.setColor('#404040')  # Darker grey color
        Trial_FB.setHeight(fontH*8)  # Made bar much larger
        trials.addOtherData('Trial.rewardType', '-1')

    # add to be presented stimuli to output
    trials.addOtherData('Total_Earned', Tot_Earn) 
            
    # keep track of which components have finished
    FeedbackComponents = [Trial_FB]
    for thisComponent in FeedbackComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED

    # -------Start Routine "Feedback"-------
    while continueRoutine and routineTimer.getTime() > 0:
        # get current time
        t = FeedbackClock.getTime()
        
        # feedback screen updates
        if t >= 0.0 and Trial_FB.status == NOT_STARTED:
            # keep track of start time/frame for later
            Trial_FB.tStart = t
            Trial_FB.setAutoDraw(True)
        frameRemains = 0.0 + study_times[3] - win.monitorFramePeriod * 0.75  # most of one frame period left
        if Trial_FB.status == STARTED and t >= frameRemains:
            Trial_FB.setAutoDraw(False)

        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in FeedbackComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # check for quit (the Esc key)
        if endExpNow or event.getKeys(keyList=[endKey]):
            core.quit()
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "Feedback"-------
    for thisComponent in FeedbackComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)

    # ------Prepare to start Routine "ITI"-------
    t = 0
    FixClock.reset()  # clock    
    DotsClock.reset()  # clock    

    # reset the non-slip timer for next routine
    routineTimer.reset()   
    continueRoutine = True
    
    # set ITI duration 
    iti_time = stimuli.iloc[trial_counter-1]['itiTime']  # Get ITI from current trial row
    nominalTime += iti_time
    routineTimer.add(iti_time)
    
    trials.addOtherData('iti_final', iti_time)  # store final adjusted ITI
    
    # keep track of which components have finished
    itiComponents = [dots]
    for thisComponent in itiComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    
    # -------Start Routine "ITI"-------
    while continueRoutine and routineTimer.getTime() > 0:
        t = DotsClock.getTime()
       # if t >= 0 and dots.status == NOT_STARTED:
        #    dots.tStart = t
        #    dots.setAutoDraw(True)
        #frameRemains = 0.0 + iti_time - win.monitorFramePeriod * 0.75  
        #if dots.status == STARTED and t >= frameRemains:
        #    dots.setAutoDraw(False)

        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in itiComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # check for quit (the Esc key)
        if endExpNow or event.getKeys(keyList=[endKey]):
            core.quit()
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "ITI"-------
    for thisComponent in itiComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)

    # add data to log file
    trials.addOtherData('time.global', globalClock.getTime())      
    trials.addOtherData('time.nominal', nominalTime)      
    trials.addOtherData('time.trial', CueClock.getTime())
    #trials.addOtherData('time.plannedfixlength', iti_time)
    #trials.addOtherData('time.actualfixlength', globalClock.getTime() - tend)
                
    # advance to next trial/line in logFile
    thisExp.nextEntry()
    
# ------Prepare to start Routine "dots"-------
t = 0
routineTimer.reset()   
continueRoutine = True

# set fixation time duration 
tend = globalClock.getTime()
total_task_duration = tend - task_start_time
dots_add = 530 - total_task_duration
dotsClock = core.Clock()

routineTimer.reset()  # Reset the routine timer to ensure it's used from the correct starting point
routineTimer.add(dots_add) 

# keep track of which components have finished
dotsComponents = [dots]
for thisComponent in dotsComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED

# -------Start Routine "fix"-------
while continueRoutine and routineTimer.getTime() > 0:
    # get current time
    t = dotsClock.getTime()
    
    # fix updates
    if t >= 0 and dots.status == NOT_STARTED:
        # keep track of start time/frame for later
        dots.tStart = t
        dots.setAutoDraw(True)
    frameRemains = 0.0 + dots_add - win.monitorFramePeriod * 0.75  
    if dots.status == STARTED and t >= frameRemains:
        dots.setAutoDraw(False)   

    # Refresh the screen (show fixation)
    if continueRoutine:
        win.flip()        

# end of study message
fend = globalClock.getTime()
final_duration = fend - task_start_time

endf.draw()
win.flip()
event.waitKeys(keyList=['2'])

# these shouldn't be strictly necessary (should auto-save)
thisExp.addData('Total_Task_Duration', final_duration)
thisExp.saveAsWideText(filename+'.csv',fileCollisionMethod = 'overwrite')
thisExp.saveAsPickle(filename, fileCollisionMethod = 'rename')
logging.flush()

# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
