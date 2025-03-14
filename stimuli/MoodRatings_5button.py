import os
import csv
import re
from psychopy import visual, event, core, gui, prefs, monitors
import pyglet
import sys
import numpy

useFullScreen = True

useDualScreen=2

#get subjID
subjDlg=gui.Dlg(title="Mood Check")
subjDlg.addField('Subject:')
subjDlg.addField('Session:', choices=['1', '2','3','4','5','6','7','8','9','10','11','12'])
subjDlg.addField("Observation", choices=['1','2','5','6'])
subjDlg.addText("Observations:\n1 = Baseline\n2 = Post Task A Run 1\n5 = Post Task A Run 2\n6 = Post Task B Run 2\n(3 and 4 are done with induction task)")
subjDlg.show()

if subjDlg.show():  # This displays the dialog
    subj_id = subjDlg.data[0]  # Extract the Subject ID
    ses = subjDlg.data[1]  # Extract the selected session
    obs = subjDlg.data[2]  # Extract the selected observation
else:
    core.quit()  # Gracefully exit if "Cancel" is pressed
    
    
def make_screen():
    """Generates screen variables"""
    platform = pyglet.canvas.get_display()
    display = pyglet.canvas.get_display()
    screens = display.get_screens()
    win_res = [screens[-2].width, screens[-2].height]
    exp_mon = monitors.Monitor('testMonitor')
    exp_mon.setSizePix(win_res)
    win = visual.Window(size=win_res, screen=useDualScreen, allowGUI=True,
                        fullscr=True, monitor='testMonitor', units='height',
                        color=(0.2, 0.2, 0.2))
    return(win_res, win)
  
[win_res, win] = make_screen()
  
xScr = float(win_res[0])/win_res[1]
yScr = 1.
fontH = yScr/25
wrapW = xScr/1.5

# Folder structure
folder_path = f"mood/logs/sub-{subj_id}/ses-{ses}"

# Check if the folder exists, and create it if not
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

# Prepare the file name based on the inputs
file_name = f"{folder_path}/sub-{subj_id}_ses-{ses}_obs-{obs}_mood.csv"

# Initialize the window
display = pyglet.canvas.get_display()
screens = display.get_screens()
sWidth = screens[-1].width
sHeight = screens[-1].height

# Create visual stimuli
question_stim = visual.TextStim(
    win,
    text="Loading...",
    color="white",
    pos=(0, .2*yScr),
    height=fontH*1.5,
    alignText="center",
)

slider_line = visual.Rect(
    win,
    width=(.75*xScr),
    height=fontH/2,
    pos=(0, 0),
    fillColor="white",
    lineColor="white",
)

slider_marker = visual.Circle(
    win,
    radius=fontH/2,
    fillColor="red",
    lineColor="red",
    pos=(0, 0),  # Will be updated dynamically
)

number_stim = visual.TextStim(
    win,
    text="50",
    color="white",
    pos=(0, -.1*yScr),
    height=fontH,
)

left_label = visual.TextStim(
    win,
    text="Not at all",
    color="white",
    pos=(-.375*xScr, -.1*yScr),  # Position on the left side of the slider
    height=fontH,
    alignText="center",
)

right_label = visual.TextStim(
    win,
    text="Extremely",
    color="white",
    pos=(.375*xScr, -.1*yScr),  # Position on the right side of the slider
    height=fontH,
    alignText="center",
)

# Function to display a slider question
def display_slider_question(question_text, initial_value=50):
    slider_value = initial_value
    slider_marker.pos = (-.375*xScr + (slider_value/100)*.75*xScr, 0)  # Update marker position

    # Update question text
    question_stim.text = question_text
    
        # Add instruction text
    instruction_stim = visual.TextStim(
        win,
        text="Index = Up 10, Middle = Up 1 \nPinky = Down 10, Ring = Down 1\nThumb = Confirm",
        color="green",
        pos=(0, -.25*yScr),
        height=fontH,
        alignText="center",
    )
    

    while True:
        # Draw all stimuli
        question_stim.draw()
        slider_line.draw()
        slider_marker.draw()
        number_stim.text = str(slider_value)  # Update displayed number
        number_stim.draw()
        left_label.draw()
        right_label.draw()
        instruction_stim.draw()
        win.flip()

        # Wait for key input
        keys = event.waitKeys(keyList=[increment_large, decrement_large, increment_small, decrement_small, finalize])

        if increment_large in keys:
            slider_value = min(slider_value + 10, 100)
        elif decrement_large in keys:
            slider_value = max(slider_value - 10, 0)
        elif increment_small in keys:
            slider_value = min(slider_value + 1, 100)
        elif decrement_small in keys:
            slider_value = max(slider_value - 1, 0)
        elif finalize in keys:
            return slider_value  # Return the final slider value

        # Update marker position based on slider value
        slider_marker.pos = (-.375*xScr + (slider_value/100)*.75*xScr, 0)


# Define button controls
increment_large = '2'  # Move slider up by 10
decrement_large = '5'  # Move slider down by 10
increment_small = '3'  # Move slider up by 1
decrement_small = '4'  # Move slider down by 1
finalize = '1'  # Finalize the answer

# Sequence of screens
# Display the first question
positive_emotions = display_slider_question(
    "To what extent are you experiencing POSITIVE emotions RIGHT NOW?\n\n(e.g., Excited, Strong, Interested, \nEnthusiastic, Determined, Inspired", initial_value=50
)

# Display the second question
negative_emotions = display_slider_question(
    "To what extent are you experiencing NEGATIVE emotions RIGHT NOW?\n\n (e.g., Distressed, Irritable, Nervous, \nUpset, Ashamed, Scared", initial_value=50
)

# Save the responses to the file
with open(file_name, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Emotion", "Response"])
    writer.writerow(["Positive", positive_emotions])
    writer.writerow(["Negative", negative_emotions])

# Close window and quit
win.close()
core.quit()