import os
import csv
import re
from psychopy import visual, event, core, gui, prefs, monitors
import pyglet
import sys
import random

# Function to generate screen variables
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

# Create window
[win_res, win] = make_screen()
xScr = float(win_res[0]) / win_res[1]
yScr = 1.
fontH = yScr / 25
wrapW = xScr / 1.5

# Initialize the window
display = pyglet.canvas.get_display()
screens = display.get_screens()
sWidth = screens[-1].width
sHeight = screens[-1].height


# Create visual stimuli
question_stim = visual.TextStim(
    win,
    text="Confirm the presented number",
    color="white",
    pos=(0, .2*yScr),
    height=fontH*1.5,
    alignText="center",
)


number_stim = visual.TextStim(
    win,
    text="50",
    color="white",
    pos=(0, -.1*yScr),
    height=fontH,
)

slider_line = visual.Rect(
    win,
    width=(0.75 * xScr),
    height=fontH / 2,
    pos=(0, 0),
    fillColor="white",
    lineColor="white",
)

slider_marker = visual.Circle(
    win,
    radius=fontH / 2,
    fillColor="red",
    lineColor="red",
    pos=(0, 0),  # Will be updated dynamically
)


        # Add instruction text
details_stim = visual.TextStim(
    win,
    text="Index = Up 10, Middle = Up 1 \nPinky = Down 10, Ring = Down 1\nThumb = Confirm",
    color="green",
    pos=(0, -.25*yScr),
    height=fontH,
    alignText="center",
)

# Define button controls
increment_large = '2'  # Move slider up by 10
decrement_large = '5'  # Move slider down by 10
increment_small = '3'  # Move slider up by 1
decrement_small = '4'  # Move slider down by 1
finalize = '1'  # Finalize the answer

# Function to display the number check question and handle slider input
def display_number_check(inputnum, initial_value=50):
    # Generate a random number within the given range
    slider_value = initial_value
    slider_marker.pos = (-0.375 * xScr + (slider_value / 100) * 0.75 * xScr, 0)  # Update marker position

        # Update question text
    question_stim.text = f"Match the slider to the number: {inputnum}"
    

    while True:
        # Draw stimuli
        question_stim.draw()
        number_stim.text = str(slider_value)  # Update displayed number
        number_stim.draw()
        slider_line.draw()
        slider_marker.draw()
        details_stim.draw()
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
            if slider_value == inputnum:
                return slider_value  # Participant successfully matched the number and pressed '1'
        slider_marker.pos = (-0.375 * xScr + (slider_value / 100) * 0.75 * xScr, 0)

            # Check if the slider value matches the random number
        



# First task: Random number between 1 and 49
rnum = random.randint(1, 49)
test1 = display_number_check(rnum)
# Second task: Random number between 51 and 99
rnum = random.randint(51, 99)
test2 = display_number_check(rnum)

# Close window and quit
win.close()
core.quit()
