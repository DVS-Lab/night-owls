from psychopy import visual, core, event
import os
import pyglet

useFullScreen = True
useDualScreen=2
# Create a window
display = pyglet.canvas.get_display()
screens = display.get_screens()
win = visual.Window([800,600], monitor="testMonitor", units="deg", fullscr=useFullScreen, allowGUI=False, screen=len(screens)-2)


# Define instruction texts
instructions = [
    "Almost done! There is no task for this part. Please relax for about the next 10 minutes.\n\n Please keep your eyes open and remain still.\n\nPress your index to advance.",
    "We're about to start the last part!"
]

# Display first instruction and wait for '2'
instruction_text = visual.TextStim(win, text=instructions[0], color='white')
instruction_text.draw()
win.flip()
while True:
    keys = event.waitKeys()
    if '2' in keys:
        break

# Display second instruction and wait for '='
instruction_text.text = instructions[1]
instruction_text.draw()
win.flip()
while True:
    keys = event.waitKeys()
    if 'equal' in keys:
        break

# Start task: Show sphere outline for 601 seconds or until 'z' is pressed
sphere = visual.Circle(win, radius=1, edges=128, lineColor='white', fillColor=None,lineWidth=4)
start_time = core.getTime()
while core.getTime() - start_time < 601:
    sphere.draw()
    win.flip()
    keys = event.getKeys()
    if 'z' in keys:
        break

# Cleanup
win.close()
core.quit()

