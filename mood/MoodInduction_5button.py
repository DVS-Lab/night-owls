import os
import csv
import textwrap
from psychopy import visual, event, sound, core, gui, prefs

#get subjID
subjDlg=gui.Dlg(title="Mood Induction")
subjDlg.addField('Subject:')
subjDlg.addField('Session:', choices=['1', '2','3','4','5','6','7','8','9','10','11','12'])
subjDlg.show()

if subjDlg.show():  # This displays the dialog
    subj_id = subjDlg.data[0]  # Extract the Subject ID
    ses = subjDlg.data[1]  # Extract the selected session
else:
    core.quit()  # Gracefully exit if "Cancel" is pressed
    
    
#### Pre Induction Mood Check ###
    
# Folder structure
folder_path = f"logs/sub-{subj_id}/ses-{ses}"

# Check if the folder exists, and create it if not
if not os.path.exists(folder_path):
    os.makedirs(folder_path)
    
# Prepare the file name based on the inputs
file_name = f"{folder_path}/sub-{subj_id}_ses-{ses}_obs-3_mood.csv"

# Initialize the window
win = visual.Window(
    size=[2400, 1800], fullscr=False, color="black", units="pix"
)

# Create visual stimuli
question_stim = visual.TextStim(
    win,
    text="Loading...",
    color="white",
    pos=(0, 200),
    height=100,
    alignText="center",
)

slider_line = visual.Rect(
    win,
    width=600,
    height=10,
    pos=(0, 0),
    fillColor="white",
    lineColor="white",
)

slider_marker = visual.Circle(
    win,
    radius=15,
    fillColor="red",
    lineColor="red",
    pos=(0, 0),  # Will be updated dynamically
)

number_stim = visual.TextStim(
    win,
    text="50",
    color="white",
    pos=(0, -50),
    height=100,
)

left_label = visual.TextStim(
    win,
    text="Not at all",
    color="white",
    pos=(-300, -75),  # Position on the left side of the slider
    height=100,
    alignText="center",
)

right_label = visual.TextStim(
    win,
    text="Extremely",
    color="white",
    pos=(300, -75),  # Position on the right side of the slider
    height=100,
    alignText="center",
)

# Function to display a slider question
def display_slider_question(question_text, initial_value=50):
    slider_value = initial_value
    slider_marker.pos = (-300 + (slider_value * 6), 0)  # Update marker position

    # Update question text
    question_stim.text = question_text

    while True:
        # Draw all stimuli
        question_stim.draw()
        slider_line.draw()
        slider_marker.draw()
        number_stim.text = str(slider_value)  # Update displayed number
        number_stim.draw()
        left_label.draw()
        right_label.draw()
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
        slider_marker.pos = (-300 + (slider_value * 6), 0)

# Define button controls
increment_large = 'right'  # Move slider up by 10
decrement_large = 'left'  # Move slider down by 10
increment_small = 'up'  # Move slider up by 1
decrement_small = 'down'  # Move slider down by 1
finalize = 'space'  # Finalize the answer

# Sequence of screens
# Display the first question
positive_emotions = display_slider_question(
    "To what extent are you experiencing POSITIVE emotions RIGHT NOW?", initial_value=50
)

# Display the second question
negative_emotions = display_slider_question(
    "To what extent are you experiencing NEGATIVE emotions RIGHT NOW?", initial_value=50
)

# Save the responses to the file
with open(file_name, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Emotion", "Response"])
    writer.writerow(["Positive", positive_emotions])
    writer.writerow(["Negative", negative_emotions])




####### Mood Induction #######

# Instruction page
instruction_text = visual.TextStim(
    win,
    text="You are about to see the memory you described earlier. Please reminisce on different parts of this memory for the next 10 minutes while listening to music. Feel free to daydream about the memory (relive it), but please keep your eyes open. Your goal is to get in as good of a mood as possible with this memory!", 
    color="white",
    pos=(0, 0),
    height=30,
    alignText="center",
)

# Show the instruction page for 20 seconds
instruction_text.draw()
win.flip()
core.wait(2) 

# Load and wrap text
def load_and_wrap_text(filename, max_width=60):
    with open(filename, 'r') as file:
        content = file.read()
    wrapped_text = textwrap.fill(content, width=max_width)
    return wrapped_text


# Load text from file and wrap
filename = "example.txt"  # Replace with text file path
max_width = 60  # Adjust as needed for text wrapping
text_content = load_and_wrap_text(filename, max_width)

# Create visual stimuli
text_stim = visual.TextStim(
    win,
    text=text_content,
    wrapWidth=1.8,
    color="white",
    alignText="center",
    height=20,
)

# Load audio files
audio_files = ["song1.wav", "song2.wav", "song3.wav", "song4.wav"]  
#audio_files = ["test1.wav", "test2.wav"]  

prefs.general['audioLib'] = ['sounddevice']
import sounddevice as sd
print(sd.query_devices())
sd.default.device[1] = 3  # Replace with preferred device ID

audio_sounds = [sound.Sound(file) for file in audio_files]

# Total time to display the text and play audio
display_duration = 15  # in seconds

# Display text while playing audio sequentially
text_stim.draw()  # Draw the text on screen
win.flip()  # Update window to show the text

# Play each audio file one after the other
for audio in audio_sounds:
    audio.play()  # Start playing the audio
    core.wait(audio.getDuration()) # Wait for the audio to finish
    
    # After all audio files are played, wait for the remaining time to reach 600 seconds
core.wait(display_duration - sum([audio.getDuration() for audio in audio_sounds]))
    

#### Post Mood Induction Check #### 

# Prepare the file name based on the inputs
file_name2 = f"{folder_path}/sub-{subj_id}_ses-{ses}_obs-4_mood.csv"
positive_emotions = display_slider_question("To what extent are you experiencing POSITIVE emotions RIGHT NOW?", initial_value=50)
negative_emotions = display_slider_question("To what extent are you experiencing NEGATIVE emotions RIGHT NOW?", initial_value=50)


# Save the responses to the file
with open(file_name2, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Emotion", "Response"])
    writer.writerow(["Positive", positive_emotions])
    writer.writerow(["Negative", negative_emotions])

# Close window and quit
win.close()
core.quit()
