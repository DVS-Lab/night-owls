from psychopy import visual, event, core

# Initialize the window
win = visual.Window(
    size=[800, 600], fullscr=False, color="black", units="pix"
)

# Create visual stimuli
question_stim = visual.TextStim(
    win,
    text="",
    color="white",
    pos=(0, 200),
    height=0.05,
    alignText="center",
)

slider_line = visual.Rect(
    win,
    width=600,
    height=2,
    pos=(0, 0),
    fillColor="white",
    lineColor="white",
)

slider_marker = visual.Circle(
    win,
    radius=10,
    fillColor="red",
    lineColor="red",
    pos=(0, 0),  # Will be updated dynamically
)

number_stim = visual.TextStim(
    win,
    text="50",
    color="white",
    pos=(0, -50),
    height=0.05,
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

# Save the responses to a file
with open("slider_responses.csv", "w") as file:
    file.write("Question,Response\n")
    file.write(f"Positive emotions,{positive_emotions}\n")
    file.write(f"Negative emotions,{negative_emotions}\n")

# Close window and quit
win.close()
core.quit()
