from psychopy import visual, core, event
win = visual.Window([800,600], monitor="testMonitor", units="deg")

button_box = visual.ImageStim(win, image = './buttons.png',size=(20,20), pos=[0,0])
index_r = visual.ImageStim(win, image = './index_right.png',size=(20,20), pos=[0,0])
middle_r = visual.ImageStim(win, image = './middle_right.png',size=(20,20), pos=[0,0])
ring_r = visual.ImageStim(win, image = './ring_right.png',size=(20,20), pos=[0,0])
thumb_r = visual.ImageStim(win, image = './thumb_right.png',size=(20,20), pos=[0,0])
pinky_r = visual.ImageStim(win, image = './pinky_right.png',size=(20,20), pos=[0,0])
index_l = visual.ImageStim(win, image = './index_left.png',size=(20,20), pos=[0,0])
middle_l = visual.ImageStim(win, image = './middle_left.png',size=(20,20), pos=[0,0])
ring_l = visual.ImageStim(win, image = './ring_left.png',size=(20,20), pos=[0,0])
thumb_l = visual.ImageStim(win, image = './thumb_left.png',size=(20,20), pos=[0,0])
pinky_l = visual.ImageStim(win, image = './pinky_left.png',size=(20,20), pos=[0,0])

responseKeys=('1','2','3','4','5','6','7','8','9','a','z')
resp = event.getKeys(keyList = responseKeys)


while core.Clock().getTime() < (60):
    button_box.draw()
    win.flip()
    resp = event.getKeys(keyList = responseKeys)
    if len(resp)>0:
        if resp[0] == 'z':
            win.close()
            core.quit()
        elif resp[0] == '1':
            thumb_r.draw()
            win.flip()
            core.wait(.5)
            button_box.draw()
            win.flip()
        elif resp[0] == '2':
            index_r.draw()
            win.flip()
            core.wait(.5)
            button_box.draw()
            win.flip()
        elif resp[0] == '3':
            middle_r.draw()
            win.flip()
            core.wait(.5)
            button_box.draw()
            win.flip()
        elif resp[0] == '4':
            ring_r.draw()
            win.flip()
            core.wait(.5)
            button_box.draw()
            win.flip()
        elif resp[0] == '5':
            pinky_r.draw()
            win.flip()
            core.wait(.5)
            button_box.draw()
            win.flip()
        elif resp[0] == '6':
            thumb_l.draw()
            win.flip()
            core.wait(.5)
            button_box.draw()
            win.flip()
        elif resp[0] == '7':
            index_l.draw()
            win.flip()
            core.wait(.5)
            button_box.draw()
            win.flip()
        elif resp[0] == '8':
            middle_l.draw()
            win.flip()
            core.wait(.5)
            button_box.draw()
            win.flip()
        elif resp[0] == '9':
            ring_l.draw()
            win.flip()
            core.wait(.5)
            button_box.draw()
            win.flip()
        elif resp[0] == 'a':
            pinky_l.draw()
            win.flip()
            core.wait(.5)
            button_box.draw()
            win.flip()


