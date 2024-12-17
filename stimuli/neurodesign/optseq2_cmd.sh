#!/bin/bash

optseq2 --ntp 212 --tr 1.7 \
        --psdwin 0 25.5 0.85 \
        --ev evt1 5.1 11 \
        --ev evt2 5.1 5 \
        --ev evt3 5.1 11 \
        --ev evt4 5.1 11 \
        --ev evt5 5.1 5 \
        --ev evt6 5.1 11 \
        --evc -1 0 1 -1 0 1 \
        --evc 1 0 -1 -1 0 1 \
        --evc -1 -1 -1 1 1 1 \
        --evc 1 -2 1 1 -2 1 \
        --tnullmin 0.85 --tnullmax 8.5 \
        --nkeep 3 \
        --o rf1-test \
        --nsearch 1000 



        # evt1	computer_loss
        # evt2	computer_neutral
        # evt3	computer_win
        # evt4	stranger_loss
        # evt5	stranger_neutral
        # evt6	stranger_win
