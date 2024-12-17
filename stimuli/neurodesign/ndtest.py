import json, os, itertools, argparse, sys, re
from os.path import dirname
from neurodesign import geneticalgorithm,generate,msequence,report

# only input is subject number
sub = sys.argv[1]


# define directories
nddir = dirname(os.path.realpath(sys.argv[0]))
stimdir = dirname(nddir)
maindir = dirname(stimdir)


EXP = geneticalgorithm.experiment(
    TR=1.75,
    duration=360,
    P = [0.1667,0.1667,0.1667,0.1667,0.1667,0.1667],
    C = [[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1],
                [1,-1,0,0,0,0],[1,0,-1,0,0,0],[1,0,0,-1,0,0],[1,0,0,0,-1,0],[1,0,0,0,0,-1],
                               [0,1,-1,0,0,0],[0,1,0,-1,0,0],[0,1,0,0,-1,0],[0,1,0,0,0,-1],
                                              [0,0,1,-1,0,0],[0,0,1,0,-1,0],[0,0,1,0,0,-1],
                                                             [0,0,0,1,-1,0],[0,0,0,1,0,-1],
                                                                            [0,0,0,0,1,-1],
    ],
    n_stimuli = 6,
    rho = 0.3,
    resolution=0.1,
    stim_duration=1,
    t_pre = 0,
    t_post = 2,
    ITImodel = "exponential",
    ITImin = 2,
    ITImean = 3,
    ITImax=6
    )

# OUTPUT DIR
subdir = 'sub-%s' %(sub)
outdir = os.path.join(nddir,'testing',subdir)
if not os.path.exists(outdir):
	os.makedirs(outdir)

POP = geneticalgorithm.population(
    experiment=EXP,
    weights=[0,0.5,0.25,0.25],
    preruncycles = 20,
    cycles = 20,
    seed=1,
    outdes=5,
    folder=outdir
    )

#########################
# run natural selection #
#########################

POP.naturalselection()
POP.download()
POP.evaluate()

################
# step by step #
################

POP.add_new_designs()
POP.to_next_generation(seed=1)
POP.to_next_generation(seed=1001)

#################
# export report #
#################
report.make_report(POP,os.path.join(outdir,'report.pdf'))
