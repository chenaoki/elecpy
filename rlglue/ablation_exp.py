# -*- coding: utf-8 -*-
"""
Simple RL glue experiment setup
"""

import numpy as np
import rlglue.RLGlue as RLGlue

epi = 1
epi_learn = 0
epi_max = 1000

print "[Exp] starting up!"
RLGlue.RL_init()

while epi_learn < epi_max:

  learning = (np.mod(epi, 10) != 0)

  # Evaluate model every 10 episodes
  if learning:
    RLGlue.RL_agent_message("unfreeze learning")
  else:
    RLGlue.RL_agent_message("freeze learning")
    epi_learn += 1

  print "[Exp] Episode {0}, Learning {1},".format(epi, learning),
  RLGlue.RL_episode(0)
  totalSteps = RLGlue.RL_num_steps()
  totalReward = RLGlue.RL_return()
  epi += 1
  print "{0} steps, total reward {1}".format(totalSteps, totalReward)

  if np.mod(epi_learn, 100) == 0 and epi_learn != 0:
    print "SAVE CURRENT MODEL"
    RLGlue.RL_agent_message("save model")

RLGlue.RL_cleanup()

print "Experiment COMPLETED @ Episode ", epi
