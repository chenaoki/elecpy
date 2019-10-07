#!/usr/local/bin/python
from rlglue.environment.Environment import Environment
from rlglue.environment import EnvironmentLoader
from rlglue.types import Action
from rlglue.types import Observation
from rlglue.types import Reward_observation_terminal

import os
import sys
import json
import h5py
import numpy as np
import scipy
import matplotlib.pyplot as plt
import chainer
from chainer import cuda
from matplotlib import animation
from optparse import OptionParser
from scipy import ndimage

from solver.PDE import PDE
from cell.luorudy.model import model as cell_model_luorudy
from util.cmap_bipolar import bipolar
from numpy.random import *

class ElecpyEnvironment(Environment):

  def __init__(self, options):

    self.options = options

    if not os.path.isdir(options.savepath):
      os.mkdir(options.savepath)
    with open (options.param_file,'r') as f:
      sim_params = json.load(f)
    with open('{0}/sim_params.json'.format(options.savepath), 'w') as f:
      json.dump(sim_params, f, indent=4)

    self.im_h       = sim_params['geometory']['height']
    self.im_w       = sim_params['geometory']['width']
    self.N          = self.im_h*self.im_w
    self.ds         = sim_params['geometory']['ds'] # Spatial discretization step (cm)
    self.udt        = sim_params['time']['udt']     # Universal time step (ms)
    self.cnt_log    = sim_params['time']['cnt_log'] # num of udt for logging
    self.Sv         = sim_params['conductivity']['Sv']  #1400                # Surface-to-volume ratio (cm^-1)
    self.Cm         = sim_params['conductivity']['Cm']  #1.0                 # Membrane capacitance (uF/cm^2)
    self.sigma_mu   = sim_params['conductivity']['sigma_mu']  #1.74/6.25=0.2784
    self.sigma_l_i  = sim_params['conductivity']['sigma_i']['long']  #1.74   # (mS/cm)
    self.sigma_t_i  = sim_params['conductivity']['sigma_i']['tran']  #0.225  # (mS/cm)
    self.param_game   = sim_params['game']
    self.param_abl = sim_params['ablation']

    self.sigma_l    = self.sigma_l_e + self.sigma_l_i
    self.pde        = PDE( self.im_h, self.im_w, self.sigma_l_i,
                           self.sigma_t_i, self.ds, 310.15,
                           np.ones((self.im_h, self.im_w))*310.15 )

    self.cells      = cell_model_luorudy((self.N))
    self.cells.create()
    self.i_ion      = np.zeros((self.N), dtype=np.float64)
    self.i_ext_e    = np.zeros((self.N), dtype=np.float64)
    self.i_ext_i    = np.zeros((self.N), dtype=np.float64)
    self.rhs_vmem   = np.zeros((self.N), dtype=np.float64)
    self.vmem       = np.copy(cells.get_param('v'))


  def calcPenalty(self, vmem, ablationArray=None):
    penalty = 0

    if not ablationArray == None:
      penalty += np.sum(1-ablationArray)/len(ablationArray.flatten())

    if np.max(vmem) > -70:
      penalty += 1

    return penalty


  def game_setup(self):

    if self.gameIndex < self.param_game["num"]:

      print "game {0}".format(self.gameIndex),
      self.savedir = self.options.savepath + '/episode{0}'.format(self.cntEpi)
      if not os.path.isdir(self.savedir):
        os.mkdir(self.savedir)
      self.savedir += '/game{0}'.format(self.gameIndex)
      if not os.path.isdir(self.savedir):
        os.mkdir(self.savedir)

      path = self.param_game["path"]
      num_start = self.param_game["num_start"]
      num_end = self.param_game["num_end"]
      start = self.param_game["start"]
      end = self.param_game["end"]
      duration = self.param_game["duration"]
      now = 0

      session_id = randint(num_start, num_end)
      now = randint(start, end)
      self.time_end = now + duration

      with h5py.File(os.path.join(path%session_id, 'out.h5'), 'r') as f:
        group_id = '{0:0>4}'.format(now)
        self.vmem = f[group_id]['vmem'].value.flatten()
        self.cells.load(f, group_id)
        #self.cnt_udt = now * self.cnt_log
        #self.t = self.cnt_udt * self.udt
        self.gameIndex += 1

      self.penalty_org = self.calcPenalty(self.vmem, None)
      np.save('{0}/vmem_org'.format(self.savedir), self.vmem)

      return self.gameIndex

    else:
      return -1


  def createObservation(self):
    obs = Observation(numDoubles = self.im_h*self.im_w)
    tmp = self.vmem
    obs.doubleArray = list(tmp.flatten())
    return obs


  def env_init(self):
    print "[Env] init ...",

    # Display
    fig = plt.figure(figsize=(5,5))
    self.im = plt.imshow(
        np.zeros((self.im_h,self.im_w),dtype=np.float32),
        vmin = -100.0, vmax = 100.0,
        cmap=bipolar(neutral=0, lutsize=1024),
        interpolation='nearest')
    plt.axis('off')
    plt.pause(.01)

    taskSpec = str(self.im_h)
    self.cntEpi = 0

    print "done"
    return taskSpec


  def env_start(self):
    print "[Env] start ...",
    self.cntEpi     += 1
    self.t          = 0.                           # Time (ms)
    self.cnt_udt    = 0                            # Count of udt
    self.dstep      = 1                            # Time step (# of udt)
    self.cnt_save   = -1
    self.gameIndex  = 0
    assert self.game_setup() >=0
    obs = self.createObservation()
    print "done"
    return obs


  def env_step(self,action):
    print "[Env] step ...",

    time_start = self.t

    # Action control
    numAction = action.intArray[0]
    print 'action {0},'.format(numAction),
    assert numAction < 10 and numAction >= 0
    if numAction > 0:
      stim_param = {}
      stim_param['name'] = 'point'
      stim_param['start'] = self.t
      stim_param['duration'] = self.param_stim['duration']
      stim_param['interval'] = 1000 # long enough
      stim_param['amplitude'] = self.param_stim['amplitude']
      stim_param['shape'] = [self.im_h, self.im_w]
      pos_y = (numAction-1) / 3
      pos_x = (numAction-1) % 3
      stim_param['size'] = [
        self.obsOrg[0] + int(pos_y * self.obsSize / 2),
        self.obsOrg[1] + int(pos_x * self.obsSize / 2),
        self.param_stim['electrodeSize']
      ]
      self.stims.append(Stimulator(**stim_param))

    # Interval event
    if len(self.stims) < self.param_stim['maxcnt']:
      t_end = time_start + self.param_stim['interval']
    else:
      t_end = self.time_end

    with h5py.File(os.path.join(self.saveDir, 'out.h5'), 'w') as outf:
      while self.t < t_end:

        self.t = self.cnt_udt * self.udt
        self.dt = self.dstep * self.udt

        # Stimulation control
        self.i_ext_e[:,:] = 0.0
        flg_st_temp = False
        for s in self.stims:
          self.i_ext_e += s.get_current(self.t)*self.Sv
          flg_st_temp = flg_st_temp or s.get_flag(self.t)

        # step.1 cell state transition
        self.cells.set_param('dt', self.dt)
        self.cells.set_param('v', cuda.to_gpu(self.vmem))
        self.cells.update()
        self.i_ion = self.cells.get_param('it')

        # step.3 vmem
        self.rhs_vmem = self.pde_i.forward(self.vmem)
        self.rhs_vmem += self.i_ext_i
        self.rhs_vmem += self.sigma_mu * self.i_ext_e
        self.rhs_vmem *= 1/(1+self.sigma_mu)
        self.rhs_vmem -= self.i_ion * self.Sv
        self.rhs_vmem *= 1/(self.Cm*self.Sv)
        self.vmem += self.dt * self.rhs_vmem

        # Logging
        cnt_save_now = int(self.t/self.udt) // self.cnt_log
        if cnt_save_now != self.cnt_save:
          self.cnt_save = cnt_save_now
          sys.stdout.write('\r----------------------{0}/{1}ms'.format(self.t, self.time_end))
          sys.stdout.flush()

          group_id = '{0:0>4}'.format(int(self.cnt_save))
          outf.create_group(group_id)
          outf[group_id].create_dataset('vmem', data=self.vmem.reshape((self.im_h, self.im_w)))
          self.cells.save(outf, group_id)
          yield self.vmem

          self.im.set_array(self.vmem)
          plt.pause(.01)

          # Error check
          flg_error = False
          for i,v in enumerate(self.vmem.flatten()):
            if v != v :
              print "[Env] error : invalid value {1} @ {0} ms, index {2}".format(self.t, v, i)
              flg_error = True
              break
          if flg_error : break

        self.cnt_udt += self.dstep

    # Reward evaluation
    reward = 0.0

    # Game stage control
    terminal = False
    if flg_error or self.t >= self.time_end:

      im_core_dst, self.penalty_dst = self.calcPenalty(self.savedir, -1)
      np.save('{0}/coremap_dst'.format(self.savedir), im_core_dst)

      reward = self.penalty_org - self.penalty_dst
      print 'reward:{0}'.format(reward),

      if self.game_setup() < 0:
        terminal = True

    obs = self.createObservation()
    rot = Reward_observation_terminal(
        reward = reward,
        theObservation = obs,
        terminal = terminal
    )

    #return self.t, rot # for display
    print "step done"
    return rot


  def env_cleanup(self):
    pass


  def env_message(self,message):
    print  "[Env] received message:" + message
    return message


if __name__ == "__main__":
  cuda.get_device(0).use()

  parser = OptionParser()
  parser.add_option(
      '-t','--test',
      dest='test', action='store_true',default=False,
      help="test mode")
  parser.add_option(
      '-p','--param_file',
      dest='param_file', action='store', type='string', default='./rlglue_param.json',
      help="json file of simulation parameters")
  parser.add_option(
      '-d','--dst',
      dest='savepath', action='store', type='string', default='./result/data',
      help="Save data path.")
  (options, args) = parser.parse_args()
  #print 'options', options

  if not options.test:
    EnvironmentLoader.loadEnvironment(ElecpyEnvironment(options))

  else:
    objEnv = ElecpyEnvironment(options)
    objEnv.env_init()
    for epi in range(3):
      print 'Episode {0}'.format(epi)
      objEnv.env_start()
      cnt_step = 0
      while True:
        cnt_step += 1
        action = Action(numInts=1)
        action.intArray = [0]
        rot = objEnv.env_step(action)
        if rot.terminal:
          break
        else:
          print cnt_step, rot.r
