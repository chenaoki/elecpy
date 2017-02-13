#!/usr/local/bin/python
from rlglue.environment.Environment import Environment
from rlglue.environment import EnvironmentLoader
from rlglue.types import Action
from rlglue.types import Observation
from rlglue.types import Reward_observation_terminal

import os
import inspect
import json
import numpy as np
import scipy
import matplotlib.pyplot as plt
import chainer
from chainer import cuda
from matplotlib import animation
from optparse import OptionParser
from scipy import ndimage

from PDE import PDE
from luorudy import luorudy
from luorudy import createCellState, loadCellState, saveCellState
from luorudy import params as lr_params
from cmap_bipolar import bipolar
from Stimulator import Stimulator
from numpy.random import *

import sys
from opmap.opmap import RawCam, VmemMap, PhaseMap, PhaseVarianceMap, CoreMap

class ElecpyEnvironment(Environment):

  def __init__(self, options):

    #super(ElecpyEnvironment, self).__init__()

    self.options = options
    #print options.param_file
    if not os.path.isdir(options.savepath) : os.mkdir(options.savepath)
    with open (options.param_file,'r') as f : sim_params = json.load(f)
    #print "sim_params", json.dumps(sim_params, indent=4)
    with open('{0}/sim_params.json'.format(options.savepath), 'w') as f : json.dump(sim_params, f, indent=4)
    
    self.im_h       = sim_params['geometory']['height']
    self.im_w       = sim_params['geometory']['width']
    self.ds         = sim_params['geometory']['ds'] # Spatial discretization step (cm)
    self.udt        = sim_params['time']['udt']     # Universal time step (ms)
    self.cnt_log    = sim_params['time']['cnt_log'] # num of udt for logging
    self.Sv         = sim_params['conductivity']['Sv']  #1400                # Surface-to-volume ratio (cm^-1) 
    self.Cm         = sim_params['conductivity']['Cm']  #1.0                 # Membrane capacitance (uF/cm^2)
    self.sigma_l_i  = sim_params['conductivity']['sigma_i']['long']  #3.75   # (mS/cm)
    self.sigma_t_i  = sim_params['conductivity']['sigma_i']['tran']  #1.5    # (mS/cm)
    self.sigma_l_e  = sim_params['conductivity']['sigma_e']['long']  #3.75   # (mS/cm) 
    self.sigma_t_e  = sim_params['conductivity']['sigma_e']['tran']  #3.00   # (mS/cm)
    self.param_game   = sim_params['game']
    self.param_stim = sim_params['stimulation']

    self.sigma_l    = self.sigma_l_e + self.sigma_l_i
    self.sigma_t    = self.sigma_t_e + self.sigma_t_i
    self.pde_i      = PDE( self.im_h, self.im_w, self.sigma_l_i, self.sigma_t_i, self.ds )
    self.pde_m      = PDE( self.im_h, self.im_w, self.sigma_l,   self.sigma_t,   self.ds )

    self.cell_state = createCellState((self.im_h,self.im_w))
    self.i_ion      = np.zeros((self.im_h,self.im_w),dtype=np.float32)
    self.phie       = np.zeros((self.im_h,self.im_w),dtype=np.float32)
    self.i_ext_e    = np.zeros((self.im_h,self.im_w),dtype=np.float32)
    self.i_ext_i    = np.zeros((self.im_h,self.im_w),dtype=np.float32)
    self.rhs_phie   = np.zeros((self.im_h,self.im_w),dtype=np.float32)
    self.rhs_vmem   = np.zeros((self.im_h,self.im_w),dtype=np.float32)
    self.vmem       = np.zeros((self.im_h,self.im_w),dtype=np.float32)
    self.vmem       = self.cell_state[lr_params.index('v')].get()

  def calcPenalty(self, path, frame_end):
    cam = RawCam(
      path='{0}'.format(path),
      cam_type='numpy',
      image_width=self.im_w,
      image_height=self.im_h,
      frame_start=0,
      frame_end=frame_end
    )
    vmem = VmemMap(cam)
    pmap = PhaseMap(vmem)
    pvmap = PhaseVarianceMap(pmap)
    coremap = CoreMap(pvmap)
    im_core, pnts_core = coremap.getFrameCore(-1)
    penalty = 0
    for p in pnts_core : penalty += ( self.im_w - p[1] )# x-axis position
    return im_core, penalty

  def game_setup(self):
    
    if self.gameIndex < self.param_game["num"]:

      print "game {0}".format(self.gameIndex),
      self.savedir    = self.options.savepath + '/episode{0}'.format(self.cntEpi) 
      if not os.path.isdir(self.savedir) : os.mkdir(self.savedir)
      self.savedir    += '/game{0}'.format(self.gameIndex) 
      if not os.path.isdir(self.savedir) : os.mkdir(self.savedir)

      path     = self.param_game["path"]
      start    = self.param_game["start"]
      end      = self.param_game["end"]
      duration = self.param_game["duration"]

      while True:
        now = randint(start, end)
        self.time_end = now + duration
        im_core_org, self.penalty_org = self.calcPenalty(
          self.param_game["path"], 
          self.time_end)
        if np.max(im_core_org) == 1: 
          np.save('{0}/coremap_org'.format(self.savedir), im_core_org)
          break      
      
      self.phie       = np.load('{0}/phie_{1:0>4}.npy'.format(path, now))
      self.vmem       = np.load('{0}/vmem_{1:0>4}.npy'.format(path, now))
      self.cell_state = loadCellState('{0}/cell_{1:0>4}'.format(path, now))
      self.cnt_udt    = now * self.cnt_log
      self.t          = self.cnt_udt * self.udt
      self.run_udt    = True
      self.flg_st     = False
      self.cnt_st_off = 0
      self.gameIndex += 1
      self.stims      = []
      elec = self.param_stim["electrodeSize"]
      self.obsOrg = np.array([
        randint(elec, self.im_h - self.obsSize - elec),
        randint(elec, self.im_w - self.obsSize - elec)])

      return self.gameIndex

    else:
      return -1

  def createObservation(self):
    obs = Observation(numDoubles = self.obsSize**2) 
    tmp = self.vmem[
      int(self.obsOrg[0]):int(self.obsOrg[0]+self.obsSize),
      int(self.obsOrg[1]):int(self.obsOrg[1]+self.obsSize)
    ]
    obs.doubleArray = list(tmp.flatten()) 
    return obs

  def env_init(self):
    print "[Env] init ...",
     
    size = int(self.param_stim["arraySize"])
    num  = int(self.param_stim["arrayNum"])
    self.obsSize = (num-1) * size

    # Display
    fig = plt.figure(figsize=(5,5))
    self.im = plt.imshow(
        np.zeros((self.im_h,self.im_w),dtype=np.float32), 
        vmin = -100.0, vmax = 100.0, 
        cmap=bipolar(neutral=0, lutsize=1024), 
        interpolation='nearest')
    plt.axis('off')
    plt.pause(.01)
    
    taskSpec = str(self.obsSize)
    self.cntEpi = 0

    print "done"
    return taskSpec 

  def env_start(self):
    print "[Env] start ...",
    self.cntEpi     += 1
    self.t          = 0.                           # Time (ms)
    self.cnt_udt    = 0                            # Count of udt
    self.dstep      = 1                            # Time step (# of udt) 
    self.run_udt    = True                         # Flag of running sim in udt
    self.flg_st     = False
    self.cnt_st_off = 0
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
      self.cell_state[lr_params.index('dt')][:,:] = self.dt 
      self.cell_state[lr_params.index('v')]=cuda.to_gpu(self.vmem)
      self.cell_state = list(luorudy().forward(tuple(self.cell_state)))
      self.i_ion = self.cell_state[lr_params.index('it')].get()

      # step.2 phie 
      self.rhs_phie = self.i_ext_e - self.i_ext_i - self.pde_i.forward(self.vmem) 
      pde_cnt, self.phie = self.pde_m.solve(self.phie, self.rhs_phie)
      self.phie -= self.phie[0,0]

      # step.3 vmem
      self.rhs_vmem = self.pde_i.forward(self.vmem)
      self.rhs_vmem += self.pde_i.forward(self.phie) 
      self.rhs_vmem -= self.i_ion * self.Sv
      self.rhs_vmem += self.i_ext_i
      self.rhs_vmem *= 1 / (self.Cm * self.Sv)
      self.vmem += self.dt * self.rhs_vmem

      # Logging
      if self.cnt_udt % self.cnt_log < self.dstep:
        cnt_save = self.cnt_udt // self.cnt_log
        print '------------------{0}ms'.format(cnt_save)
        #print '+' if self.flg_st else '-',
        #print '+' if self.run_udt else '-'
        np.save      ('{0}/phie_{1:0>4}'.format(self.savedir,cnt_save), self.phie)
        np.save      ('{0}/vmem_{1:0>4}'.format(self.savedir,cnt_save), self.vmem)
        #saveCellState('{0}/cell_{1:0>4}'.format(self.savedir,cnt_save), self.cell_state)
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

      # Time step control
      if flg_st_temp is False:
        self.cnt_st_off = 0 if self.flg_st else self.cnt_st_off + 1
      self.flg_st = flg_st_temp
      if self.run_udt:
        if self.cnt_st_off >= 3 and self.cnt_udt % 10 == 0:
          self.dstep = 2
          self.run_udt = False
      else:
        if pde_cnt > 5:
          self.dstep = 1
          self.run_udt = True

      self.cnt_udt += self.dstep

    # Reward evaluation
    reward = 0.0

    # Game stage control
    terminal = False
    if flg_error or self.t >= self.time_end:

      cnt_save = self.cnt_udt // self.cnt_log
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

