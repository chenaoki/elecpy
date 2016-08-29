#!/usr/local/bin/python
from rlglue.environment.Environment import Environment
from rlglue.environment import EnvironmentLoader
from rlglue.types import Action
from rlglue.types import Observation
from rlglue.types import Reward_observation_terminal

import os
import json
import numpy as np
import scipy
import matplotlib.pyplot as plt
import chainer
from chainer import cuda
from matplotlib import animation
from optparse import OptionParser

from PDE import PDE
from luorudy import luorudy
from luorudy import createCellState, loadCellState, saveCellState
from luorudy import params as lr_params
from cmap_bipolar import bipolar
from Stimulator import Stimulator

class ElecpyEnvironment(Environment):
  obsImageSie = 50

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
    self.gameList   = sim_params['games']

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
    self.stims      = []
    

  def game_setup(self):
    
    if self.gameIndex < len(self.gameList):

      game = self.gameList[self.gameIndex]
      # print "[Env] game", self.gameIndex, json.dumps( game, indent=4),

      self.phie       = np.load('{path}/phie_{start:0>4}.npy'.format(**game))
      self.vmem       = np.load('{path}/vmem_{start:0>4}.npy'.format(**game))
      self.cell_state = loadCellState('{path}/cell_{start:0>4}'.format(**game))
      self.time_end   = game["end"]
      self.time_int   = game["interval"]
      self.cnt_udt    = game["start"] * self.cnt_log
      self.t          = self.cnt_udt * self.udt
      self.run_udt    = True
      self.flg_st     = False
      self.cnt_st_off = 0
      self.gameIndex += 1
      return self.gameIndex

    else:
      return -1

  def createObservation(self):
    obs = Observation(
        numDoubles = self.obsImageSie * self.obsImageSie) 
    tmp = scipy.misc.imresize(self.vmem, (self.obsImageSie, self.obsImageSie)) 
    obs.doubleArray = list(tmp.flatten()) 
    return obs

  def env_init(self):
    print "[Env] init ...",
    assert self.obsImageSie > 0
    taskSpec = str(self.obsImageSie)
    print "done"
    return taskSpec 

  def env_start(self):
    print "[Env] start ...",
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

    # Interval event
    while self.t < time_start + self.time_int:

      self.t = self.cnt_udt * self.udt
      self.dt = self.dstep * self.udt
      
      # Stimulation control
      self.i_ext_e[:,:] = 0.0
      flg_st_temp = False
      for s in self.stims:
        self.i_ext_e += s.get_current(t)*self.Sv
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

      # Logging & error check
      flg_error = False
      if self.cnt_udt % self.cnt_log < self.dstep:
        cnt_save = self.cnt_udt // self.cnt_log
        #print '------------------{0}ms'.format(cnt_save)
        #print '+' if self.flg_st else '-',
        #print '+' if self.run_udt else '-'
        np.save      ('{0}/phie_{1:0>4}'.format(self.options.savepath,cnt_save), self.phie)
        np.save      ('{0}/vmem_{1:0>4}'.format(self.options.savepath,cnt_save), self.vmem)
        saveCellState('{0}/cell_{1:0>4}'.format(self.options.savepath,cnt_save), self.cell_state)

      # Error check
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
      if self.game_setup() < 0:
        terminal = True

    obs = self.createObservation() 
    rot = Reward_observation_terminal(
        reward = reward,
        theObservation = obs,
        terminal = terminal
    )

    #return self.t, rot # for display
    print "time {0},".format(self.t),
    print "done"
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
    
    fig = plt.figure(figsize=(5,5))
    im = plt.imshow(
        np.zeros((50,50),dtype=np.float32), 
        vmin = -100.0, vmax = 100.0, 
        cmap=bipolar(neutral=0, lutsize=1024), 
        interpolation='nearest')
    plt.axis('off')

    def sim():
      objEnv = ElecpyEnvironment(options)
      objEnv.env_init()
      for epi in range(3):
        print 'Episode {0}'.format(epi)
        objEnv.env_start()
        while True:
          rot = objEnv.env_step(0)
          if rot.terminal:
            break
          else:
            yield [rot]

      exit()

    def draw(data):
      rot = data[0]
      tmp = np.array(rot.o.doubleArray)
      tmp = tmp.reshape((50,50))
      im.set_array(tmp)
      return [im]

    ani = animation.FuncAnimation(fig, draw, sim, blit=False, interval=200, repeat=False)
    plt.show()
