import numpy as np
import time

import chainer
import chainer.functions as F
from chainer import cuda, Function, FunctionSet, Variable
from chainer.functions.connection.convolution_2d import convolution_2d
xp = cuda.cupy

from const import const_d
import os

params = [ 'v',  'dt',  'h',  'j',  'm',  'nai',  'naiss',  'cai',  'caiss',  'ki',  'kiss',  'ltypeCzero',  'ltypeCone',  'ltypeCtwo',  'ltypeCthree',  'ltypeIVf',  'ltypeIVs',  'ltypeO',  'fcasc',  'fmode0',  'ical',  'ilcana',  'ilcak',  'irel',  'itr',  'b',  'g',  'xr',  'xs1',  'xs2',  'inacass',  'nsr',  'ryrCone',  'ryrCtwo',  'ryrCthree',  'ryrCfour',  'ryrOone',  'ryrIone',  'ryrItwo',  'ryrIthree',  'ryrIfour',  'ryrIfive',  'csqn',  'jsr', 'it', 'st']

string_args = ''
string_rets = ''
for i, key in enumerate(params) : 
  string_args += 'T {0}, '.format(key)
  string_rets += 'T _{0}, '.format(key)
string_args = string_args.rstrip(', ')
string_rets = string_rets.rstrip(', ')
elemwise = open('ElementwiseKernel.c').read().format(**const_d)

def createCellState(shape):
  state = [ xp.asarray( np.ones(shape, dtype=np.float32) * const_d[param+'_'] ) for param in params ]
  return state

def loadCellState(path):
  state = [ xp.asarray( np.load(path+'/'+param+'.npy')) for param in params ]
  return state

def saveCellState(path, state):
  os.system('mkdir '+path)
  for i, param in enumerate(params):
    np.save(path+'/'+param, state[i].get()) 
  pass

class luorudy(Function):
  def forward(self,x):
    state = x 
    _state = cuda.elementwise(
      string_args, 
      string_rets,
      elemwise,
      'luorudy')(*state) 
    return _state

if __name__ == '__main__':

  print 'LuoRudy model on chainer'
  print 'Unit test : pacing'

  #print const_d
  #print params 
  #print string_args 
  #print string_rets 
  #print elemwise

  t = 0.                                   # Time (ms)
  dt = const_d['dt_']                      # Time step (ms)
  st_train  = 1                            # Number of Beats
  st_start = 10.                           # Time to begin stim (ms)
  st_amp = -100                            # Stim amplitude (uA/cm^2)
  st_inter = 250.                          # Basic Cycle Length (ms)
  st_dur = .5                              # Stim duration (ms)
  st_time = 0.0                            # Past time of stimulation (ms)
  st_on = False                            # Stimulation flag 
  steps = int((st_inter*st_train)/dt)      # Number of loop
  
  # Cell state initialization
  cell_state = createCellState((50,50))
  #print '-------------'
  #for j, (key, param) in enumerate(zip(params, cell_state)) : 
  #  print j, key, param[0,0]

  log_interval = 1000
  log_v = np.zeros((steps/log_interval))
  
  start = time.time()
  
  for step in range(steps):
  #for step in range(1):

    # Stimulation setting
    if int( (t - st_start)/dt ) % int(st_inter/dt) == 0:
      st_on = True
    if st_on:
      cell_state[params.index('st')] = st_amp 
      st_time += dt 
      if st_time > st_dur:
        st_on = False
        st_time = 0.0
    else: # st_on
      cell_state[params.index('st')] = 0.0 

    # State transition
    cell_state = list(luorudy().forward(tuple(cell_state)))

    # Membrane voltage update
    cell_state[params.index('v')] -= cell_state[params.index('it')] * dt 

    if step % log_interval == 0:
      print '-------------{0}(ms)'.format(int(t))
      print 'v:',cell_state[params.index('v')][0,0] 
      print 'it:',cell_state[params.index('it')][0,0]
      print 'st:',cell_state[params.index('st')][0,0]
      log_v[step/log_interval] = cell_state[params.index('v')][0,0] 

    t += dt 

  elapsed_time = time.time() - start
  print 'elapsed_time:', elapsed_time
  np.save('log_v', log_v)
