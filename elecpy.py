#!/usr/local/bin/python

import os
from optparse import OptionParser
import json

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

import chainer
from chainer import cuda, Function, FunctionSet, Variable
import chainer.functions as F
from chainer.functions.connection.convolution_2d import convolution_2d

#from diffusion_2d import IsolatedBoundary, Parabolic
from PDE import PDE
from luorudy import luorudy
from luorudy import createCellState, loadCellState, saveCellState
from luorudy import params as lr_params

from const import const_d

from cmap_bipolar import bipolar

from Stimulator import Stimulator

im_h      = 120#180#30#
im_w      = 180#30#90#

flag_restart   = True  # flag_restart [bool]
time_restart     = -1  # time_restart [ms]. 
time_end       = 1000
savepath       = './data'
srcpath        = './data'
stim_param_file = './stim_params.json'
debug = False
isolated = False

def sim( ):   
  
  print "elecpy start"

  ds        = 1.5e-2                    # Spatial discretization step (cm)
  sigma_l_e = 3.75e-3                   # (mS/cm) 
  sigma_l_i = 3.75e-3                   # (mS/cm)
  sigma_t_e = 2.14e-3                   # (mS/cm)
  sigma_t_i = 3.75e-4                   # (mS/cm)
  sigma_l   = sigma_l_e + sigma_l_i     # (mS/cm)
  sigma_t   = sigma_t_e + sigma_t_i     # (mS/cm)
  Cm        = 1.0                       # Membrane capacitance (uF/cm^2)
  Sv        = .3                        # Surface-to-volume ratio (cm^-1) 

  t         = 0.                           # Time (ms)
  udt       = 1.0/1000                     # Universal time step (ms)
  cnt_udt   = 0                            # Count of udt
  dstep     = 1                            # Time step (# of udt) 
  run_udt   = True                         # Flag of running sim in udt

  # Stimulators
  stims = []
  with open (stim_param_file,'r') as f : stim_params = json.load(f)
  print stim_params
  for param in stim_params:
    print param
    stim = Stimulator(**param)
    print stim.shape, im_h, im_w
    assert tuple(stim.shape) == (im_h, im_w)
    stims.append(stim)
    
  log_interval = 1 if debug else 1000

  print "Allocating data...",
  # Cell state initialization
  cell_state = createCellState((im_h,im_w))
  # Other variables
  i_ion              = np.zeros((im_h,im_w),dtype=np.float32)
  phie               = np.zeros((im_h,im_w),dtype=np.float32)
  i_ext_e            = np.zeros((im_h,im_w),dtype=np.float32)
  i_ext_i            = np.zeros((im_h,im_w),dtype=np.float32)
  rhs_phie           = np.zeros((im_h,im_w),dtype=np.float32)
  rhs_vmem           = np.zeros((im_h,im_w),dtype=np.float32)
  vmem               = np.zeros((im_h,im_w),dtype=np.float32)
  vmem               = cell_state[lr_params.index('v')].get()
  print "done"

  # PDE
  print 'Building PDE system ...',
  pde_phie_m = PDE( im_h, im_w, sigma_l,   sigma_t,   ds )
  pde_phie_i = PDE( im_h, im_w, sigma_l_i, sigma_t_i, ds )
  pde_vmem   = PDE( im_h, im_w, sigma_l_i, sigma_t_i, ds )
  print 'done'

  if not os.path.isdir(savepath) : os.mkdir(savepath)

  if flag_restart:
    pfx = '_{0:0>4}'.format(time_restart) if time_restart > 0 else ''
    phie = np.load('{0}/phie{1}.npy'.format(srcpath,pfx))
    vmem = np.load('{0}/vmem{1}.npy'.format(srcpath,pfx))
    cell_state = loadCellState('{0}/cell{1}'.format(srcpath,pfx))
    cnt_udt = time_restart if debug else 1000 * time_restart 

  while t < time_end:

    t = cnt_udt * udt
    dt = dstep * udt

    # Stimulation control
    i_ext_e[:,:] = 0.0
    for s in stims:
      i_ext_e += s.get_current(t)

    # step.1 cell state transition
    cell_state[lr_params.index('dt')][:,:] = dt 
    cell_state[lr_params.index('v')]=cuda.to_gpu(vmem)
    cell_state = list(luorudy().forward(tuple(cell_state)))
    i_ion = cell_state[lr_params.index('it')].get()

    if isolated:
      print 'Running without diffusion...',
      cell_state[lr_params.index('v')] -= cell_state[lr_params.index('it')] * dt 
      vmem = cell_state[lr_params.index('v')].get() 
      print 'done'

    else:
      # step.2 phie 
      rhs_phie = i_ext_e - i_ext_i - pde_vmem.forward(vmem) 
      pde_cnt, phie = pde_phie_m.solve(phie, rhs_phie)
      #print 'pde iteration {0} @ {1} ms'.format(cnt, t)

      # step.3 vmem
      rhs_vmem = pde_vmem.forward(vmem)
      rhs_vmem += pde_phie_i.forward(phie) 
      rhs_vmem -= i_ion  * Sv
      rhs_vmem += i_ext_i
      rhs_vmem *= 1 / (Cm * Sv)  
      vmem += dt * rhs_vmem

    if cnt_udt%log_interval<dstep:
      t_int = cnt_udt if debug else cnt_udt // log_interval
      print '------------------{0}{1}'.format(t_int, 'us' if debug else 'ms')
      np.save(      '{0}/phie_{1:0>4}'.format(savepath,t_int), phie)
      np.save(      '{0}/vmem_{1:0>4}'.format(savepath,t_int), vmem)
      saveCellState('{0}/cell_{1:0>4}'.format(savepath,t_int), cell_state)
      yield vmem, t # for display

      flg = False
      for i,v in enumerate(vmem.flatten()):
        if v != v :
          print "error : invalid value {1} @ {0} ms, index {2}".format(t, v, i)
          flg = True
          break
      if flg is True:
        break

    if debug:
      #i_debug = int(im_h * im_w / 2 )
      i_debug = 188
      print '@', i_debug,
      print 't:', t,
      print 'i_ion:',i_ion.flatten()[i_debug],
      print 'vmem:',vmem.flatten()[i_debug],
      print 'phie:',phie.flatten()[i_debug],
      print 'rhs_vmem:',rhs_vmem.flatten()[i_debug]

    if run_udt:
      if pde_cnt < 10 and cnt_udt % 5 == 0:
        dstep = 5
        run_udt = False
    else:
      if pde_cnt > 10:
        dstep = 1
        run_udt = True

    cnt_udt += dstep

    #print '{0:.3f}'.format(t), pde_cnt, dstep

  np.save('{0}/phie'.format(savepath), phie) 
  np.save('{0}/vmem'.format(savepath), vmem) 
  saveCellState('{0}/cell'.format(savepath), cell_state)
  with open('{0}/stims.json'.format(savepath), 'w') as f : json.dump(stims, f, indent=4)
  print "elecpy done"
  exit()

fig = plt.figure(figsize=(10,10))
im = plt.imshow(
    np.zeros((im_h,im_w),dtype=np.float32), 
    vmin = -100.0, vmax = 100.0, 
    cmap=bipolar(neutral=0, lutsize=1024), 
    interpolation='nearest')
plt.axis('off')

def draw(data):
  img, time = data[0], data[1]
  im.set_array(img)
  return [im]

if __name__ == '__main__':

  cuda.get_device(0).use()

  parser = OptionParser()
  parser.add_option(
      '-e','--end', 
      dest='time_end', action='store', type='int', 
      help="time(ms) to end.")
  parser.add_option(
      '-r','--restart', 
      dest='time_restart', action='store', type='int', 
      help="time(ms) to restart. Set -1 to choose the last frame")
  parser.add_option(
      '-d','--dst', 
      dest='savepath', action='store', type='string', 
      help="Save data path.")
  parser.add_option(
      '-s','--src', 
      dest='srcpath', action='store', type='string', 
      help="Source data path for restart.")
  parser.add_option(
      '-v','--verbose', 
      dest='debug', action='store_true', 
      help="Debug mode")
  parser.add_option(
      '-i','--isolated', 
      dest='isolated', action='store_true', 
      help="Isolated mode")
  parser.add_option(
      '-p','--stim_param_file', 
      dest='stim_param_file', action='store', type='string',
      help="json file containing stimlation parameters")

  (options, args) = parser.parse_args()
  isolated        = True                    if options.isolated        else False
  debug           = True                    if options.debug           else False
  flag_restart    = True                    if options.time_restart    else False
  time_restart    = options.time_restart    if options.time_restart    else 0
  time_end        = options.time_end        if options.time_end        else 1000
  savepath        = options.savepath        if options.savepath        else './result/data'
  srcpath         = options.srcpath         if options.srcpath         else './result/data'
  stim_param_file = options.stim_param_file if options.stim_param_file else './stim_params.json'

  ani = animation.FuncAnimation(fig, draw, sim, blit=False, interval=200, repeat=False)
  plt.show()
