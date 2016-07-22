#!/usr/local/bin/python

import os
import json
import numpy as np
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

cuda.get_device(0).use()

parser = OptionParser()
parser.add_option(
    '-p','--param_file', 
    dest='param_file', action='store', type='string', default='./sim_params.json',
    help="json file of simulation parameters")
parser.add_option(
    '-r','--restart', 
    dest='cnt_restart', action='store', type='int', default=-1, 
    help="time(ms) to restart. Set -1 to choose the last frame")
parser.add_option(
    '-d','--dst', 
    dest='savepath', action='store', type='string', default='./result/data',
    help="Save data path.")
parser.add_option(
    '-s','--src', 
    dest='srcpath', action='store', type='string', default='./result/data',
    help="Source data path for restart.")
parser.add_option(
    '-i','--isolated', 
    dest='isolated', action='store_true', default=False,
    help="Isolated mode")

(options, args) = parser.parse_args()
print options.param_file

with open (options.param_file,'r') as f : sim_params = json.load(f)
with open('{0}/sim_params.json'.format(options.savepath), 'w') as f : json.dump(sim_params, f, indent=4)
im_h         = sim_params['geometory']['height']
im_w         = sim_params['geometory']['width']
ds           = sim_params['geometory']['ds'] # Spatial discretization step (cm)
Sv           = sim_params['geometory']['Sv'] # Surface-to-volume ratio (cm^-1) 
udt          = sim_params['time']['udt']     # Universal time step (ms)
cnt_log      = sim_params['time']['cnt_log'] # num of udt for logging
time_end     = sim_params['time']['end']
stims = []
for param in sim_params['stimulation']:
  stim = Stimulator(**param)
  assert tuple(stim.shape) == (im_h, im_w)
  stims.append(stim)

Cm           = 1.0                       # Membrane capacitance (uF/cm^2)
sigma_l_i    = 3.75e-3                   # (mS/cm)
sigma_t_i    = 3.75e-4                   # (mS/cm)
sigma_l_e    = 3.75e-3                   # (mS/cm) 
sigma_t_e    = 2.14e-3                   # (mS/cm)

fig = plt.figure(figsize=(10,10))
im = plt.imshow(
    np.zeros((im_h,im_w),dtype=np.float32), 
    vmin = -100.0, vmax = 100.0, 
    cmap=bipolar(neutral=0, lutsize=1024), 
    interpolation='nearest')
plt.axis('off')

def sim( ):   
  
  print "elecpy start"
  t         = 0.                           # Time (ms)
  cnt_udt   = 0                            # Count of udt
  dstep     = 1                            # Time step (# of udt) 
  run_udt   = True                         # Flag of running sim in udt

  # Stimulators    
  print "Allocating data...",
  cell_state = createCellState((im_h,im_w))
  i_ion              = np.zeros((im_h,im_w),dtype=np.float32)
  phie               = np.zeros((im_h,im_w),dtype=np.float32)
  i_ext_e            = np.zeros((im_h,im_w),dtype=np.float32)
  i_ext_i            = np.zeros((im_h,im_w),dtype=np.float32)
  rhs_phie           = np.zeros((im_h,im_w),dtype=np.float32)
  rhs_vmem           = np.zeros((im_h,im_w),dtype=np.float32)
  vmem               = np.zeros((im_h,im_w),dtype=np.float32)
  vmem               = cell_state[lr_params.index('v')].get()
  if options.cnt_restart >= 0:
    pfx = '_{0:0>4}'.format(options.cnt_restart)
    phie = np.load('{0}/phie{1}.npy'.format(options.srcpath,pfx))
    vmem = np.load('{0}/vmem{1}.npy'.format(options.srcpath,pfx))
    cell_state = loadCellState('{0}/cell{1}'.format(options.srcpath,pfx))
    cnt_udt = options.cnt_restart * cnt_log
  print "...done"

  # PDE
  print 'Building PDE system ...',
  sigma_l      = sigma_l_e + sigma_l_i
  sigma_t      = sigma_t_e + sigma_t_i
  pde_i = PDE( im_h, im_w, sigma_l_i, sigma_t_i, ds )
  pde_m = PDE( im_h, im_w, sigma_l,   sigma_t,   ds )
  print '...done'

  if not os.path.isdir(options.savepath) : os.mkdir(options.savepath)

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

    if options.isolated:
      # Run without diffusion
      cell_state[lr_params.index('v')] -= cell_state[lr_params.index('it')] * dt 
      vmem = cell_state[lr_params.index('v')].get() 
    else:

      # step.2 phie 
      rhs_phie = i_ext_e - i_ext_i - pde_i.forward(vmem) 
      pde_cnt, phie = pde_m.solve(phie, rhs_phie)

      # step.3 vmem
      rhs_vmem = pde_i.forward(vmem)
      rhs_vmem += pde_i.forward(phie) 
      rhs_vmem -= i_ion * Sv
      rhs_vmem += i_ext_i
      rhs_vmem *= 1 / (Cm * Sv)
      vmem += dt * rhs_vmem

    # Logging & error check
    if cnt_udt%cnt_log<dstep:
      cnt_save = cnt_udt // cnt_log
      print '------------------{0}ms'.format(cnt_save)
      np.save      ('{0}/phie_{1:0>4}'.format(options.savepath,cnt_save), phie)
      np.save      ('{0}/vmem_{1:0>4}'.format(options.savepath,cnt_save), vmem)
      saveCellState('{0}/cell_{1:0>4}'.format(options.savepath,cnt_save), cell_state)
      yield vmem, t # for display

      flg = False
      for i,v in enumerate(vmem.flatten()):
        if v != v :
          print "error : invalid value {1} @ {0} ms, index {2}".format(t, v, i)
          flg = True
          break
      if flg is True:
        break

    # Time step control
    if run_udt:
      if pde_cnt < 10 and cnt_udt % 5 == 0:
        dstep = 5
        run_udt = False
    else:
      if pde_cnt > 10:
        dstep = 1
        run_udt = True

    cnt_udt += dstep

  print "elecpy done"
  exit()

def draw(data):
  img, time = data[0], data[1]
  im.set_array(img)
  return [im]

ani = animation.FuncAnimation(fig, draw, sim, blit=False, interval=200, repeat=False)
plt.show()
