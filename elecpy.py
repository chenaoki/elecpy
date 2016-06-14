import os
from optparse import OptionParser

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

import chainer
from chainer import cuda, Function, FunctionSet, Variable
import chainer.functions as F
from chainer.functions.connection.convolution_2d import convolution_2d
xp = cuda.cupy

#from diffusion_2d import IsolatedBoundary, Parabolic
from PDE import PDE
from luorudy import luorudy
from luorudy import createCellState, loadCellState, saveCellState
from luorudy import params as lr_params

from const import const_d

from cmap_bipolar import bipolar

im_w      = 20
im_h      = 20

flag_restart   = True  # flag_restart [bool]
time_restart     = -1  # time_restart [ms]. 
savepath       = './data'
srcpath        = './data'

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
  Sv        = .3                         # Surface-to-volume ratio (cm^-1) 

  t         = 0.                           # Time (ms)
  dt        = 1.0/1000                     # Time step (ms)
  st_train  = 1                            # Number of Beats
  st_start  = 0.000                        # Time to begin stim (ms)
  st_amp    = 100.0                        # Stim amplitude (uA/cm^2)
  st_inter  = 250.                         # Basic Cycle Length (ms)
  st_dur    = 1.0                          # Stim duration (ms)
  st_time   = 0.0                          # Past time of stimulation (ms)
  st_on     = False                        # Stimulation flag 
  t_end     = st_inter*st_train            # loop end time of loop

  log_interval = 1000

  # Cell state initialization
  cell_state = createCellState((im_h,im_w))

  # Other variables
  i_ion              = xp.zeros((im_h,im_w),dtype=xp.float32)
  phie               = xp.zeros((im_h,im_w),dtype=xp.float32)
  i_ext_e            = xp.zeros((im_h,im_w),dtype=xp.float32)
  i_ext_i            = xp.zeros((im_h,im_w),dtype=xp.float32)
  rhs_phie           = xp.zeros((im_h,im_w),dtype=xp.float32)
  rhs_vmem           = xp.zeros((im_h,im_w),dtype=xp.float32)
  vmem               = xp.zeros((im_h,im_w),dtype=xp.float32)
  vmem               = cell_state[lr_params.index('v')]

  # PDE
  pde_phie_m = PDE( im_h, im_w, sigma_l,   sigma_t,   ds )
  pde_phie_i = PDE( im_h, im_w, sigma_l_i, sigma_t_i, ds )
  pde_vmem   = PDE( im_h, im_w, sigma_l_i, sigma_t_i, ds )

  os.system('mkdir {0}'.format(savepath))
  
  if flag_restart:
    pfx = '_{0:0>4}'.format(time_restart) if time_restart > 0 else ''
    phie = xp.asarray(np.load('{0}/phie{1}.npy'.format(srcpath,pfx)))
    vmem = xp.asarray(np.load('{0}/vmem{1}.npy'.format(srcpath,pfx)))
    cell_state = loadCellState('{0}/cell{1}'.format(srcpath,pfx))
    t = float(time_restart) + dt # to avoid first pacing stimulation

  while t < t_end:

    # Stimulation control
    if int( (t - st_start)/dt ) % int(st_inter/dt) == 0:
      st_on = True
      print 'stim on'
      i_ext_e[1   ,1:-1] = -st_amp
      i_ext_e[-2  ,1:-1] = -st_amp
      i_ext_e[1:-1,1   ] = -st_amp
      i_ext_e[1:-1,-2  ] = -st_amp
      i_ext_e[im_h/2-2:im_h/2+2,im_w/2-2:im_w/2+2] = +(st_amp * (im_h+im_w-4) / 8 )
      #cell_state[lr_params.index('st')][im_h/2-2:im_h/2+2,im_w/2-2:im_w/2+2] = st_amp
    if st_on:
      st_time += dt 
      if st_time > st_dur:
        print 'stim off'
        st_on = False
        i_ext_e[:,:] = 0.0 
        #cell_state[lr_params.index('st')][:,:] = 0.0
        st_time = 0.0

    # step.1 cell state transition
    cell_state[lr_params.index('v')]=vmem 
    cell_state = list(luorudy().forward(tuple(cell_state)))

    if False:
      cell_state[lr_params.index('v')] -= cell_state[lr_params.index('it')] * dt 
      vmem = cell_state[lr_params.index('v')] 
    else:
      i_ion = cell_state[lr_params.index('it')]

      # step.2 phie 
      rhs_phie = i_ext_e - i_ext_i - pde_vmem.forward(vmem) 
      cnt, phie = pde_phie_m.solve(phie, rhs_phie, tol=1e-4)
      #print 'pde iteratino {0}'.format(cnt)

      # step.3 vmem
      rhs_vmem = pde_vmem.forward(vmem)
      rhs_vmem += pde_phie_i.forward(phie) 
      rhs_vmem -= i_ion  * Sv
      rhs_vmem += i_ext_i
      rhs_vmem *= 1 / (Cm * Sv)  
      vmem += dt * rhs_vmem

    #flg = False
    #for i,v in enumerate(vmem.get().flatten()):
    #  if v != v :
    #    print "error : invalid value {1} @ {0} ms, index {2}".format(t, v, i)
    #    flg = True
    #    break
    #if flg is True:
    #  break
    
    step = int(t / dt)
    if step%log_interval==0:
      t_int = int(step/log_interval) * int( dt * log_interval )
      print '------------------{0}ms'.format(t_int)
      np.save(      '{0}/phie_{1:0>4}'.format(savepath,t_int), phie.get())
      np.save(      '{0}/vmem_{1:0>4}'.format(savepath,t_int), vmem.get())
      saveCellState('{0}/cell_{1:0>4}'.format(savepath,t_int), cell_state)
      yield vmem.get(), t # for display

    t += dt

  np.save('{0}/phie'.format(savepath), phie.get()) 
  np.save('{0}/vmem'.format(savepath), vmem.get()) 
  saveCellState('{0}/cell'.format(savepath), cell_state)
  print "elecpy done"

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
      '-r','--restart', 
      dest='time_restart', action='store', type='int', 
      help="time(ms) to restart. Set -1 to choose the last frame")
  parser.add_option(
      '-d','--dst', 
      dest='savepath', action='store', type='string', 
      help="Save data path. (default:./data)")
  parser.add_option(
      '-s','--src', 
      dest='srcpath', action='store', type='string', 
      help="Source data path for restart. (default:./data)")

  (options, args) = parser.parse_args()
  flag_restart = True                 if options.time_restart else False
  time_restart = options.time_restart if options.time_restart else 0
  savepath     = options.savepath     if options.savepath     else './data'
  srcpath      = options.srcpath      if options.srcpath      else './data'

  ani = animation.FuncAnimation(fig, draw, sim, blit=False, interval=1, repeat=False)
  plt.show()
