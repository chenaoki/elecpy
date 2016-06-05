import numpy as np

import chainer
from chainer import cuda, Function, FunctionSet, Variable
import chainer.functions as F
from chainer.functions.connection.convolution_2d import convolution_2d
xp = cuda.cupy

#from diffusion_2d import IsolatedBoundary, Parabolic
from PDE import PDE
from luorudy import luorudy, createCellState
from luorudy import params as lr_params

import matplotlib.pyplot as plt
from matplotlib import animation

from const import const_d

# sim params
ds        = 1.5e-2                    # Spatial discretization step (cm)
sigma_l_e = 3.75e-3                   # (mS/cm) 
sigma_l_i = 3.75e-3                   # (mS/cm)
sigma_t_e = 2.14e-3                   # (mS/cm)
sigma_t_i = 3.75e-4                   # (mS/cm)
sigma_l   = sigma_l_e + sigma_l_i     # (mS/cm)
sigma_t   = sigma_t_e + sigma_t_i     # (mS/cm)
Cm        = 1.0                       # Membrane capacitance (uF/cm^2)
Sv        = 3000                      # Surface-to-volume ratio (cm^-1) 

im_w      = 50
im_h      = 50

cuda.get_device(0).use()
fig = plt.figure(figsize=(10,10))
im = plt.imshow(
    np.zeros((im_h,im_w),dtype=np.float32), 
    vmin = -100.0, vmax = 100.0, cmap='hot', interpolation='nearest')
plt.axis('off')

def sim():

  print "elecpy start"

  # Cell state initialization
  cell_state = createCellState((im_h,im_w))

  # Other variables
  i_ion              = xp.zeros((im_h,im_w),dtype=xp.float32)
  phie               = xp.zeros((im_h,im_w),dtype=xp.float32)
  i_ext              = xp.zeros((im_h,im_w),dtype=xp.float32)
  rhs_phie           = xp.zeros((im_h,im_w),dtype=xp.float32)
  rhs_vmem           = xp.zeros((im_h,im_w),dtype=xp.float32)
  vmem               = xp.zeros((im_h,im_w),dtype=xp.float32)
  vmem               = cell_state[lr_params.index('v')]

  # PDE
  pde_phie_m = PDE( im_h, im_w, sigma_l,   sigma_t,   ds )
  pde_phie_i = PDE( im_h, im_w, sigma_l_i, sigma_t_i, ds )
  pde_vmem   = PDE( im_h, im_w, sigma_l_i, sigma_t_i, ds )

  t         = 0.                           # Time (ms)
  dt        = const_d['dt_']               # Time step (ms)
  st_train  = 1                            # Number of Beats
  st_start  = 0.005                        # Time to begin stim (ms)
  st_amp    = -10.0e-0                     # Stim amplitude (uA/cm^2)
  st_inter  = 250.                         # Basic Cycle Length (ms)
  st_dur    = 5.0e0                        # Stim duration (ms)
  st_time   = 0.0                          # Past time of stimulation (ms)
  st_on     = False                        # Stimulation flag 
  steps     = int((st_inter*st_train)/dt)  # Number of loop

  log_interval = 1000
  for step in range(steps):

    # Stimulation control
    if int( (t - st_start)/dt ) % int(st_inter/dt) == 0:
      st_on = True
      print 'stim on'
      #i_ext[im_h/2-2:im_h/2+2,im_w/2-2:im_w/2+2] = st_amp
      cell_state[lr_params.index('st')][im_h/2-2:im_h/2+2,im_w/2-2:im_w/2+2] = st_amp
    if st_on:
      st_time += dt 
      if st_time > st_dur:
        print 'stim off'
        st_on = False
        #i_ext[:,:] = 0.0 
        cell_state[lr_params.index('st')][:,:] = 0.0
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
      rhs_phie = i_ext - pde_vmem.forward(vmem) 
      #if step%log_interval==0:
      #  np.save('y', rhs_phie.get())
      #  np.save('x', phie.get())
      #  np.save('R', pde_phie_m.R.get())
      #  np.save('D', pde_phie_m.D.get())
      cnt, phie = pde_phie_m.solve(phie, rhs_phie, tol=1e-5)

      # step.3 vmem
      rhs_vmem = pde_vmem.forward(vmem)
      rhs_vmem += pde_phie_i.forward(phie) 
      rhs_vmem -= i_ion  * Sv
      rhs_vmem *= (1 / Cm / Sv)  
      vmem += dt * rhs_vmem

    if step%log_interval==0:
      print '------------------{0}ms'.format(t)
      #print 'i_ext'
      #print i_ext
      #print 'rhs_phie'
      #print rhs_phie
      #print 'phie'
      #print phie
      #print 'rhs_vmem'
      #print rhs_vmem
      #print 'vmem_pst'
      #print vmem
      yield phie.get(), vmem.get(), t

    t += dt

  print "elecpy done"

def draw(data):
  phie, vmem, t = data[0], data[1], data[2]
  im.set_array(vmem)
  return [im]

ani = animation.FuncAnimation(fig, draw, sim, blit=False, interval=1, repeat=False)
plt.show()
