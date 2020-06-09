#!/usr/local/bin/python

import os, sys
import json
import numpy as np
import matplotlib.pyplot as plt
import chainer
import h5py
import time
from chainer import cuda
from matplotlib import animation
from optparse import OptionParser

from .solver.PDE_dev import PDE
#from .solver.PDE import PDE
from .stim.MembraneStimulator import MembraneStimulator
from .cell.ohararudy.model import model as cell_model_ohararudy
from .cell.luorudy.model import model as cell_model_luorudy
from .cell.mahajan.model import model as cell_model_mahajan
from .cell.courtmanche.model import model as cell_model_courtmanche
from .util.cmap_bipolar import bipolar

# global variables
class MonodomainSimulator(object):

    def __init__(self, sim_params):
        self.sim_params = sim_params 

    def conv_cntUdt2time(self, cnt_udt):
        udt          = self.sim_params['time']['udt']     # Universal time step (ms)
        return cnt_udt * udt

    def conv_time2cntUdt(self, t):
        udt          = self.sim_params['time']['udt']     # Universal time step (ms)
        return int(t/udt)

    def conv_time2cntSave(self, t):
        udt          = self.sim_params['time']['udt']     # Universal time step (ms)
        cnt_log      = self.sim_params['log']['cnt']      # num of udt for logging
        return self.conv_time2cntUdt(t) // cnt_log 

    def step(self):

        sim_params = self.sim_params
        assert sim_params is not None

        #print("elecpy simulation start!")

        if 'gpu_id' in sim_params.keys():
            cuda.get_device(sim_params['gpu_id']).use()
        else:
            cuda.get_device(0).use()

        # Constants
        Sv           = 1400                  # Surface-to-volume ratio (cm^-1)
        Cm           = 1.0                   # Membrane capacitance (uF/cm^2)
        sigma_l_i    = 1.74                  # (mS/cm)
        sigma_t_i    = 0.225                 # (mS/cm)
        sigma_mu     = 1.74/6.25             # equal anisotropy (sigma_i = mu * sigma_e) 

        # Geometory settings
        im_h         = sim_params['geometory']['height']
        im_w         = sim_params['geometory']['width']
        ds           = sim_params['geometory']['ds'] # Spatial discretization step (cm)
        N            = im_h*im_w

        # Time settings
        udt          = sim_params['time']['udt']     # Universal time step (ms)
        time_end     = sim_params['time']['end']

        # Logging settings
        cnt_log      = sim_params['log']['cnt']      # num of udt for logging
        savepath     = sim_params['log']['path']
        save_all     = sim_params['log']['save_all']
    
        # Create result folder
        if not os.path.isdir(savepath) :
            os.mkdir(savepath)
        with open('{0}/sim_params.json'.format(savepath), 'w') as f:
            json.dump(sim_params, f, indent=4)

        # Cell model settings
        if sim_params['cell_type'] == 'ohararudy':
            cells = cell_model_ohararudy((N))
        if sim_params['cell_type'] == 'luorudy':
            cells = cell_model_luorudy((N))
        if sim_params['cell_type'] == 'mahajan':
            cells = cell_model_mahajan((N))
        if sim_params['cell_type'] == 'courtmanche':
            cells = cell_model_courtmanche((N))
        assert cells is not None

        #print("Stimulation settings")
        stims_ext = []
        stims_mem = []
        if 'stimulation' in sim_params.keys():
            stim_param = sim_params['stimulation']
            if 'membrane' in stim_param:
                for param in stim_param['membrane']:
                    stim = MembraneStimulator(**param)
                    assert tuple(stim.shape) == (im_h, im_w)
                    stims_mem.append(stim)
        #print("...done")

        #print("Allocating data...")
        cells.create()
        i_ion              = np.zeros((N),dtype=np.float64)
        phie               = np.zeros((N),dtype=np.float64)
        i_ext_e            = np.zeros((N),dtype=np.float64)
        i_ext_i            = np.zeros((N),dtype=np.float64)
        rhs_vmem           = np.zeros((N),dtype=np.float64)
        vmem               = np.copy(cells.get_param('v'))
        #print("...done")

        #print("Initializing data...")
        if 'restart' in sim_params.keys():
            cnt_restart = sim_params['restart']['count']
            srcpath = sim_params['restart']['source']
            with h5py.File(os.path.join(srcpath, 'out.h5'), 'r') as f:
                group_id = '{0:0>4}'.format(cnt_restart)
                vmem = f[group_id]['vmem'][()].flatten()
                cells.load(f, group_id)
            cnt_udt = cnt_restart * cnt_log
        #print("...done")

        #print("Mask settings...")
        if 'mask' in sim_params.keys():
            mask_param = sim_params['mask']
            for key in mask_param.keys():
                array = np.load(mask_param[key])
                assert array.shape == (im_h, im_w)
                cells.set_param(key, array.flatten())
        #print("...done")

        #print("Building PDE system ...")
        N_all = (im_h+2)*(im_w+2)
        sigma_l_i_array = np.ones((N_all),dtype=np.float64)*sigma_l_i
        sigma_t_i_array = np.ones((N_all),dtype=np.float64)*sigma_t_i
        if 'fiber' in sim_params.keys():
            fiber_param = sim_params['fiber']
            fiber_angle = np.load(fiber_param).flatten()
        else:
            fiber_angle = np.ones((N_all),dtype=np.float64)*0

        sw_it = cells.get_param('sw_it')
        thickness = cells.get_param('thickness')

        vmem = vmem*(sw_it==1) + np.ones_like(vmem)*(sw_it==0)*(-86.24)

        pde_i = PDE( im_h, im_w, sigma_l_i_array, sigma_t_i_array, ds, fiber_angle, sw_it )
        #pde_i = PDE( im_h, im_w, sigma_l_i, sigma_t_i, ds )
        #print("...done")

        # Initialization
        t         = 0.                       # Time (ms)
        cnt_udt   = 0                        # Count of udt
        dstep     = 1                        # Time step (# of udt)
        cnt_save  = -1

        #print("Main loop start!")
        with h5py.File(os.path.join(savepath, 'out.h5'),'w') as outf:

            while t < time_end:

                t = self.conv_cntUdt2time(cnt_udt)
                dt = dstep * udt
                cycle_num = dt / 0.001

                # Stimulation control
                i_ext_e[:] = 0.0
                flg_st_temp = False
                for s in stims_ext:
                    i_ext_e += s.get_current(t)*Sv/thickness
                    flg_st_temp = flg_st_temp or s.get_flag(t)
                
                array_mem = np.zeros(im_h*im_w)
                for s in stims_mem:
                    #cells.set_param('st', s.get_current(t)) 
                    array_mem += s.get_current(t)
                mask_mem = (np.abs(array_mem) > 0.0)*1.0
                vmem = mask_mem * array_mem + (1.0-mask_mem)*vmem

                # step.1 cell state transition
                cells.set_param('dt', dt / cycle_num)
                cells.set_param('v', cuda.to_gpu(vmem, device=sim_params['gpu_id']) )
                for i in range(int(cycle_num)):
                    cells.update()
                i_ion = cells.get_param('it')

                # step.3 vmem
                rhs_vmem = pde_i.forward(vmem) / thickness
                rhs_vmem += i_ext_i
                rhs_vmem += sigma_mu * i_ext_e
                rhs_vmem *= 1 / (1+sigma_mu)
                rhs_vmem -= i_ion * Sv / thickness
                rhs_vmem *= 1 / (Cm * Sv)
                vmem += dt * rhs_vmem * sw_it

                # Logging & error check
                cnt_save_now = self.conv_time2cntSave(t)
                if cnt_save_now != cnt_save:
                    cnt_save = cnt_save_now
                    sys.stdout.write('\r------------------{0}/{1}ms'.format(t, time_end))
                    sys.stdout.flush()

                    group_id = '{0:0>4}'.format(int(cnt_save))
                    outf.create_group(group_id)
                    outf[group_id].create_dataset('vmem', data = vmem.reshape((im_h, im_w)))
                    outf[group_id].create_dataset('phie', data = phie.reshape((im_h, im_w)))
                    cells.save(outf, group_id)
                    if save_all is True:
                        cells.save(outf, group_id)
                    yield vmem

                    flg = False
                    for i,v in enumerate(vmem):
                        if v != v :
                            print("error : invalid value {1} @ {0} ms, index {2}".format(t, v, i))
                            flg = True
                            break
                    if flg is True:
                        break

                cnt_udt += dstep

            #print("elecpy done")
            yield False

