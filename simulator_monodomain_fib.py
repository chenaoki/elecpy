#!/usr/local/bin/python

import os, sys
import json
import numpy as np
import matplotlib.pyplot as plt
import chainer
import h5py
from chainer import cuda
from matplotlib import animation
from optparse import OptionParser

from solver.PDE import PDE
from stim.MembraneStimulator import MembraneStimulator
from cell.ohararudy.model import model as cell_model_ohararudy
from cell.luorudy.model import model as cell_model_luorudy
from cell.mahajan.model import model as cell_model_mahajan
from cell.maccannell.model import model as cell_model_maccannell
from util.cmap_bipolar import bipolar

# global variables
class MonodomainWithFibroblastSimulator(object):

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

        print "elecpy simulation start!"

        cuda.get_device(0).use()

        # Constants
        Sv           = 1400                  # Surface-to-volume ratio (cm^-1)
        Cm           = 1.0                   # Membrane capacitance (uF/cm^2)
        sigma_l_i    = 1.74                  # (mS/cm)
        sigma_t_i    = 0.225                 # (mS/cm)
        sigma_mu     = 1.74/6.25             # equal anisotropy (sigma_i = mu * sigma_e) 

        # Myocyte <=> Fibroblast settings
        S_mf         = 6.3e-6                # fibroblast surface area (cm^2)
        S_m          = 1.0e-4                # myocyte membrane surface area (cm^2)
        C_mf         = 1.0                   # (myo)fibroblast membrane capacitance (uF/cm^2)
        G_gap        = 3.0e-6                # gap junction conductance between myocyte and (myo)fibroblast (mS)

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
        assert cells is not None

        # fibroblast cells
        cells_fib = cell_model_maccannell((N)) 

        print "Stimulation settings",
        stims_ext = []
        stims_mem = []
        if 'stimulation' in sim_params.keys():
            stim_param = sim_params['stimulation']
            if 'membrane' in stim_param:
                for param in stim_param['membrane']:
                    stim = MembraneStimulator(**param)
                    assert tuple(stim.shape) == (im_h, im_w)
                    stims_mem.append(stim)
        print "...done"

        print "Allocating data...",
        cells.create()
        v_myo              = np.copy(cells.get_param('v'))
        v_fib              = np.copy(cells_fib.get_param('v'))
        i_ion              = np.zeros((N),dtype=np.float64) # ionic current of myocyte
        i_ion_fib          = np.zeros((N),dtype=np.float64) # ionic current of fibroblast
        i_ion_fm           = np.zeros((N),dtype=np.float64) # ionic current from fibroblasts into myocyte
        i_ext_e            = np.zeros((N),dtype=np.float64)
        i_ext_i            = np.zeros((N),dtype=np.float64)
        dist_fib           = np.zeros((N),dtype=np.float64)
        rhs_v_myo          = np.zeros((N),dtype=np.float64)
        print "...done"

        print "Initializing data...",
        if 'restart' in sim_params.keys():
            cnt_restart = sim_params['restart']['count']
            srcpath = sim_params['restart']['source']
            with h5py.File(os.path.join(srcpath, 'out.h5'), 'r') as f:
                group_id = '{0:0>4}'.format(cnt_restart)
                v_myo = f[group_id]['vmem'].value.flatten()
                if 'vfib' in f[group_id].keys():
                    v_fib = f[group_id]['vfib'].value.flatten()
                cells.load(f, group_id)
            cnt_udt = cnt_restart * cnt_log
        print "...done"

        print "Mask settings...",
        if 'mask' in sim_params.keys():
            mask_param = sim_params['mask']
            for key in mask_param.keys():
                array = np.load(mask_param[key])
                assert array.shape == (im_h, im_w)
                cells.set_param(key, array.flatten())
        print "...done"

        print "Fibroblast settings...",
        assert 'fibroblast' in sim_params.keys()
        fib_params = sim_params['fibroblast']
        assert 'type' in fib_params.keys() 
        fib_type = fib_params['type']
        assert fib_type in ('fib' ,'mfib')
        fib_scale = {'fib':1.0, 'mfib':8.0}[fib_type] 
        assert 'distribution' in fib_params.keys() 
        array = np.load(fib_params['distribution'])
        assert array.shape == (im_h, im_w)
        dist_fib  = np.copy(array).flatten()
        print "...done"

        print 'Building PDE system ...',
        pde_i = PDE( im_h, im_w, sigma_l_i, sigma_t_i, ds )
        print '...done'

        # Initialization
        t         = 0.                       # Time (ms)
        cnt_udt   = 0                        # Count of udt
        dstep     = 1                        # Time step (# of udt)
        cnt_save  = -1

        print 'Main loop start!'
        with h5py.File(os.path.join(savepath, 'out.h5'),'w') as outf:
            
            while t < time_end:

                t = self.conv_cntUdt2time(cnt_udt)
                dt = dstep * udt

                # Stimulation control
                i_ext_e[:] = 0.0
                flg_st_temp = False
                for s in stims_ext:
                    i_ext_e += s.get_current(t)*Sv
                    flg_st_temp = flg_st_temp or s.get_flag(t)
                    
                array_mem = np.zeros(im_h*im_w)
                for s in stims_mem:
                    #cells.set_param('st', s.get_current(t)) 
                    array_mem += s.get_current(t)
                mask_mem = (np.abs(array_mem) > 0.0)*1.0
                v_myo = mask_mem * array_mem + (1.0-mask_mem)*v_myo

                # step.1 fibroblast update
                cells_fib.set_param('dt', dt)
                cells_fib.set_param('v', cuda.to_gpu(v_fib) )
                cells_fib.update()
                i_ion_fib = cells_fib.get_param('it')
                dv_fib = -dt*(1.0/(C_mf*fib_scale))*(i_ion_fib + G_gap*(v_fib - v_myo)/S_mf)
                v_fib += dv_fib
                i_ion_fm = dist_fib*G_gap*(v_myo - v_fib)/S_m

                # step.2 myocyte state transition
                cells.set_param('dt', dt )
                cells.set_param('v', cuda.to_gpu(v_myo) )
                cells.update()
                i_ion = cells.get_param('it')

                # step.3 v_myo
                rhs_v_myo = pde_i.forward(v_myo)
                rhs_v_myo += i_ext_i
                rhs_v_myo += sigma_mu * i_ext_e
                rhs_v_myo *= 1 / (1+sigma_mu)
                rhs_v_myo -= ( i_ion + i_ion_fm )* Sv
                rhs_v_myo *= 1 / (Cm * Sv)
                v_myo += dt * rhs_v_myo

                # Logging & error check
                cnt_save_now = self.conv_time2cntSave(t)
                if cnt_save_now != cnt_save:
                    cnt_save = cnt_save_now
                    sys.stdout.write('\r------------------{0}/{1}ms'.format(t, time_end))
                    sys.stdout.flush()

                    group_id = '{0:0>4}'.format(cnt_save)
                    outf.create_group(group_id)
                    outf[group_id].create_dataset('vmem', data = v_myo.reshape((im_h, im_w)))
                    outf[group_id].create_dataset('vfib', data = v_fib.reshape((im_h, im_w)))
                    cells.save(outf, group_id)
                    yield v_myo

                    flg = False
                    for i,v in enumerate(v_myo):
                        if v != v :
                            print "error : invalid value {1} @ {0} ms, index {2}".format(t, v, i)
                            flg = True
                            break
                    if flg is True:
                        break

                cnt_udt += dstep

            print "elecpy done"
            yield False

