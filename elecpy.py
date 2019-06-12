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
from stim.ExtracellularStimulator import ExtracellularStimulator
from stim.MembraneStimulator import MembraneStimulator
from cell.ohararudy.model import model as cell_model_ohararudy
from cell.luorudy.model import model as cell_model_luorudy
from cell.mahajan.model import model as cell_model_mahajan
from util.cmap_bipolar import bipolar

# global variables
sim_params = None
cells = None
stims_ext = []
stims_mem = []
i_ion    = None
phie     = None
i_ext_e  = None
i_ext_i  = None
rhs_phie = None 
rhs_vmem = None 
vmem     = None 

# Functions
def conv_cntSave2time(cnt_save):
    global sim_params
    udt          = sim_params['time']['udt']     # Universal time step (ms)
    cnt_log      = sim_params['log']['cnt']      # num of udt for logging
    return udt*cnt_log

def conv_cntUdt2time(cnt_udt):
    global sim_params
    udt          = sim_params['time']['udt']     # Universal time step (ms)
    return cnt_udt * udt

def conv_time2cntUdt(t):
    global sim_params
    udt          = sim_params['time']['udt']     # Universal time step (ms)
    return int(t/udt)

def conv_time2cntSave(t):
    global sim_params
    udt          = sim_params['time']['udt']     # Universal time step (ms)
    cnt_log      = sim_params['log']['cnt']      # num of udt for logging
    return conv_time2cntUdt(t) // cnt_log 

def sim_generator( params ):

    global sim_params, cells, stims_ext, stims_mem
    global i_ion, phie, i_ext_e, i_ext_i, rhs_phie, rhs_vmem, vmem
    
    sim_params = params

    assert sim_params is not None

    print "elecpy simulation start!"

    cuda.get_device(0).use()

    # Constants
    Sv           = 1400                  # Surface-to-volume ratio (cm^-1)
    Cm           = 1.0                   # Membrane capacitance (uF/cm^2)
    sigma_l_i    = 1.74                  # (mS/cm)
    sigma_t_i    = 0.19                  # (mS/cm)
    sigma_l_e    = 6.25                  # (mS/cm)
    sigma_t_e    = 2.36                  # (mS/cm)

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
    cells.create()

    print "Stimulation settings",
    stims_ext = []
    stims_mem = []
    if 'stimulation' in sim_params.keys():
        stim_param = sim_params['stimulation']
        if 'extracellular' in stim_param:
            for param in stim_param['extracellular']:
                stim = ExtracellularStimulator(**param)
                assert tuple(stim.shape) == (im_h, im_w)
                stims_ext.append(stim)
        if 'membrane' in stim_param:
            for param in stim_param['membrane']:
                stim = MembraneStimulator(**param)
                assert tuple(stim.shape) == (im_h, im_w)
                stims_mem.append(stim)
    print "...done"
    

    print "Allocating data...",
    i_ion              = np.zeros((N),dtype=np.float64)
    phie               = np.zeros((N),dtype=np.float64)
    i_ext_e            = np.zeros((N),dtype=np.float64)
    i_ext_i            = np.zeros((N),dtype=np.float64)
    rhs_phie           = np.zeros((N),dtype=np.float64)
    rhs_vmem           = np.zeros((N),dtype=np.float64)
    tone               = np.zeros((N),dtype=np.float64)
    vmem               = np.copy(cells.get_param('v'))
    print "...done"

    print "Initializing data...",
    if 'restart' in sim_params.keys():
        cnt_restart = sim_params['restart']['count']
        srcpath = sim_params['restart']['source']
        with h5py.File(os.path.join(srcpath, 'out.h5'), 'r') as f:
            group_id = '{0:0>4}'.format(cnt_restart)
            phie = f[group_id]['phie'].value.flatten()
            vmem = f[group_id]['vmem'].value.flatten()
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

    print 'Building PDE system ...',
    sigma_l      = sigma_l_e + sigma_l_i
    sigma_t      = sigma_t_e + sigma_t_i
    pde_i = PDE( im_h, im_w, sigma_l_i, sigma_t_i, ds )
    pde_m = PDE( im_h, im_w, sigma_l,   sigma_t,   ds )
    print '...done'

    # Initialization
    t         = 0.                       # Time (ms)
    cnt_udt   = 0                        # Count of udt
    dstep     = 1                        # Time step (# of udt)
    cnt_save  = -1
    
    run_udt   = True                     # Flag of running sim in udt
    flg_st    = False                    # Flaf of stimulation
    cnt_st_off = 0

    print 'Main loop start!'
    with h5py.File(os.path.join(savepath, 'out.h5'),'w') as outf:
        while t < time_end:

            t = conv_cntUdt2time(cnt_udt)
            dt = dstep * udt

            # Stimulation control
            i_ext_e[:] = 0.0
            flg_st_temp = False
            for s in stims_ext:
                i_ext_e += s.get_current(t)*Sv
                flg_st_temp = flg_st_temp or s.get_flag(t)
            for s in stims_mem:
                cells.set_param('st', s.get_current(t)) 

            # step.1 cell state transition
            cells.set_param('dt', dt )
            cells.set_param('v', cuda.to_gpu(vmem) )
            cells.update()
            i_ion = cells.get_param('it')


            # step.2 phie
            rhs_phie = i_ext_e - i_ext_i - pde_i.forward(vmem)
            pde_cnt, phie = pde_m.solve(phie, rhs_phie, tol=1e-2, maxcnt=1e5)
            phie -= phie[0]

            # step.3 vmem
            rhs_vmem = pde_i.forward(vmem)
            rhs_vmem += pde_i.forward(phie)
            tone     = ( rhs_vmem * dt ) / (Cm * Sv)
            rhs_vmem -= i_ion * Sv
            rhs_vmem += i_ext_i
            rhs_vmem *= 1 / (Cm * Sv)
            vmem += dt * rhs_vmem

            # Logging & error check
            cnt_save_now = conv_time2cntSave(t)
            if cnt_save_now != cnt_save:
                cnt_save = cnt_save_now
                sys.stdout.write('\r------------------{0}/{1}ms'.format(t, time_end))
                sys.stdout.flush()

                group_id = '{0:0>4}'.format(cnt_save)
                outf.create_group(group_id)
                outf[group_id].create_dataset('vmem', data = vmem.reshape((im_h, im_w)))
                outf[group_id].create_dataset('phie', data = phie.reshape((im_h, im_w)))
                outf[group_id].create_dataset('tone', data = tone.reshape((im_h, im_w)))
                cells.save(outf, group_id)
                yield vmem

                flg = False
                for i,v in enumerate(vmem):
                    if v != v :
                        print "error : invalid value {1} @ {0} ms, index {2}".format(t, v, i)
                        flg = True
                        break
                if flg is True:
                    break

            cnt_udt += dstep

    print "elecpy done"
    yield False

