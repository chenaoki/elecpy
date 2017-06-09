import matplotlib.pylab as plt
from matplotlib import animation
import json
import numpy as np
from optparse import OptionParser
from MembraneStimulator import MembraneStimulator
from ExtracellularStimulator import ExtracellularStimulator

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option(
        '-p','--param_file', 
        dest='param_file', action='store', type='string', default='./sim_params_test.json',
        help="json file of simulation parameters")
    (options, args) = parser.parse_args()
    print options.param_file
    
    with open (options.param_file,'r') as f : sim_params = json.load(f)
    im_h      = sim_params['geometory']['height']
    im_w      = sim_params['geometory']['width']
    t_end     = sim_params['time']['end']

    stim_params = sim_params['stimulation']
    assert len(stim_params) > 0
    stims_ext = []
    stims_mem = []
    for param in stim_params:
        stims_ext.append(ExtracellularStimulator(**param))
        stims_mem.append(MembraneStimulator(**param))

    fig = plt.figure(figsize=(10,10))
    im = plt.imshow(
        np.zeros((im_h,im_w),dtype=np.float32), 
        vmin = -200.0, vmax = 200.0, 
        cmap='hot', 
        interpolation='nearest')
    plt.axis('off')

    def sim():
        i_ext = np.zeros((im_h,im_w),dtype=np.float32)
        i_mem = np.zeros((im_h,im_w),dtype=np.float32)

        t = 0.0
        while t < t_end:
            i_ext[:,:] = 0.0
            for s in stims_ext:
                i_ext += s.get_current(t)
            print t
            yield i_ext
            np.save('./stim_pattern/ext_{0:0>4}'.format(int(t)), i_ext)
            t += 1
        
        t = 0.0
        while t < t_end:
            i_mem[:,:] = 0.0
            for s in stims_ext:
                i_ext += s.get_current(t)
            print t
            yield i_ext
            np.save('./stim_pattern/mem_{0:0>4}'.format(int(t)), i_ext)
            t += 1
            
        exit()

    def draw(data):
        im.set_array(data)
        return [im]

    ani = animation.FuncAnimation(fig, draw, sim, blit=False, interval=200, repeat=False)
    plt.show()