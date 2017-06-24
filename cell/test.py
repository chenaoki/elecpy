import os, json, time
import numpy as np
from optparse import OptionParser

from luorudy.model import model as luorudy
from ohararudy.model import model as ohararudy
from mahajan.model import model as mahajan

def test(pacing_params, savepath):

    print 'cell models pacing test.'

    shape = (1,1)
    if not os.path.isdir(savepath) :
        os.mkdir(savepath)
    with open('{0}/pacing_params.json'.format(savepath), 'w') as f:
        json.dump(pacing_params, f, indent=4)

    dt           = pacing_params["time"]["dt"]
    log_cnt      = pacing_params["time"]["log_cnt"]
    st_start     = pacing_params["stim"]["start"]      # Time to begin stim (ms)
    st_dur       = pacing_params["stim"]["duration"]   # Stim duration (ms)
    st_inter     = pacing_params["stim"]["interval"]   # Cycle Length (ms)
    st_train     = pacing_params["stim"]["train"]      # Number of Beats
    st_amp       = pacing_params["stim"]["amplitude"]  # Stim amplitude (uA/cm^2)

    # Cell state initialization
    if pacing_params["cell_type"] == 'luorudy':
        model = luorudy(shape)
    if pacing_params["cell_type"] == 'ohararudy':
        model = ohararudy(shape)
    if pacing_params["cell_type"] == 'mahajan':
        model = mahajan(shape)

    model.set_param('dt', dt)
    for param in pacing_params["cell_param"].keys():
        value = pacing_params["cell_param"][param]
        model.set_param(param, value)

    # Initialization
    t = 0.                                   # Time (ms)
    st_time = 0.0                            # Past time of stimulation (ms)
    st_on = False                            # Stimulation flag
    steps = int((st_inter*st_train)/dt)      # Number of loop
    start = time.time()
    cnt = 0

    for step in range(steps):
    #for step in range(1):
    
        # Stimulation setting
        if int( (t - st_start)/dt ) % int(st_inter/dt) == 0:
            st_on = True
        if st_on:
            model.set_param('st', st_amp)
            st_time += dt
            if st_time > st_dur:
                st_on = False
                st_time = 0.0
        else: # st_on
            model.set_param('st', 0.)

        # State transition
        model.update()

        if step % log_cnt == 0:
            cnt += 1
            model.save(os.path.join(savepath, '{0:0>5}'.format(cnt)))
            if cnt % 100 == 0:
                print '-------------{0}'.format(cnt)
                print 'v:', model.get_param('v')[0,0]
                print 'it:',model.get_param('it')[0,0]
                print 'st:',model.get_param('st')[0,0]

        t += dt

    elapsed_time = time.time() - start
    print 'elapsed_time:', elapsed_time
    
if __name__ == "__main__":
    
    parser = OptionParser()
    parser.add_option(
        '-p','--param_file',
        dest='param_file', action='store', type='string', default='../temp/pacing_params.json',
        help="json file of simulation parameters")
    parser.add_option(
        '-d','--dst',
        dest='savepath', action='store', type='string', default='../temp/result/',
        help="Save data path.")

    (options, args) = parser.parse_args()

    with open (options.param_file,'r') as f:
        pacing_params = json.load(f)

    test(pacing_params, options.savepath)
