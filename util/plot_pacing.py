import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import opmap
from opmap.VideoData import VideoData

from elecpy.cell.mahajan.model import model as model_mj
from elecpy.cell.ohararudy.model import model as model_ord
from elecpy.cell.luorudy.model import model as model_lrd

def get_pacing_series(dir_path, cell_type):

    #dir_path = './elecpy/temp/result/'
    #cell_type = 'mahajan'
    #param_name = 'v'

    dirs = glob.glob(dir_path+'*/')[:-10]
    shape =  [long(len(dirs))]
    shape.extend(np.load(os.path.join( dirs[0], 'v.npy')).shape)
    shape = tuple(shape)

    if cell_type == 'mahajan':
        params = model_mj(0).params
    if cell_type == 'luorudy':
        params = model_lrd(0).params
    if cell_type == 'ohararudy':
        params = model_ord(0).params

    series = {}
    for p in params:
        series[p] = np.zeros(shape, np.float64)

        for i, d in enumerate(dirs):
            f = os.path.join(d, p+'.npy')
            series[p][i, :, :] = np.load(f)
    return series

def save_pacing_plot(series, save_path):
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    for p in series.keys():
        plt.clf()
        plt.plot(series[p][:,0,0])
        plt.savefig(os.path.join(save_path, '{0}.png'.format(p)))
