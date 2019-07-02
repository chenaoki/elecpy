import numpy as np
import os, glob
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import h5py

from scipy.ndimage import convolve
from numba.decorators import autojit

class Loader(object):
    
    def __init__(self, path, keys=None, frames=None):
        
        self.path = path
        self.L = None
        self.shape = None
        self.data = {}
        
        if keys is None:
            self.keys = {'vmem', 'phie'}
        else:
            self.keys = keys

        with h5py.File( path, 'r') as f:
            for key in self.keys:
                self.data[key] = np.array([f[frame][key].value for frame in f.keys()])
                if frames is not None:
                    self.data[key] = self.data[key][frames,:,:]
                self.L = self.data[key].shape[0]
        pass

                
    def setRange(self, x_min = None, x_max = None, y_min = None, y_max = None, f_min=None, f_max=None):
        if x_min is None: x_min = 0
        if y_min is None: y_min = 0
        if y_min is None: f_min = 0
        if x_min is None: x_max = self.data['vmem'].shape[2]
        if y_min is None: y_max = self.data['vmem'].shape[1]
        if y_min is None: f_max = self.data['vmem'].shape[0]
        
        for key in self.data.keys():
            self.data[key] = self.data[key][f_min:f_max, y_min:y_max, x_min:x_max]
            
            
    def getNormalized(self):
            
        ret = Loader (self.path)
        
        def pixelwise_normalize(X):
            L, M, N = X.shape
            ret = np.zeros_like(X)
            for i in range(M):
                for j in range(N):
                    ts = X[:,i,j]
                    ret[:, i, j] = (ts-ts.min())/(abs(ts.max() - ts.min())+1.0e-30)
            return ret
                    
        _func = autojit(pixelwise_normalize)
        for key in self.data.keys():
            ret.data[key] = _func(self.data[key])
        
        return ret
    
    def saveAnimation(self, save_dir, keys=None, time_range=None, ext='mp4', cmap='gray'):
        
        if keys is None : keys = self.data.keys()
        if time_range is None : time_range = np.arange(self.L)
            
        for key in keys:

            fig = plt.figure()
            plt.axis('off')
            data = self.data[key]
            vmin = data.min()
            vmax = data.max()

            ims = []
            for i, img in enumerate(data[time_range]):
                im = plt.imshow(
                    img,
                    vmin = vmin, vmax = vmax,
                    cmap=cmap,
                    interpolation='nearest')
                ims.append([im])

            ani = animation.ArtistAnimation(fig, ims, interval=30)
            
            if ext is 'gif':
                ani.save(os.path.join(save_dir, '{0}.gif'.format(key)), writer="imagemagick")
            if ext is 'mp4':
                ani.save(os.path.join(save_dir, '{0}.mp4'.format(key)), writer="ffmpeg")
            
    def pseudoECG(self, x, y, z, key='vmem'):
        
        vmem = self.data[key]                
        L,M,N = vmem.shape

        def pixelwise_distance(x,y,z):
            dist = np.zeros((M, N))
            for i in range(M):
                for j in range(N):
                    dist[i, j] = np.sqrt( z**2 + (y - i)**2 + (x - j)**2 )
            return dist                    

        im_dist = autojit(pixelwise_distance)(x,y,z)
        R_y, R_x = np.gradient(1/im_dist)

        def framewise_ecg():
            ecg = np.zeros((L))
            for f in range(L):
                V_y, V_x = np.gradient(vmem[f])
                ecg[f] = -np.sum( R_x*V_x + R_y*V_y )
            return ecg
        
        ecg = autojit(framewise_ecg)()
        return ecg
                
