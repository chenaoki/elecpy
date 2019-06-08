import numpy as np
import os, glob
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from scipy.ndimage import convolve
from numba.decorators import autojit

class Loader(object):
    
    def __init__(self, path, keys=None):
        
        self.path = path
        self.L = None
        self.shape = None
        self.data = {}
        
        if keys is None:
            self.keys = {'vmem', 'phie', 'cell/m', 'cell/h', 'cell/j', 'cell/xina'}
        else:
            self.keys = keys
        
        for key in self.keys:
            
            if 'cell' in key:
                files = sorted( glob.glob(os.path.join(self.path, 'cell*/{0}.npy'.format(key.replace('cell', '')))))
            else:
                files = sorted( glob.glob(os.path.join(self.path, '{0}*.npy'.format(key))))
            
            if self.L is None:
                self.L = len(files)
            else:
                assert self.L == len(files)
                
            img = np.load( files[0])
            if 'cell' in key:
                img = img.reshape(self.shape)
                
            if self.shape is None:
                self.shape = img.shape
            else:
                assert self.shape == img.shape
        
            self.data[key] = np.zeros( np.concatenate(([self.L], img.shape)), dtype=img.dtype)            
            for i, _f in enumerate(files): 
                img = np.load( _f)
                if 'cell' in key:
                    img = img.reshape(self.shape)
                self.data[key][i,:,:] = img
                
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
            
    def pseudoECG(self, i, j):
        
        eta = np.ones((3,3), dtype=np.float64)
        
        k = np.ones((3,3), dtype=np.float64)
        k[1,1] = -8.
        
        k = k*eta
        
        phie = self.data['phie']                
        L,M,N = phie.shape

        def pixelwise_distance(pos_y, pos_x, pos_z = 0.):
            dist = np.zeros((M, N))
            for i_ in range(M):
                for j_ in range(N):
                    dist[i_, j_] = np.sqrt( pos_z**2 + (pos_y - i_)**2 + (pos_x -  j_)**2 )
            return dist                    
        im_dist = autojit(pixelwise_distance)(i, j)
        im_dist[i, j] = 1.0 ## to avoid zero division

        def framewise_ecg():
            ecg = np.zeros((L))
            for f in range(L):
                im_conv = convolve(phie[f,:,:], k, mode='constant', cval=0.)
                im_mul = im_conv/im_dist
                im_mul[i, j] = 0.0
                ecg[f] = np.sum(im_mul)
            return ecg
        ecg = autojit(framewise_ecg)()
        return ecg
                
