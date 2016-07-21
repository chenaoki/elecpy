#!/usr/local/bin/python

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
import glob
from optparse import OptionParser
from cmap_bipolar import bipolar

if __name__ == '__main__':

  parser = OptionParser()
  parser.add_option(
    '-s','--src', 
    dest='srcpath', action='store', type='string', default='./result/data', 
    help="Source data path [default : %default]")
  parser.add_option(
    '-p','--prefix', 
    dest='prefix', action='store', type='string', default='vmem_', 
    help="File prefix [default : %default]")
  parser.add_option(
    '-v','--visualize', 
    dest='visualize', action='store_true', 
    help="Visualize animation")
  (options, args) = parser.parse_args()

  files = glob.glob(options.srcpath+'/'+options.prefix+'*.npy')
  files = sorted(files)
  assert len(files) is not 0
  
  img_sample = np.load(files[0])
  im_h, im_w = img_sample.shape
  assert im_h > 0 and im_w > 0
  V = np.zeros((len(files), im_h, im_w), dtype=np.float32)
  for i, path in enumerate(files) : V[i, :, :] = np.load(path)
  vmax = np.max(V)
  vmin = np.min(V)

  if options.visualize:

    fig = plt.figure(figsize=(5,5))
    im = plt.imshow(
      np.zeros((im_h,im_w),dtype=np.float32), 
      vmin = vmin, vmax = vmax, 
      cmap=bipolar(neutral=0, lutsize=1024), 
      interpolation='nearest')
    plt.axis('off')

    def yield_data():
      for i in range(V.shape[0]):
        print files[i]
        yield V[i, :, :]
      exit()

    def draw(data):
      img = data
      im.set_array(img)
      return [im]

    ani = animation.FuncAnimation(fig, draw, yield_data, blit=False, interval=1, repeat=False)
    plt.show()

  else: # save mode

    matplotlib.use('Agg')

    for i in range(V.shape[0]):
      fig = plt.gcf()
      savepath = files[i].replace('.npy', '.jpg').replace('vmem', 'disp')
      print savepath
      plt.clf()
      plt.imshow(
        V[i,:,:], 
        vmin = vmin, vmax = vmax, 
        cmap=bipolar(neutral=0, lutsize=1024), 
        interpolation='nearest')
      #plt.show()
      plt.savefig(savepath)
