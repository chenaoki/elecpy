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
  assert len(files) is not 0
  
  img_sample = np.load(files[0])
  im_h, im_w = img_sample.shape
  assert im_h > 0 and im_w > 0
  V = np.zeros((len(files), im_h, im_w), dtype=np.float32)
  for i, path in enumerate(files) : V[i, :, :] = np.load(path)
  vmax = np.max(V)
  vmin = np.min(V)

  def yield_data():
    for i in range(V.shape[0]):
      print files[i]
      yield V[i, :, :]
    exit()

  def yield_frame():
    from Getch import _Getch
    getch = _Getch()
    frame=0
    while True:
      print 'frame {0:0>4}'.format(frame)
      yield V[frame,:,:]
      #c = raw_input('>>>' )
      c = str(getch())
      print c
      if c is 'q':
        break 
      if c is 'j':
        frame -= 1
      if c is 'k':
        frame += 1
      if c is 'h':
        frame -= 5
      if c is 'l':
        frame += 5
      frame = frame if frame >= 0 else 0
      frame = frame if frame < V.shape[0] else V.shape[0]-1

  fig = plt.figure(figsize=(5,5))
  im = plt.imshow(
    np.zeros((im_h,im_w),dtype=np.float32), 
    vmin = vmin, vmax = vmax, 
    cmap=bipolar(neutral=0, lutsize=1024), 
    interpolation='nearest')
  plt.axis('off')

  def draw(data):
    img = data
    im.set_array(img)
    return [im]

  if options.visualize:
    ani = animation.FuncAnimation(fig, draw, yield_data, blit=False, interval=1, repeat=False)
    #ani = animation.FuncAnimation(fig, draw, yield_frame, blit=False, interval=1, repeat=False)
    plt.show()
  else:
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=30, metadata=dict(artist='Elecpy'), bitrate=1800)
    ani.save('{0}/out.mp4'.format(options.srcpath), writer=writer) 
