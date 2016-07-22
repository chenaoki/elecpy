# -*- coding:utf-8 -*-

import numpy as np

class Stimulator(object):

  def __init__(self, shape, size, start, interval, duration, amplitude, is_anode=True, train=1, name="noname"): 
    self.shape        = shape
    self.size         = size
    self.interval     = interval
    self.duration     = duration
    self.amplitude    = amplitude  # (uA/cm^2)
    self.is_anode     = is_anode
    self.start        = start
    self.train        = train
    self.name         = name
    self._map_on      = np.zeros( self.shape, dtype=np.float32 )
    self._map_off     = np.zeros( self.shape, dtype=np.float32 )
    self._map_anode   = np.zeros( self.shape, dtype=np.float32 )
    self._map_cathode = np.zeros( self.shape, dtype=np.float32 )

    if self.name == "area_L" : self.set_area_L(self.size) 
    if self.name == "area_R" : self.set_area_R(self.size) 
    if self.name == "area_T" : self.set_area_T(self.size) 
    if self.name == "area_B" : self.set_area_B(self.size) 
    if self.name == "edge_L" : self.set_edge_L(self.size) 
    if self.name == "edge_R" : self.set_edge_R(self.size) 
    if self.name == "edge_T" : self.set_edge_T(self.size) 
    if self.name == "edge_B" : self.set_edge_B(self.size) 
    if self.name == "point"  : self.set_point(self.size) 

  def get_current(self, t):
    assert self._map_on is not None
    assert self._map_off is not None
    if t >= self.start and t < self.start + self.train * self.interval and (t - self.start) % self.interval < self.duration:
      return self._map_on
    else:
      return self._map_off

  def set_location(self):
    anode   = (( self._map_anode > 0 )*1).astype(np.float32)
    cathode = (( self._map_cathode > 0 )*1).astype(np.float32)
    self._map_on *= 0.0 # reset
    if self.is_anode:
      self._map_on += anode*self.amplitude
      self._map_on += cathode*( -self.amplitude * np.sum(anode)/ float(np.sum(cathode)) ) 
    else:
      self._map_on += cathode*self.amplitude
      self._map_on += anode*( -self.amplitude * np.sum(cathode)/ float(np.sum(anode)) ) 

  def set_area_L(self, border):
    assert border > 0 and border < self.shape[1]
    self._map_anode[:,:border] = 1.0
    self._map_cathode[:,border:] = 1.0
    self.set_location()
  
  def set_area_R(self, border):
    assert border > 0 and border < self.shape[1]
    self._map_cathode[:,:border] = 1.0
    self._map_anode[:,border:] = 1.0
    self.set_location()

  def set_area_T(self, border):
    assert border > 0 and border < self.shape[0]
    self._map_anode[:border,:] = 1.0
    self._map_cathode[border:,:] = 1.0
    self.set_location()

  def set_area_B(self, border):
    assert border > 0 and border < self.shape[0]
    self._map_cathode[:border,:] = 1.0
    self._map_anode[border:,:] = 1.0
    self.set_location()

  def set_edge_L(self, width = 1):
    assert width > 0 and width < self.shape[1]
    self._map_anode[:,:width] = 1.0
    self._map_cathode[:,width:2*width] = 1.0
    self.set_location()

  def set_edge_R(self, width = 1):
    assert width > 0 and width < self.shape[1]
    self._map_anode[:,-width:] = 1.0
    self._map_cathode[:,-2*width:-width] = 1.0
    self.set_location()

  def set_edge_T(self, width = 1):
    assert width > 0 and width < self.shape[0]
    self._map_anode[:width,:] = 1.0
    self._map_cathode[width:2*width,:] = 1.0
    self.set_location()

  def set_edge_B(self, width = 1):
    assert width > 0 and width < self.shape[0]
    self._map_anode[-width:,:] = 1.0
    self._map_cathode[-2*width:-width,:] = 1.0
    self.set_location()

  def set_point( self, ( y, x, rad) ):
    assert x >= rad and x + rad <= self.shape[0]
    assert y >= rad and y + rad <= self.shape[1]
    if self.is_anode:
      self._map_anode[x-rad:x+rad,y-rad:y+rad] = 1.0
      self._map_cathode[0,:] = 1.0
      self._map_cathode[-1,:] = 1.0
    else:
      self._map_cathode[x-rad:x+rad,y-rad:y+rad] = 1.0
      self._map_anode[0,:] = 1.0
      self._map_anode[-1,:] = 1.0
    self.set_location()

if __name__ == '__main__':

  import matplotlib.pylab as plt
  from matplotlib import animation
  import json

  with open ('./sim_params.json','r') as f : sim_params = json.load(f)
  stim_params = sim_params['stimulation']
  assert len(stim_params) > 0
  stims = []
  for param in stim_params:
    stims.append(Stimulator(**param))

  im_h, im_w = stims[0].shape

  fig = plt.figure(figsize=(10,10))
  im = plt.imshow(
      np.zeros((im_h,im_w),dtype=np.float32), 
      vmin = -100.0, vmax = 100.0, 
      cmap='hot', 
      interpolation='nearest')
  plt.axis('off')

  def sim():
    t = 0.0
    t_end = 2000.0
    dt =10
    i_ext_e = np.zeros((im_h,im_w),dtype=np.float32)
    while t < t_end:
      i_ext_e[:,:] = 0.0
      for s in stims:
        i_ext_e += s.get_current(t)
      print t
      yield i_ext_e
      np.save('./result/stim_pattern/{0:0>4}'.format(int(t)), i_ext_e)
      t += dt
    exit()

  def draw(data):
    img = data
    im.set_array(img)
    return [im]

  ani = animation.FuncAnimation(fig, draw, sim, blit=False, interval=200, repeat=False)
  plt.show()
