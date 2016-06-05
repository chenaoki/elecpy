import numpy as np
import chainer
from chainer import cuda, Function, FunctionSet, Variable
import chainer.functions as F

from chainer.functions.connection.convolution_2d import convolution_2d

import matplotlib.pyplot as plt
from matplotlib import animation

''' 
dudt = a * Laplacian(u) + u - u^3 - v + k
tau * dvdt = b * Laplacian(v) + u - v
'''

dx = 0.02
dt = 0.001
a = 3e-4
b = 5e-3
k = 0.005
tau = 0.1
dif_u = a/dx/dx
dif_v = b/dx/dx/tau

class Laplacian:
  def __init__(self):
    self.W = Variable(xp.asarray([[[[0,1,0],[1,-4,1],[0,1,0]]]],dtype=np.float32), volatile=True)
    self.b = Variable(xp.zeros((1,), dtype=np.float32), volatile=True)
  def forward(self, x):
    return convolution_2d(x, self.W, self.b, pad=1)


class Reaction(Function):
  def forward(self, x):
    u = x[0]
    v = x[1]
    du, dv = cuda.elementwise(
        'T u, T v, T k, T tau', 'T du, T dv',
        '''
          du = u - u*u*u - v - k;
          dv = (u -v)/tau;
        ''',
        'reaction')(u,v,k,tau)
    return du, dv

class Periodic(Function):
  def forward(self, x):
    x = x[0]
    w = x.shape[2]
    h = x.shape[3]
    x[0,0,0,:] = x[0,0,w-2,:]
    x[0,0,w-1,:] = x[0,0,1,:]
    x[0,0,:,0] = x[0,0,:,h-2]
    x[0,0,:,h-1] = x[0,0,:,1]
    return x,

xp = cuda.cupy
cuda.get_device(0).use()

print "reac_diff start"

lap = Laplacian()

fig = plt.figure(figsize=(10,10))
im = plt.imshow(np.zeros((250,250),dtype=np.float32), vmin = -0.6, vmax = 0.6, cmap=plt.cm.gray)
plt.axis('off')

def sim():
  u = np.random.uniform(-0.2, 0.2, (1,1,250,250))
  v = np.random.uniform(-0.2, 0.2, (1,1,250,250))
  u = Variable(xp.asarray(u, dtype=np.float32), volatile=True)
  v = Variable(xp.asarray(v, dtype=np.float32), volatile=True)

  for i in range(10000):
    du, dv = Reaction()(u,v)
    du += dif_u * lap.forward(u)
    dv += dif_v * lap.forward(v)
    u_ = u+0.5*dt*du
    v_ = v+0.5*dt*dv
    u_ = Periodic()(u_)
    v_ = Periodic()(v_)

    du, dv = Reaction()(u_,v_)
    du += dif_u * lap.forward(u_)
    dv += dif_v * lap.forward(v_)
    u = u+dt*du
    v = v+dt*dv
    u = Periodic()(u)
    v = Periodic()(v)

    if i%100==0:
      yield u, v, i

def draw(data):
  u, v, i = data[0], data[1], data[2]
  print i
  im.set_array(u.data.get()[0,0,:,:])
  return [im]

ani = animation.FuncAnimation(fig, draw, sim, blit=False, interval=1, repeat=False)
plt.show()
print "reac_diff done"
