## Import
import numpy as np
import chainer
from chainer import cuda, Function, FunctionSet, Variable
import chainer.functions as F

from chainer.functions.connection.convolution_2d import convolution_2d

xp = cuda.cupy

## Functions

class Parabolic:
  def __init__(self, xx, yy, ds):
    self.W = Variable(xp.asarray(
      [[[
        [0.0,        yy/(ds**2),           0.0       ],
        [xx/(ds**2), -(2*xx+2*yy)/(ds**2), xx/(ds**2)],
        [0.0,        yy/(ds**2),           0.0       ]
      ]]]
      , dtype=np.float32), volatile=True)
    self.b = Variable( xp.zeros((1,), dtype=np.float32), volatile=True)
  def forward(self, x):
    return convolution_2d( x, self.W, self.b, pad=0 ) 

class IsolatedBoundary(Function):
  def forward(self,x):
    org   = x[0]
    #dst   = x[1]
    w     = org.shape[2]
    h     = org.shape[3]
    dst   = xp.asarray(np.zeros((1,1,h+2,w+2),dtype=org.dtype))
    #---#
    dst[0,0,1:-1,1:-1] = org
    #---#
    dst[0,0,0,1]       =     org[ 0, 0, 0, 0] -     org[ 0, 0, 0, 1]                           + org[ 0, 0, 1, 1]
    dst[0,0,1,0]       =     org[ 0, 0, 0, 0]                        -    org[ 0, 0, 1, 0]     + org[ 0, 0, 1, 1]
    dst[0,0,0,0]       = 2 * org[ 0, 0, 0, 0] - 2 * org[ 0, 0, 0, 1] - 2* org[ 0, 0, 1, 0] + 3 * org[ 0, 0, 1, 1]
    #---#
    dst[0,0,-1,1]      =     org[ 0, 0,-1, 0] -     org[ 0, 0,-1, 1]                           + org[ 0, 0,-2, 1]
    dst[0,0,-2,0]      =     org[ 0, 0,-1, 0]                        -    org[ 0, 0,-2, 0]     + org[ 0, 0,-2, 1]
    dst[0,0,-1,0]      = 2 * org[ 0, 0,-1, 0] - 2 * org[ 0, 0,-1, 1] - 2* org[ 0, 0,-2, 0] + 3 * org[ 0, 0,-2, 1]
    #---#
    dst[0,0,0,-2]      =     org[ 0, 0, 0,-1] -     org[ 0, 0, 0,-2]                           + org[ 0, 0, 1,-2]
    dst[0,0,1,-1]      =     org[ 0, 0, 0,-1]                        -    org[ 0, 0, 1,-1]     + org[ 0, 0, 1,-2]
    dst[0,0,0,-1]      = 2 * org[ 0, 0, 0,-1] - 2 * org[ 0, 0, 0,-2] - 2* org[ 0, 0, 1,-1] + 3 * org[ 0, 0, 1,-2]
    #---#
    dst[0,0,-1,-2]     =     org[ 0, 0,-1,-1] -     org[ 0, 0,-1,-2]                           + org[ 0, 0,-2,-2]
    dst[0,0,-2,-1]     =     org[ 0, 0,-1,-1]                        -    org[ 0, 0,-2,-1]     + org[ 0, 0,-2,-2]
    dst[0,0,-1,-1]     = 2 * org[ 0, 0,-1,-1] - 2 * org[ 0, 0,-1,-2] - 2* org[ 0, 0,-2,-1] + 3 * org[ 0, 0,-2,-2]
    #---#
    for m in range(1,w-1):
      dst[ 0, 0, m+1, 0]  = org[0,0,m,0]  
      dst[ 0, 0, m+1,-1]  = org[0,0,m,-1] 
    for n in range(1,h-1):
      dst[ 0, 0, 0, n+1]  = org[0,0,0,n]
      dst[ 0, 0,-1, n+1]  = org[0,0,-1,n]
    return dst,

def pde_jacobi(x, rhs, pde, boundary, num_iter=15, alpha=0.2):
  for i in range(num_iter):
    lhs = pde.forward(x)
    x_ = x.data[0,0,1:-1,1:-1] + alpha*(rhs - lhs)
    x = boundary( x_, x )
  return x

if __name__ == '__main__':

  import matplotlib.pyplot as plt
  from matplotlib import animation

  cuda.get_device(0).use()

  sigma_l = 8.0e-6
  sigma_t = 1.0e-9
  ds = 2.0e-3
  im_w = 15
  im_h = 15

  fig = plt.figure(figsize=(5,5))
  im = plt.imshow(np.zeros((250,250),dtype=np.float32), vmin = -10., vmax = 10., cmap='gray', interpolation='none')
  plt.axis('off')

  def sim():
    parab = Parabolic(sigma_l, sigma_t, ds)
    img = Variable(
        xp.asarray( np.random.uniform(0, 1.0, (1,1,im_w,im_h)), dtype=np.float32), 
        volatile=True)
    rhs = Variable( 
        parab.forward(img).data + xp.array( np.random.uniform(-0.1, 0.1,(1,1,im_w-2,im_h-2)),dtype=np.float32),
        volatile=True)
    bndry   = IsolatedBoundary()

    #img = pde_jacobi(img, rhs, parab, bndry)
    alpha=0.2
    for i in range(15):
      lhs = parab.forward(img)
      #img -= 1.0*(rhs - 0.25*lhs)
      img_ = img.data[0,0,1:-1,1:-1] - alpha*(rhs - lhs)
      img = bndry(img_)
      if i%1==0:
        #yield img, i
        yield lhs, i

    yield rhs, -1

  def draw(data):
    img, i = data[0], data[1]
    print i
    print img.data
    im.set_array(img.data.get()[0,0,:,:])
    return [im]

  ani = animation.FuncAnimation(fig, draw, sim, blit=False, interval=100, repeat=False)
  plt.show()

