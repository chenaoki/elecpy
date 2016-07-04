import numpy as np
import scipy.sparse as sp
from scipy.linalg import norm as spnorm

def getPDEMatrix(h, w, xx, yy, ds):

  size_real = h*w
  size_all  = (h+2)*(w+2)

  def ivec_real( y, x ):
    _x = x if x >=0 else w+x 
    _y = y if y >=0 else h+y 
    return w*_y + _x
  def ivec_all( y, x ):
    _x = x if x >=0 else (w+2)+x 
    _y = y if y >=0 else (h+2)+y 
    return (w+2)*_y + _x


  # Isolated Boundary condition : A1(size_real, size_all)
  A1 = sp.lil_matrix((size_all, size_real), dtype=np.float32)

  # - Real points
  for n in range(h):
    for m in range(w):
      A1[ ivec_all(n+1, m+1), ivec_real(n,m)] = 1.0

  # - 4 corners
  A1[ ivec_all( 0, 0), ivec_real( 0, 0)] =  2.0 
  A1[ ivec_all( 0, 0), ivec_real( 1, 0)] = -2.0 
  A1[ ivec_all( 0, 0), ivec_real( 0, 1)] = -2.0
  A1[ ivec_all( 0, 0), ivec_real( 1, 1)] =  3.0
  A1[ ivec_all( 1, 0), ivec_real( 0, 0)] =  1.0
  A1[ ivec_all( 1, 0), ivec_real( 1, 0)] = -1.0
  A1[ ivec_all( 1, 0), ivec_real( 1, 1)] =  1.0
  A1[ ivec_all( 0, 1), ivec_real( 0, 0)] =  1.0
  A1[ ivec_all( 0, 1), ivec_real( 0, 1)] = -1.0
  A1[ ivec_all( 0, 1), ivec_real( 1, 1)] =  1.0

  A1[ ivec_all( -1, 0), ivec_real( -1, 0)] =  2.0 
  A1[ ivec_all( -1, 0), ivec_real( -2, 0)] = -2.0 
  A1[ ivec_all( -1, 0), ivec_real( -1, 1)] = -2.0
  A1[ ivec_all( -1, 0), ivec_real( -2, 1)] =  3.0
  A1[ ivec_all( -2, 0), ivec_real( -1, 0)] =  1.0
  A1[ ivec_all( -2, 0), ivec_real( -2, 0)] = -1.0
  A1[ ivec_all( -2, 0), ivec_real( -2, 1)] =  1.0
  A1[ ivec_all( -1, 1), ivec_real( -1, 0)] =  1.0
  A1[ ivec_all( -1, 1), ivec_real( -1, 1)] = -1.0
  A1[ ivec_all( -1, 1), ivec_real( -2, 1)] =  1.0

  A1[ ivec_all( 0, -1), ivec_real( 0, -1)] =  2.0 
  A1[ ivec_all( 0, -1), ivec_real( 1, -1)] = -2.0 
  A1[ ivec_all( 0, -1), ivec_real( 0, -2)] = -2.0
  A1[ ivec_all( 0, -1), ivec_real( 1, -2)] =  3.0
  A1[ ivec_all( 1, -1), ivec_real( 0, -1)] =  1.0
  A1[ ivec_all( 1, -1), ivec_real( 1, -1)] = -1.0
  A1[ ivec_all( 1, -1), ivec_real( 1, -2)] =  1.0
  A1[ ivec_all( 0, -2), ivec_real( 0, -1)] =  1.0
  A1[ ivec_all( 0, -2), ivec_real( 0, -2)] = -1.0
  A1[ ivec_all( 0, -2), ivec_real( 1, -2)] =  1.0

  A1[ ivec_all( -1, -1), ivec_real( -1, -1)] =  2.0 
  A1[ ivec_all( -1, -1), ivec_real( -2, -1)] = -2.0 
  A1[ ivec_all( -1, -1), ivec_real( -1, -2)] = -2.0
  A1[ ivec_all( -1, -1), ivec_real( -2, -2)] =  3.0
  A1[ ivec_all( -2, -1), ivec_real( -1, -1)] =  1.0
  A1[ ivec_all( -2, -1), ivec_real( -2, -1)] = -1.0
  A1[ ivec_all( -2, -1), ivec_real( -2, -2)] =  1.0
  A1[ ivec_all( -1, -2), ivec_real( -1, -1)] =  1.0
  A1[ ivec_all( -1, -2), ivec_real( -1, -2)] = -1.0
  A1[ ivec_all( -1, -2), ivec_real( -2, -2)] =  1.0

  # - 4 edges
  for m in range(1,w-1):
    A1[ ivec_all(  0, m+1), ivec_real(  0, m)] =  1.0
    A1[ ivec_all( -1, m+1), ivec_real( -1, m)] =  1.0

  for n in range(1,h-1):
    A1[ ivec_all( n+1,  0), ivec_real( n,  0)] =  1.0
    A1[ ivec_all( n+1, -1), ivec_real( n, -1)] =  1.0
  
  A1=sp.csr_matrix(A1)

  # Parabolic PDE coefficients
  # A2 = np.zeros((size_real, size_all), dtype=np.float32)
  A2 = sp.lil_matrix((size_real, size_all), dtype=np.float32)

  for y in range(h):
    for x in range(w):
      A2[ ivec_real(y,x), ivec_all(y+1, x+1)] = -(2*xx+2*yy)/(ds**2)
      A2[ ivec_real(y,x), ivec_all(y  , x+1)] = yy/(ds**2)
      A2[ ivec_real(y,x), ivec_all(y+2, x+1)] = yy/(ds**2)
      A2[ ivec_real(y,x), ivec_all(y+1, x  )] = xx/(ds**2)
      A2[ ivec_real(y,x), ivec_all(y+1, x+2)] = xx/(ds**2)
  A2=sp.csr_matrix(A2)

  A = A2.dot(A1)
  Dinv = sp.csr_matrix((1/np.diag(A.todense()))[:, np.newaxis])
  R = A - sp.csr_matrix(np.diag(np.diag(A.todense())))

  #print 'A1', A1
  #print 'A2', A2
  #return cuda.to_gpu(A), cuda.to_gpu(R), cuda.to_gpu(D)
  return A, R, Dinv

class PDE(object):
  def __init__(self, h, w, xx, yy, ds):
    self.h = h
    self.w = w
    self.A, self.R, self.Dinv = getPDEMatrix(h,w,xx,yy,ds) 

  def forward(self, x):
    assert x.shape == (self.h, self.w) 
    x_ = sp.csr_matrix(x.flatten()[:, np.newaxis])
    return np.array(self.A.dot(x_).todense()).reshape((self.h,self.w))

  def solve(self, x, y, tol=1e-4, maxcnt = 300):
    assert x.shape == (self.h, self.w) 
    assert y.shape == (self.h, self.w)
    x_old = sp.csr_matrix(x.flatten()[:, np.newaxis])
    x_    = sp.csr_matrix(x.flatten()[:, np.newaxis])
    y_    = sp.csr_matrix(y.flatten()[:, np.newaxis])
    error = 1e12
    cnt = 0
    while error > tol:
      x_ = ( y_ - self.R.dot(x_old) ).multiply( self.Dinv )
      error = np.linalg.norm((x_ - x_old).todense()) / np.linalg.norm((x_).todense())
      x_old = x_
      cnt += 1
      if cnt >= maxcnt : break
    return cnt, np.array(x_.todense()).reshape((self.h,self.w))

if __name__ == '__main__':

  import matplotlib.pyplot as plt
  from matplotlib import animation

  fig = plt.figure(figsize=(5,5))
  im = plt.imshow(np.zeros((250,250),dtype=np.float32), vmin = -10., vmax = 10., interpolation='nearest')
  plt.axis('off')

  def sim():
    w = 30
    h = 30

    x = np.ones((h,w), dtype=np.float32)
    y = np.zeros((h,w), dtype=np.float32)
    y[h/4:3*h/4, w/4:3*w/4] = 1.0

    ''' 
    A, R, Dinv = getPDEMatrix(h,w,3.75e-3, 3.75e-3, 1.5e-2) 

    #np.save( 'A', A.get() )
    #np.save( 'R', R.get() )
    #np.save( 'D', D.get() )
    #np.save( 'x', x.get() )
    #np.save( 'y', y.get() )

    x_ = x.flatten()[:,np.newaxis]
    y_ = y.flatten()[:,np.newaxis]

    for i in range(1000):
      print 'loop@{0}'.format(i)
      x_ = ( y_ - xp.dot(R, x_ )) / D[:, np.newaxis]
      yield i, x_.get().reshape((h,w))
    
    ''' 

    pde_solver = PDE( h, w, 3.75e-3, 3.75e-3, 1.5e-2) 
    cnt, x = pde_solver.solve(x, y, tol=1e-3)
    yield cnt, x


  def draw(data):
    i, img = tuple(data)
    print i
    im.set_array(img)
    return [im]

  ani = animation.FuncAnimation(fig, draw, sim, blit=False, interval=100, repeat=False)
  plt.show()
