import numpy as np
import scipy.sparse as sp

import pyculib.sparse as cusp

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
    def coord_real( ivec ):
        return ( ivec // w, ivec % w )
    def coord_all( ivec ):
        return ( ivec // (w+2), ivec % (w+2) )
    
    # Isolated Boundary condition : A1(size_real, size_all)
    A1 = sp.lil_matrix((size_all, size_real), dtype=np.float64)

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
    # A2 = np.zeros((size_real, size_all), dtype=np.float64)
    A2 = sp.lil_matrix((size_real, size_all), dtype=np.float64)

    for ivec in range(size_real):
        y, x = coord_real(ivec)
        A2[ ivec, ivec_all(y+1, x+1)] = -(2*xx+2*yy)/(ds**2)
        A2[ ivec, ivec_all(y  , x+1)] = yy/(ds**2)
        A2[ ivec, ivec_all(y+2, x+1)] = yy/(ds**2)
        A2[ ivec, ivec_all(y+1, x  )] = xx/(ds**2)
        A2[ ivec, ivec_all(y+1, x+2)] = xx/(ds**2)
    A2=sp.csr_matrix(A2)

    A = A2.dot(A1)
    A_diag = A.diagonal()
    
    Dinv = 1/A_diag
    R = cusp.csr_matrix( A - sp.diags( [ A_diag ], [0] ) )
    A = cusp.csr_matrix(A)

    return A, R, Dinv

class PDE(object):

    def __init__(self, h, w, xx, yy, ds):
        self.h = h
        self.w = w
        self.shape = (h*w,)
        self.A, self.R, self.Dinv = getPDEMatrix(h,w,xx,yy,ds) 
        self.handl = cusp.Sparse()         
        self.descr = self.handl.matdescr() # matrix descriptor

    def forward(self, x):
        assert x.shape == self.shape 
        x_ = np.copy(x)
        self.handl.csrmv( 
            trans='N', # transform, N: do nothing
            m=self.A.shape[0], n=self.A.shape[1], 
            nnz = self.A.nnz, 
            alpha=1.0, 
            descr=self.descr, 
            csrVal=self.A.data, 
            csrRowPtr=self.A.indptr, 
            csrColInd=self.A.indices,  
            x=x_, 
            beta=0.0, 
            y=x_)
        return x_

    def solve(self, x, y, tol=1e-4, cdcnt = 3, maxcnt = 300):
        assert x.shape == self.shape
        assert y.shape == self.shape
        
        x_    = np.copy(x)
        x_old = np.copy(x_)
        norm_y = np.linalg.norm(y)
        
        error = 1e12
        cnt = 0
        while error > tol:
            
            # x_ = ( y - self.R.dot(x_old) ).multiply( self.Dinv )
            self.handl.csrmv( 
                trans='N', # transform, N: do nothing
                m=self.R.shape[0], n=self.R.shape[1], 
                nnz = self.R.nnz, 
                alpha=1.0, 
                descr=self.descr, 
                csrVal=self.R.data, 
                csrRowPtr=self.R.indptr, 
                csrColInd=self.R.indices,  
                x=x_old, 
                beta=0.0, 
                y=x_)
            x_ = ( y - x_ )*( self.Dinv )
            
            # convergence determination
            if cnt % cdcnt == 0:
                error = np.linalg.norm(y - self.forward(x_)) / norm_y
            
            x_old = x_
            cnt += 1
            if cnt >= maxcnt : break
            
        return cnt, x_

if __name__ == '__main__':

    import matplotlib.pyplot as plt
    from matplotlib import animation

    fig = plt.figure(figsize=(5,5))
    im = plt.imshow(np.zeros((250,250),dtype=np.float64), vmin = -10., vmax = 10., interpolation='nearest')
    plt.axis('off')

    def sim():
        w = 30
        h = 30

        x = np.ones((h,w), dtype=np.float64)
        y = np.zeros((h,w), dtype=np.float64)
        y[h/4:3*h/4, w/4:3*w/4] = 1.0

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
