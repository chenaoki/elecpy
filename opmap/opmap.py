import numpy as np
from glob import glob
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from cmap_bipolar import bipolar
import scipy
from scipy import signal, ndimage
from scipy.interpolate import interp1d, splrep, splev
import pickle
import cv2
from pylab import tight_layout 

cam_dtype={
    'sa4':np.ushort,
    'mini':np.ushort,
    'max':np.uint8,
    'max10':np.ushort
}

class VideoData(object):

    def __init__(self, length, height, width ):
        self.data = np.zeros((length, height, width), dtype=np.float32)
        self.roi = np.ones( ( height, width), dtype=np.float32 )
        self.vmin =  0.0
        self.vmax = 1.0
        self.cmap = 'hot'

    def showFrame(self, frame):
        assert frame >= 0 and frame < self.data.shape[0]
        plt.imshow(self.data[frame, :, :], vmin=self.vmin, vmax=self.vmax, cmap=self.cmap)


    def saveMovie(self, path, fps=30, dpi=100, interval=1):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_aspect('equal')
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        im = ax.imshow(self.data[0,:,:],vmin=self.vmin, vmax=self.vmax, cmap=self.cmap)
        fig.set_size_inches([5,5])
        tight_layout()
        def update_img(n):
          im.set_data(self.data[ (interval * n ) % self.data.shape[0],:,:])
          return im
        ani = animation.FuncAnimation(fig, update_img, self.data.shape[0] / interval)
        writer = animation.writers['ffmpeg'](fps=fps)
        ani.save(path, writer=writer, dpi=dpi)

    def setROI(self, top=None, bottom=None, left=None, right=None):
        if top is not None :
            assert top >= 0 and top < self.roi.shape[0]
            self.roi[:top, :] = 0
        if bottom is not None :
            assert bottom >= 0 and bottom < self.roi.shape[0]
            self.roi[bottom:, :] = 0
        if left is not None :
            assert left >= 0 and left < self.roi.shape[1]
            self.roi[:, :top] = 0
        if right is not None :
            assert right >= 0 and right < self.roi.shape[1]
            self.roi[:, right:] = 0
        self.data *= self.roi

    def showROI(self):
        plt.imshow(self.roi, vmin=0.0, vmax=1.0, cmap='gray')

    def saveImage(self, savedir, img_type = 'png'):
        if not os.path.exists(savedir):
            os.makedirs(savedir)
        for frame in range(self.data.shape[0]):
            plt.imsave(
                '{0}/{1:0>6}.{2}'.format(savedir, frame, img_type),
                self.data[frame, :, :], vmin=self.vmin, vmax=self.vmax, cmap=self.cmap
            )
        plt.imsave(
         '{0}/roi.{1}'.format(savedir, img_type),
            self.roi, vmin=0.0, vmax=1.0, cmap='gray'
        )

    def plot(self, points, start=None, end=None, filter_size=None, savepath = None):
        if start is None : start = 0
        if end is None : end = self.data.shape[0]
        fig = plt.figure()
        for i, p in enumerate(points):
            assert len(p) == 2
            x, y = p
            assert y >= 0 and y < self.data.shape[1]
            assert x >= 0 and x < self.data.shape[2]
            ax = fig.add_subplot(len(points),1,i+1)
            ts = self.data[start:end, y, x]
            if filter_size is not None:
                ts = signal.savgol_filter(ts, filter_size, 3)
            ax.plot(ts)
        if savepath is None:
            plt.show()
        else:
            plt.savefig(savepath)

class RawCam( VideoData ):

	def __init__(self, path, cam_type, image_width, image_height, frame_start, frame_end):

		if "numpy" == cam_type:

			self.files = sorted(glob(path+"/vmem_*.npy"))
			assert len(self.files) > 0
			if frame_end < 0 : frame_end = len(self.files) + frame_end + 1
			self.files = self.files[frame_start:frame_end]

			super(RawCam, self).__init__(len(self.files), image_height, image_width)

			for i, f in enumerate(self.files):
				im = np.load(f)
				self.data[i, :,:] = im

		else:

			self.files = sorted(glob(path+"/*.raw*"))
			assert len(self.files) > 0
			self.files = self.files[frame_start:frame_end]

			super(RawCam, self).__init__(len(self.files), image_height, image_width)

			for i, f in enumerate(self.files):
				im = np.fromfile(f, dtype=cam_dtype[cam_type])
				im = im.reshape(image_height, image_width)
				self.data[i, :,:] = im

		self.vmin = np.min(self.data)
		self.vmax = np.max(self.data)
		self.cmap = 'gray'

		return

	def selectPoints(self, savepath="./selectPoints.png"):

		points = []
		def onClick(event, x, y, flag, params):
			wname, img = params
			if event == cv2.EVENT_LBUTTONDOWN:
				img_disp = np.copy(img)
				points.append((x, y))
				for i, p in enumerate(points):
					cv2.circle(img_disp, p, 2, (255,0,0))
					cv2.putText(img_disp,str(i),(p[0]-5, p[1]-5),cv2.FONT_HERSHEY_PLAIN, 0.6,(255,0,0))
					cv2.imshow(wname, img_disp)
					cv2.imwrite(savepath, img_disp)

		wname = "selectPoints"
		img = self.data[0,:,:]

		assert self.vmin < self.vmax
		img = (255 * (img - self.vmin) / float(self.vmax - self.vmin))
		img = img.astype(np.uint8)
		img = cv2.cvtColor(img, cv2.COLOR_GRAY2RGB)

		cv2.namedWindow(wname)
		cv2.setMouseCallback(wname, onClick, [wname, img] )
		cv2.imshow(wname, img)
		while cv2.waitKey(0) != 27 : pass
		cv2.destroyWindow(wname)
		return np.array(points)


class VmemMap( VideoData ):

    def __init__(self, rawcam):

    	shape = rawcam.data.shape
    	super(VmemMap, self).__init__(shape[0], shape[1], shape[2])

        im_max = np.max(rawcam.data, axis=0)
        im_min = np.min(rawcam.data, axis=0)
        self.im_range = (im_max - im_min) + (im_max == im_min) * 1
        self.data_org = 2.0 * (im_max - rawcam.data ) / self.im_range - 1.0
        self.data = np.copy(self.data_org)
        self.diff_max = rawcam.vmax - rawcam.vmin

        self.vmin = -1.0
        self.vmax = 1.0            
        self.cmap = bipolar(neutral=0, lutsize=1024)
        return

    def setDiffRange(self, diff_min=None, diff_max=None):
        self.roi = 1.0 # reset
        if diff_min is None :
          diff_min = 0
        if diff_max is None :
          diff_max = self.diff_max
        self.roi *= (self.im_range>=diff_min)*1
        self.roi *= (self.im_range<=diff_max)*1
        self.data = self.data_org*self.roi
        return

    def smooth(self, size):
        assert size > 0
        for frame in range( self.data.shape[0]):
            self.data[frame,:,:] = ndimage.gaussian_filter(self.data[frame,:,:], sigma = size)
        return

class PhaseMap( VideoData ):

    def __init__(self, vmem, shrink = 4, fs = 1000.0, cutoff = 20.0):

        self.shrink = shrink
        size_org = vmem.data.shape

        super(PhaseMap, self).__init__(size_org[0],size_org[1]/shrink, size_org[2]/shrink)

        nyq = fs/2.0
        fe = cutoff / nyq   # Cut off frequency : 20Hz
        numtaps = 15  # Filter size
        b = signal.firwin(numtaps, fe) # Low pass filter

        for n in range(self.data.shape[1]):
            for m in range(self.data.shape[2]):
                n_ = n*shrink
                m_ = m*shrink
                try:
                    assert vmem.roi[n_, m_] > 0
                    data = vmem.data[ :, n_, m_]
                    data_an = signal.lfilter(b, 1, data)
                    data_max = np.max(data_an)
                    data_min = np.min(data_an)
                    assert data_max > data_min
                    data_an  = 2.0 * (data_an - data_min) / (data_max - data_min) - 1.0
                    self.data[:, n, m] = np.angle(signal.hilbert(data_an))
                    self.roi[n, m] = 1.0
                except:
                    self.data[:, n, m] = 0.0
                    self.roi[n, m] = 0.0

        self.vmin = -np.pi
        self.vmax = np.pi
        self.cmap = 'jet'
        return

    def smooth(self, size = 4):
        assert size > 0
        def phaseComplement(value):
            value -= (value > np.pi)*2*np.pi
            value += (value < - np.pi)*2*np.pi
            return value
        new_data = np.zeros_like(self.data)
        for frame in range( self.data.shape[0]):
            if frame % 10 == 0 : print frame
            for n in range(self.data.shape[1])[size:-size]:
                for m in range(self.data.shape[2])[size:-size]:
                    base = self.data[frame, n, m]
                    target = self.data[frame, n-size:n+size+1, m-size:m+size+1]
                    difference = phaseComplement(target-base)
                    diff = signal.medfilt2d(difference, kernel_size=size*2+1)[size,size]
                    new_data[frame, n, m] = phaseComplement( base +  diff)
                    #new_data[frame, n, m] = phaseComplement( base + np.mean(difference.flatten()) )
        self.data = new_data

class PhaseVarianceMap( VideoData ):

    def __init__(self, phasemap, size = 9):
        assert size > 0

        super(PhaseVarianceMap, self).__init__(*phasemap.data.shape)
        kernel = np.ones((size, size), dtype=np.float32)
        kernel /= np.sum(kernel)

        for frame in range(self.data.shape[0]):
            im_cos = np.cos(phasemap.data[frame,:,:])
            im_sin = np.sin(phasemap.data[frame,:,:])
            im_cos = signal.convolve2d(im_cos, kernel, mode = 'same', boundary = 'fill')
            im_sin = signal.convolve2d(im_sin, kernel, mode = 'same', boundary = 'fill')
            self.data[frame, :, :] = 1.0 - np.abs( im_cos + 1j * im_sin )

        self.roi = scipy.ndimage.binary_closing(phasemap.roi, structure=np.ones((size/2,size/2))).astype(phasemap.roi.dtype)
        self.roi = scipy.ndimage.binary_erosion(self.roi, structure=np.ones((size,size))).astype(phasemap.roi.dtype)
        self.data *= self.roi

        self.vmin = 0.0
        self.vmax = 1.0
        self.cmap = 'gray'

class CoreMap( VideoData ):

  def __init__(self, pvmap, threshold=0.8):
    assert threshold >= 0.0 and threshold <= 1.0
    self.threshold = threshold

    super(CoreMap, self).__init__(*pvmap.data.shape)

    self.roi = pvmap.roi
    self.data, self.coreNum = ndimage.label((pvmap.data>self.threshold)*1) 

    self.vmin = 0
    self.vmax = np.max(self.data) 
    self.cmap = plt.cm.spectral 

  def getFrameCore(self, frame):
    im = (self.data[frame,:,:]>0)*1
    im_label, n = ndimage.label(im) 
    return im_label, ndimage.measurements.center_of_mass(im, im_label, range(1,n+1))
