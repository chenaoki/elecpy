import numpy as np

class Stimulator(object):
    
    def __init__(self, shape, start, interval, duration, amplitude, train):
        self.shape        = shape
        self.interval     = interval
        self.duration     = duration
        self.start        = start
        self.train        = train
        self.amplitude    = abs(amplitude)  # (uA/cm^2)
        self._map_on      = np.zeros( self.shape, dtype=np.float64 )
        self._map_off     = np.zeros( self.shape, dtype=np.float64 )
    
    def get_flag(self, t):
        if t >= self.start and t < self.start + self.train * self.interval and (t - self.start) % self.interval < self.duration:
            return True
        else:
            return False

    def get_current(self, t):
        assert self._map_on is not None
        assert self._map_off is not None
        return self._map_on if self.get_flag(t) else self._map_off
