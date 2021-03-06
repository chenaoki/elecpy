# -*- coding:utf-8 -*-

from Stimulator import Stimulator
import numpy as np

class ExtracellularStimulator(Stimulator):

    def __init__(self, shape, start, interval, duration, amplitude, size, train=1, name="noname"):
        
        super(ExtracellularStimulator,self).__init__(shape, start, interval, duration, amplitude, train)
        
        self.size = size
        self.name = name
        self.is_anode     = True if amplitude > 0 else False
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

    def set_location(self):
        anode   = (( self._map_anode > 0 )*1).astype(np.float32)
        cathode = (( self._map_cathode > 0 )*1).astype(np.float32)
        self._map_on *= 0.0 # reset
        if self.is_anode:
            self._map_on += anode*self.amplitude
            self._map_on += cathode*( -self.amplitude * np.sum(anode)/ float(np.sum(cathode)) ) 
        else:
            self._map_on += cathode*-self.amplitude
            self._map_on += anode*( self.amplitude * np.sum(cathode)/ float(np.sum(anode)) ) 

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
            self._map_cathode[:, 0] = 1.0
            self._map_cathode[:, -1] = 1.0
        else:
            self._map_cathode[x-rad:x+rad,y-rad:y+rad] = 1.0
            self._map_anode[:, 0] = 1.0
            self._map_anode[:, -1] = 1.0
        self.set_location()