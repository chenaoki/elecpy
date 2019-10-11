from .Stimulator import Stimulator

class MembraneStimulator(Stimulator):
    
    def __init__(self, shape, start, interval, duration, amplitude, size, train=1, name="noname"):
        
        super(MembraneStimulator,self).__init__(shape, start, interval, duration, amplitude, train)

        self.size = size
        self.name = name
        
        if self.name == "area_L" : self.set_area_L(self.size) 
        if self.name == "area_R" : self.set_area_R(self.size) 
        if self.name == "area_T" : self.set_area_T(self.size) 
        if self.name == "area_B" : self.set_area_B(self.size) 
        if self.name == "edge_L" : self.set_edge_L(self.size) 
        if self.name == "edge_R" : self.set_edge_R(self.size) 
        if self.name == "edge_T" : self.set_edge_T(self.size) 
        if self.name == "edge_B" : self.set_edge_B(self.size) 
        if self.name == "point"  : self.set_point(self.size) 

    def set_area_L(self, border):
        assert border > 0 and border < self.shape[1]
        self._map_on[:,:border] = self.amplitude

    def set_area_R(self, border):
        assert border > 0 and border < self.shape[1]
        self._map_on[:,border:] = self.amplitude

    def set_area_T(self, border):
        assert border > 0 and border < self.shape[0]
        self._map_on[:border,:] = self.amplitude

    def set_area_B(self, border):
        assert border > 0 and border < self.shape[0]
        self._map_on[border:,:] = self.amplitude

    def set_edge_L(self, width = 1):
        assert width > 0 and width < self.shape[1]
        self._map_on[:,:width] = self.amplitude

    def set_edge_R(self, width = 1):
        assert width > 0 and width < self.shape[1]
        self._map_on[:,-width:] = self.amplitude

    def set_edge_T(self, width = 1):
        assert width > 0 and width < self.shape[0]
        self._map_on[:width,:] = self.amplitude

    def set_edge_B(self, width = 1):
        assert width > 0 and width < self.shape[0]
        self._map_on[-width:,:] = self.amplitude

    def set_point( self, y, x, rad ):
        assert x >= rad and x + rad <= self.shape[0]
        assert y >= rad and y + rad <= self.shape[1]
        self._map_on[x-rad:x+rad,y-rad:y+rad] = self.amplitude
    
