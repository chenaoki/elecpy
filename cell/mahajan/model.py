import os, sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from cellmodel import cellmodel
from .const import const_d

class model(cellmodel):

    def __init__(self, shape):

        super(model, self).__init__()

        self.shape = shape
        self.const = const_d
        self.params = [
            'v',
            'temp',
            'dt',
            'm',
            'h',
            'j',
            'dyad',
            'c1',
            'c2',
            'xi1ca',
            'xi1ba',
            'xi2ca',
            'xi2ba',
            'xr',
            'cai',
            'xs1',
            'xs2',
            'xtos',
            'xtof',
            'ytos',
            'ytof',
            'nai',
            'submem',
            'nsr',
            'jsr',
            'xir',
            'tropi',
            'trops',
            'it',
            'st',
            'xina',
            'xik1',
            'xikr',
            'xiks',
            'xito',
            'xitof',
            'xitos',
            'xiNaCa',
            'xica',
            'xiNaK',
            'c_sodium',
            'c_brugada',
            'sw_it'
        ]
        self.create()
        path = os.path.dirname(__file__)
        self.elemwise_path = os.path.join(path, 'kernel.c')
        self.set_function()
