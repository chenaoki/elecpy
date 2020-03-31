import os, sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from cellmodel import cellmodel
from const import const_d

class model(cellmodel):

    def __init__(self, shape):

        super(model, self).__init__()

        self.shape = shape
        self.const = const_d
        self.params = [
            'v',
            'it',
           # 't',
           # 'udt',
            'dt',
            'st',
            'nai',
            'ki',
            'cai',
            'ina',
            'inab',
            'iki',
            'inak',
            'inaca',
            'ikr',
            'iks',
            'ito',
            'ikur',
            'ikach',
            'ilcatot',
            'm',
            'h',
            'j',
            'd',
            'f',
            'xs',
            'xr',
            'ato',
            'iito',
            'uakur',
            'uikur',
            'fca',
            'ireljsrol',
            'jsr',
            'nsr',
            'trpn',
            'cmdn',
            'fibro',
           # 'csqn',
           # 'utsc',
            'urel',
            'vrel',
            'wrel',
            'yach',
            'APD_change_h',
            'APD_change_j',
            'APD_change_m',
            'APD_change_gna',
            'APD_change_inak',
            'APD_change_inab',
            'APD_change_ina',
            'APD_change_ikur',
            'APD_change_ilca',
        ]
        self.create()
        path = os.path.dirname(__file__)
        self.elemwise_path = os.path.join(path, 'kernel.c')
        self.set_function()
