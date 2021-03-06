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
            'dt',
            'temp',
            'st',
            'k_i',
            'na_i',
            'r',
            's',
            'ikv',
            'ik1',
            'inak',
            'inab',
            'ik_net',
            'ina_net',
            'sw_it',
        ]
        self.create()
        path = os.path.dirname(__file__)
        self.elemwise_path = os.path.join(path, 'kernel.c')
        self.set_function()
