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
            'temp',
            'dt',
            'h',
            'j',
            'm',
            'nai',
            'naiss',
            'cai',
            'caiss',
            'ki',
            'kiss',
            'ltypeCzero',
            'ltypeCone',
            'ltypeCtwo',
            'ltypeCthree',
            'ltypeIVf',
            'ltypeIVs',
            'ltypeO',
            'fcasc',
            'fmode0',
            'ilca',
            'ilcana',
            'ilcak',
            'irel',
            'itr',
            'b',
            'g',
            'xr',
            'xs1',
            'xs2',
            'inacass',
            'nsr',
            'ryrCone',
            'ryrCtwo',
            'ryrCthree',
            'ryrCfour',
            'ryrOone',
            'ryrIone',
            'ryrItwo',
            'ryrIthree',
            'ryrIfour',
            'ryrIfive',
            'csqn',
            'jsr',
            'it',
            'st'
        ]
        self.create()
        path = os.path.dirname(__file__)
        self.elemwise_path = os.path.join(path, 'kernel.c')
        self.set_function()
