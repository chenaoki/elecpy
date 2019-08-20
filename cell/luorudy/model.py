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
            'nai',
            'dt',
            'm',
            'h',
            'j',
            'ltypeCzero',
            'ltypeCone',
            'ltypeCtwo',
            'ltypeCthree',
            'ltypeIVf',
            'ltypeIVs',
            'ltypeO',
            'caiss',
            'ilca',
            'inacass',
            'irel',
            'naiss',
            'ilcana',
            'kiss',
            'ki',
            'ilcak',
            'cai',
            'fcasc',
            'fmode0',
            'b',
            'g',
            'xr',
            'xs1',
            'xs2',
            'st',
            'itr',
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
            'sponrel',
            'jsr',
            'it', 
            'ina',
            'c_brugada',
            'sw_it'
        ]
        self.create()
        path = os.path.dirname(__file__)
        self.elemwise_path = os.path.join(path, 'kernel.c')
        self.set_function()
