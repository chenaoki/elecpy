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
            'CaMKt',
            'nai',
            'nass',
            'ki',
            'kss',
            'cai',
            'cass',
            'm',
            'hf',
            'hs',
            'j',
            'jp',
            'hsp',
            'mL',
            'hL',
            'hLp',
            'a',
            'iF',
            'iS',
            'ap',
            'iFp',
            'iSp',
            'd',
            'ff',
            'fs',
            'fcaf',
            'fcas',
            'jca',
            'ffp',
            'fcafp',
            'nca',
            'xrf',
            'xrs',
            'xs1',
            'xs2',
            'xk1',
            'Jrelnp',
            'Jrelp',
            'cansr',
            'cajsr',
            'it',
            'st'
        ]
        self.create()
        path = os.path.dirname(__file__)
        self.elemwise_path = os.path.join(path, 'kernel.c')
        self.set_function()
