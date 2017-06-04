import os
import numpy as np
import time
from optparse import OptionParser

import chainer
import chainer.functions as F
from chainer import cuda
from chainer import Function
xp = cuda.cupy

class cellmodel(object):

    def __init__(self):
        self.shape = None
        self.const = None
        self.params = None
        self.func = None
        self.state = None
        self.elemwise_path = None

    def set_function(self):
        assert self.params is not None
        assert self.const is not None
        assert self.elemwise_path is not None

        string_args = ''
        string_rets = ''
        for i, key in enumerate(self.params) :
          string_args += 'T {0}, '.format(key)
          string_rets += 'T _{0}, '.format(key)
        string_args = string_args.rstrip(', ')
        string_rets = string_rets.rstrip(', ')
        elemwise = open(self.elemwise_path).read().format(**self.const)

        class f_model(Function):
          def forward(self,x):
            state = x
            _state = cuda.elementwise(
              string_args,
              string_rets,
              elemwise,
              'luorudy')(*state)
            return _state

        self.func = f_model()

    def create(self):
        assert self.shape is not None
        assert self.params is not None
        self.state = [ xp.asarray( np.ones(self.shape, dtype=np.float32) * self.const[param+'_'] ) for param in self.params ]
        pass

    def save(self, path):
        assert self.params is not None
        if not os.path.isdir(path) : os.mkdir(path)
        for i, param in enumerate(self.params):
            np.save(path+'/'+param, self.state[i])
        pass

    def load(self, path):
        assert self.params is not None
        self.state = [ xp.asarray( np.load(path+'/'+param+'.npy')) for param in self.params ]
        pass

    def set_param(self, param_name, param_value):
        assert param_name in self.params
        print self.shape, param_value.shape
        assert param_value.shape == self.shape
        self.state[self.params.index(param_name)] = param_value

    def get_param(self, param_name):
        assert param_name in self.params
        return cell_state[model_params.index(param_name)]

    def update(self):
        assert self.func is not None
        self.state = list( self.func.forward(tuple(self.state)) )
        pass
