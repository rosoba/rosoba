'''
Created on Apr 23, 2012

@author: rch
'''

import numpy as np
from enthought.traits.api import HasTraits, Array, Tuple, Constant
import sympy as sp

class FETS2D4Q(HasTraits):
    '''Quadratic finite element.
    '''
    
    dN_idx_map = Constant(value = ([0, 1, 1, 0],
                                    slice(None)))
    
    B_idx_map = Constant(value = ([0, 1, 2, 2],
                                  slice(None),
                                  [0, 1, 0, 1]))

    def get_dN_mtx(self, points):
        '''Return the matrix of global shape functions.
        '''
        return np.array([[11, 21, 31, 41],
                         [12, 22, 32, 42]],
                        dtype = 'f')
        
    def get_B_mtx(self, points):
        '''Return kinematic matrix at the given points.
        '''
        B_mtx = np.zeros ((3, 4, 2), dtype = 'f')
        dN_mtx = self.get_dN_mtx(points)
        B_mtx[self.B_idx_map] = dN_mtx[self.dN_idx_map]
        return B_mtx.reshape(3, 8)

if __name__ == '__main__':

    fets = FETS2D4Q()
    
    print "B_mtx_map = \n", fets.get_B_mtx(None)

    fets.configure_traits()
