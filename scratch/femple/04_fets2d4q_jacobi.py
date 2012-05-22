'''
Created on Apr 23, 2012

@author: rch
'''

import numpy as np
from enthought.traits.api import \
    HasTraits, Array, Tuple, Constant, Property, cached_property, Float
import sympy as sp

x_, y_, z_, t_ = sp.symbols('x,y,z,t')

class FETS2D4Q(HasTraits):
    '''Quadratic finite element.
    '''

    #===========================================================================
    # Specification of nodal points
    #===========================================================================
    geo_r = Array(Float, value = [[-1, -1], [1, -1], [1, 1], [-1, 1]])

    #===========================================================================
    # Kinematic mapping
    #===========================================================================
    dN_idx_map = Constant(value = ([0, 1, 1, 0],
                                    slice(None)))
    
    B_idx_map = Constant(value = ([0, 1, 2, 2],
                                  slice(None),
                                  [0, 1, 0, 1]))

    #===========================================================================
    # Generation of the shape function operator
    #===========================================================================
    dN_fn = Property
    def _get_dN_fn(self):
        
        P = [1, x_, y_, x_ * y_]
        PX = sp.lambdify([x_, y_], P)
        
        C = np.array([ PX(xi[0], xi[1]) for xi in self.geo_r ], dtype = 'f')
        
        C1 = np.linalg.inv(C)
        
        P_arr = np.array(P).reshape(1, 4)
        
        N_mtx = sp.Matrix(np.dot(P_arr, C1))
        N_mtx.simplify()
        
        dN_mtx = N_mtx.jacobian([x_, y_]).T
        dN_mtx.simplify()
        
        return sp.lambdify([x_, y_], dN_mtx)        
        
    def get_B_mtx(self, r, s):
        '''Return kinematic matrix at the given points.
        '''
        B_mtx = np.zeros ((3, 4, 2), dtype = 'f')
        dN_mtx = self.dN_fn(r, s)
        B_mtx[self.B_idx_map] = dN_mtx[self.dN_idx_map]
        return B_mtx.reshape(3, 8)

if __name__ == '__main__':

    fets = FETS2D4Q()
    
    gauss_points = np.array(
                            [[-0.5777, -0.5777],
                             [0.5777, -0.5777],
                             [0.5777, 0.5777],
                             [-0.5777, 0.5777],
                    ], dtype = 'f')
    print 'gp', gauss_points[:, 0]
    print 'dN_fn', fets.dN_fn(0., 0.)
    print "B_mtx_map = \n", fets.get_B_mtx(gauss_points[:, 0], gauss_points[:, 1])

    fets.configure_traits()
