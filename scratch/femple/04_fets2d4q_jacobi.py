'''
Created on Apr 23, 2012

@author: rch
'''

import numpy as np 
from numpy import \
     array, zeros, int_, float_, ix_, dot, linspace, hstack, vstack, arange, \
     identity
from traits.api import \
    HasTraits, Array, Tuple, Constant, Property, cached_property, Float
import sympy as sp
from scipy.linalg import \
     inv, det

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

    #===========================================================================
    # Specification of global points
    #===========================================================================
    geo_x = Array(Float, value = [[0, 0], [3, 1], [2, 4], [-1, 2]])
    #===========================================================================
    # Generation of the B_Matrix
    #===========================================================================
    b_mtx = Property
    
    def _get_b_mtx(self):
        gaus= np.array([[-5.777,-5.777],[5.777,-5.777],[5.777,5.777],[-5.777,5.777]], dtype ='f')
        b_mtx = zeros ((2,4), dtype = 'f')

        for i in range (0,3):
            dN = self.dN_fn(gaus[i,0],gaus[i,1])
            J_mtx = dot(dN,self.geo_x)
            J_mtx_inv = inv(J_mtx)
            b_mtx = b_mtx + dot(J_mtx_inv,dN)
            
        B = zeros ((3, 4, 2), dtype = 'f')


        B[self.B_idx_map] = b_mtx[self.dN_idx_map]

        return B.reshape(3, 8)
    
    
    
    #===========================================================================
    # Kinematic mapping
    #===========================================================================


if __name__ == '__main__':

    fets = FETS2D4Q()
    
   
    
    
    

print "B_mtx_map = \n", fets.b_mtx

fets.configure_traits()
