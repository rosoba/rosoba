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
    
    #===========================================================================
    # Specification of global points
    #===========================================================================
    geo_x = Array(Float, value = [[0, 0], [3, 1], [2, 4], [-1, 2]])
    #===========================================================================
    
    def get_B_mtx(self, geo_x,r, s):
        '''Return kinematic matrix at the given points.
        '''
        B_mtx = np.zeros ((3, 4, 2), dtype = 'f')
        
        dN_mtx = self.dN_fn(r, s)
        
        J_mtx = dot(dN_mtx, geo_x)
        
        J_mtx_inv = inv(J_mtx)
        
        dNx_mtx = dot(J_mtx_inv,dN_mtx)
        
        B_mtx[self.B_idx_map] = dNx_mtx[self.dN_idx_map]
        
        return B_mtx.reshape(3, 8)

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
            b_mtx += dot(J_mtx_inv,dN)
            
        B = zeros ((3, 4, 2), dtype = 'f')


        B[self.B_idx_map] = b_mtx[self.dN_idx_map]
        
        return B.reshape(3, 8)
    
    
    
    #===========================================================================
    # Kinematic mapping
    #===========================================================================


if __name__ == '__main__':

    fets = FETS2D4Q ()
 
    
   
    from ibvpy.api import FEGrid
    from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
    
    fets_sample = FETS2D4Q()

    fe_grid = FEGrid( coord_max = ( 2., 3., ),
                              shape = ( 2, 1 ),
                              fets_eval = fets_sample )
    
    X_map = fe_grid.elem_X_map
    
    print X_map [0,:]
    print X_map
    
    gauss= np.array([[-5.777,-5.777],[5.777,-5.777],[5.777,5.777],[-5.777,5.777]], dtype ='f')
    E= 100000.00
    v= 0.2
    D_mtx = (E/(1-v**2))* np.array ([[1,v,0],[v,1,0],[0,0,((1-v)/2)]], dtype = 'f')
    
    K_mtx = zeros ((8,8), dtype ='f')
    K_mtx_2 = zeros ((8,8), dtype ='f')
    
    for i in range (0,2):
    
        X_map = fe_grid.elem_X_map[i,:]
        
        for j in range (0,4):
        
            B_mtx = fets.get_B_mtx (X_map,gauss[j,0],gauss[j,1])
        
            K_mtx += dot(dot(B_mtx.T,D_mtx),B_mtx)
        
        
        print "K_mtx_", i+1, " = \n", K_mtx 
    
    
    #fets = FETS2D4Q()
    
    print fe_grid.elem_dof_map
    print B_mtx
    
    

    #print "B_mtx_map = \n", fets.b_mtx

    #fets.configure_traits()

