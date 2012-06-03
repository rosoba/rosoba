
import sys
print sys.path

from etsproxy.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property

import etsproxy.traits.has_traits
etsproxy.traits.has_traits.CHECK_INTERFACES = 2

from etsproxy.traits.ui.api import \
     Item, View, HGroup, ListEditor, VGroup, Group

from etsproxy.traits.ui.menu import \
     NoButtons, OKButton, CancelButton, Action, CloseAction, Menu, \
     MenuBar, Separator

import numpy as np
import math

from scipy.linalg import \
     inv

from ibvpy.fets.fets_eval import FETSEval

#-----------------------------------------------------------------------------------
# FETS2D4Q - 4 nodes iso-parametric quadrilateral element (2D, linear, Lagrange family)    
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
# Element Information: 
#-----------------------------------------------------------------------------------
#
# Here an isoparametric element formulation is applied.
# The implemented shape functions are derived based on the 
# ordering of the nodes of the parent element defined in 
# '_node_coord_map' (see below)
#
#-----------------------------------------------------------------------------------

class FETS2D4Q(FETSEval):

    debug_on = True

    # Dimensional mapping
    dim_slice = slice(0, 2)

    # Order of node positions for the formulation of shape function
    #
    dof_r = Array(value = [[-1, -1], [1, -1], [1, 1], [-1, 1]])
    geo_r = Array(value = [[-1, -1], [1, -1], [1, 1], [-1, 1]])

    n_e_dofs = Int(8)
    t = Float(1.0, label = 'thickness')

    # Integration parameters
    #
    ngp_r = Int(2)
    ngp_s = Int(2)

    # Corner nodes are used for visualization 
    vtk_r = Array(value = [[-1., -1.], [ 1., -1.], [ 1., 1.], [-1., 1.]])
    vtk_cells = [[0, 1, 2, 3]]
    vtk_cell_types = 'Quad'

    #vtk_point_ip_map = [0,1,3,2]
    n_nodal_dofs = Int(2)

    #---------------------------------------------------------------------
    # Method required to represent the element geometry
    #---------------------------------------------------------------------
    def get_N_geo_mtx(self, r_pnt):
        '''
        Return the value of shape functions for the specified local coordinate r
        '''
        cx = np.array(self.geo_r, dtype = 'float_')
        Nr = np.array([[1 / 4. * (1 + r_pnt[0] * cx[i, 0]) * (1 + r_pnt[1] * cx[i, 1])
                      for i in range(0, 4) ]])
        return Nr

    def get_dNr_geo_mtx(self, r_pnt):
        '''
        Return the matrix of shape function derivatives.
        Used for the conrcution of the Jacobi matrix.

        @TODO - the B matrix is used
        just for uniaxial bar here with a trivial differential
        operator.
        '''
        #cx = self._node_coord_map
        cx = np.array(self.geo_r, dtype = 'float_')
        dNr_geo = np.array([[ 1 / 4. * cx[i, 0] * (1 + r_pnt[1] * cx[i, 1]) for i in range(0, 4) ],
                          [ 1 / 4. * cx[i, 1] * (1 + r_pnt[0] * cx[i, 0]) for i in range(0, 4) ]])
        return dNr_geo

    #---------------------------------------------------------------------
    # Method delivering the shape functions for the field variables and their derivatives
    #---------------------------------------------------------------------
    def get_N_mtx(self, r_pnt):
        '''
        Returns the matrix of the shape functions used for the field approximation
        containing zero entries. The number of rows corresponds to the number of nodal
        dofs. The matrix is evaluated for the specified local coordinate r.
        '''
        Nr_geo = self.get_N_geo_mtx(r_pnt)
        I_mtx = np.identity(self.n_nodal_dofs, float)
        N_mtx_list = [I_mtx * Nr_geo[0, i] for i in range(0, Nr_geo.shape[1])]
        N_mtx = np.hstack(N_mtx_list)
        return N_mtx

    def get_dNr_mtx(self, r_pnt):
        '''
        Return the derivatives of the shape functions
        '''
        return self.get_dNr_geo_mtx(r_pnt)

    def get_B_mtx(self, r_pnt, X_mtx):
        J_mtx = self.get_J_mtx(r_pnt, X_mtx)
        dNr_mtx = self.get_dNr_mtx(r_pnt)
        dNx_mtx = np.dot(inv(J_mtx), dNr_mtx)
        Bx_mtx = np.zeros((3, 8), dtype = 'float_')
        for i in range(0, 4):
            Bx_mtx[0, i * 2] = dNx_mtx[0, i]
            Bx_mtx[1, i * 2 + 1] = dNx_mtx[1, i]
            Bx_mtx[2, i * 2] = dNx_mtx[1, i]
            Bx_mtx[2, i * 2 + 1] = dNx_mtx[0, i]
        return Bx_mtx


#----------------------- example --------------------

def example_with_new_domain():
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainListField, \
        RTraceDomainListInteg, TLoop, TLine, BCSlice, FEGrid, FERefinementGrid, \
        FEDomain
    from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic

    fets_eval = FETS2D4Q(mats_eval = MATS2DElastic())

    from mathkit.mfn import MFnLineArray

    #===========================================================================
    # Geometry
    #===========================================================================
    L1 = 0.2
    L2 = 0.2
    alpha = math.pi / 2.0 / 3.0
    d = 0.01

    def gt1(points):
        x, y = points.T
        x_ = x * L2
        y_ = y * d
        x_ -= L2
        T = np.array([[ math.cos(alpha), math.sin(alpha)],
                      [ -math.sin(alpha), math.cos(alpha)]], dtype = 'f')
        
        points = np.dot(np.c_[ x_, y_ ], T)
        return points
     

    def gt2(points):
        x, y = points.T
        x_ = x * L2
        y_ = y * d
        return np.c_[x_, y_]
     

    fe_domain = FEDomain()
    
    fe_rg1 = FERefinementGrid(name = 'rg1',
                              fets_eval = fets_eval,
                              domain = fe_domain)

    # Discretization
    fe_grid1 = FEGrid(coord_max = (1., 1.),
                      shape = (10, 5),
                      fets_eval = fets_eval,
                     geo_transform = gt1,
                     level = fe_rg1)


    fe_rg2 = FERefinementGrid(name = 'rg2',
                              fets_eval = fets_eval,
                              domain = fe_domain)

    # Discretization
    fe_grid2 = FEGrid(coord_max = (1., 1.),
                     shape = (10, 5),
                     fets_eval = fets_eval,
                     geo_transform = gt2,
                     level = fe_rg2)

    print 'count dofs', fe_domain.n_dofs

    bc_fixed = BCSlice(var = 'u', value = 0., dims = [0, 1],
                       slice = fe_grid1[0, :, 0, :])
    bc_link12 = BCSlice(var = 'u',
                        value = 0.,
                        dims = [0, 1],
                        slice = fe_grid1[-1, :, -1, :],
                        link_coeffs = [1.0, 1.0],
                        link_dims = [0, 1],
                        link_slice = fe_grid2[0, :, 0, :])
    bc_load2 = BCSlice(var = 'u', value = -0.01, dims = [1],
                      slice = fe_grid2[-1, -1, -1, -1 ])

    mf = MFnLineArray(#xdata = arange(10),
                       ydata = np.array([0, 1, 2, 3]))

    right_dof = 2
    tstepper = TS(sdomain = fe_domain,
                   bcond_list = [bc_fixed,
                                 bc_link12,
                                 bc_load2,
                                  ],
         rtrace_list = [
                     RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                               var_y = 'F_int', idx_y = right_dof,
                               var_x = 'U_k', idx_x = right_dof,
                               record_on = 'update'),
                         RTraceDomainListField(name = 'Stress' ,
                         var = 'sig_app', idx = 0,
                         #position = 'int_pnts',
                         record_on = 'update'),
                     RTraceDomainListField(name = 'Displacement' ,
                                    var = 'u', idx = 0,
                                    record_on = 'update',
                                    warp = True),
                     RTraceDomainListField(name = 'Strain energy' ,
                                    var = 'strain_energy', idx = 0,
                                    record_on = 'update',
                                    warp = False),
                     RTraceDomainListInteg(name = 'Integ strain energy' ,
                                    var = 'strain_energy', idx = 0,
                                    record_on = 'update',
                                    warp = False),
                ]
            )

    # Add the time-loop control
    tloop = TLoop(tstepper = tstepper, KMAX = 300, tolerance = 1e-4,
                   tline = TLine(min = 0.0, step = 1.0, max = 1.0))

    tloop.eval()

    # Put the whole thing into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp(ibv_resource = tloop)
    app.main()

if __name__ == '__main__':
    example_with_new_domain()
