
import sys
print sys.path

from etsproxy.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property, \
     Constant

import etsproxy.traits.has_traits
etsproxy.traits.has_traits.CHECK_INTERFACES = 2

from etsproxy.traits.ui.api import \
     Item, View, HGroup, ListEditor, VGroup, Group

from etsproxy.traits.ui.menu import \
     NoButtons, OKButton, CancelButton, Action, CloseAction, Menu, \
     MenuBar, Separator

import numpy as np

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

    #===========================================================================
    # Kinematic mapping
    #===========================================================================
    dN_idx_map = Constant(value = ([0, 1, 1, 0],
                                    slice(None)))
    
    B_idx_map = Constant(value = ([0, 1, 2, 2],
                                  slice(None),
                                  [0, 1, 0, 1]))
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
        B_mtx = np.zeros ((3, 4, 2), dtype = 'f')
        
        B_mtx[self.B_idx_map] = dNx_mtx[self.dN_idx_map]
   
        return B_mtx.reshape(3, 8)


#----------------------- example --------------------

def run():
    from ibvpy.api import \
        TStepper as TS, RTraceGraph, RTraceDomainListField, \
        RTraceDomainListInteg, TLoop, \
        TLine, BCDof, IBVPSolve as IS, DOTSEval
    from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
    from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage

    from ibvpy.api import BCDofGroup

    mats_eval = MATS2DElastic()
    fets_eval = FETS2D4Q(mats_eval = mats_eval)
    #fets_eval = FETS2D4Q(mats_eval = MATS2DScalarDamage()) 

    print fets_eval.vtk_node_cell_data

    from ibvpy.mesh.fe_grid import FEGrid
    from ibvpy.mesh.fe_refinement_grid import FERefinementGrid
    from ibvpy.mesh.fe_domain import FEDomain
    from mathkit.mfn import MFnLineArray

    # Discretization
    fe_grid = FEGrid(coord_max = (10., 4., 0.),
                      shape = (10, 3),
                      fets_eval = fets_eval)

    bcg = BCDofGroup(var = 'u', value = 0., dims = [0],
                   get_dof_method = fe_grid.get_left_dofs)
    bcg.setup(None)
    print 'labels', bcg._get_labels()
    print 'points', bcg._get_mvpoints()

    mf = MFnLineArray(#xdata = arange(10),
                       ydata = np.array([0, 1, 2, 3]))

    right_dof = 2
    tstepper = TS(sdomain = fe_grid,
                   bcond_list = [ BCDofGroup(var = 'u', value = 0., dims = [0, 1],
                                               get_dof_method = fe_grid.get_left_dofs),
#                                   BCDofGroup( var='u', value = 0., dims = [1],
#                                  get_dof_method = fe_grid.get_bottom_dofs ),                                  
                         BCDofGroup(var = 'u', value = .005, dims = [1],
                                  time_function = mf.get_value,
                                  get_dof_method = fe_grid.get_right_dofs) ],
         rtrace_list = [
                     RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                               var_y = 'F_int', idx_y = right_dof,
                               var_x = 'U_k', idx_x = right_dof,
                               record_on = 'update'),
                         RTraceDomainListField(name = 'Stress' ,
                         var = 'sig_app', idx = 0,
                         position = 'int_pnts',
                         record_on = 'update'),
#                     RTraceDomainListField(name = 'Damage' ,
#                                    var = 'omega', idx = 0,
#                
#                    record_on = 'update',
#                                    warp = True),
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
                                    #                    RTraceDomainListField(name = 'N0' ,
#                                      var = 'N_mtx', idx = 0,
#                                      record_on = 'update')
                ]
            )

    # Add the time-loop control
    #global tloop
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
    run()
