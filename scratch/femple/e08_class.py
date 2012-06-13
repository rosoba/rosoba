'''
Created on 08.06.2012

@author: demian
'''

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

from ibvpy.api import IBVModel, FEDomain, FERefinementGrid
from ibvpy.mesh.fe_grid import \
     FEGrid


import numpy as np

from scipy.linalg import \
     inv
from ibvpy.fets.fets2D5.fets2D58h20u import \
    FETS2D58H20U

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
from ibvpy.fets.fets_eval import FETSEval
import math

class FoldedBondTest(IBVModel):
    '''Idealization of the test for the characterization of
    bond behavior within crease line of a folded plate. 
    '''
    #===========================================================================
    # Geometry
    #===========================================================================
    L1 = Float(0.2, desc = 'Length of the left part')
    L2 = Float(0.2, desc = 'Length of the right part')
    alpha = Float(math.pi / 2.0 / 3.0, desc = 'Fold angle')
    d = Float(0.01, desc = 'thickness of the plate')
    h = Float(0.01, desc = 'width of the plate')

    geo_r = Float([[-1, -1], [1, -1], [1, 1], [-1, 1]])
    
    
    N_transform = Property
    def _get_N_transform(self,r, X):
        cx = Float(self.geo_r)
        Nr = np.array([1 / 4. * (1 + r[:, 0] * cx[i, 0]) * (1 + r[:, 1] * cx[i, 1])
                      for i in range(0, 4) ])
        return np.dot(Nr.T, X)
    
    gt1 = Property
    def _get_gt1(self,points):
        X1 = np.array([[-self.L2, 0], [0, 0], [0, self.d], [-self.L2, self.d]], dtype = 'f')
        T = np.array([[ math.cos(self.alpha), math.sin(self.alpha)],
                      [ -math.sin(self.alpha), math.cos(self.alpha)]], dtype = 'f')
        X1 = np.dot(X1, T)
        return self.N_transform(points, X1)
 
    
    
    
    #===========================================================================
    # Discretization parameters
    #===========================================================================
    n_z = Int(4, desc = 'number of elements in the thickness direction')
    n_x = Int(5, desc = 'number of elements in the length direction of a plate')

    view = View(Item('L1', label = 'length of fixed part'),
                Item('L2', label = 'length of the loaded part'),
                Item('alpha'),
                width = 0.2,
                height = 0.3
                )

    #----------------------------------------------------------------------------------
    # mats_eval
    #----------------------------------------------------------------------------------

    # age of the plate at the time of testing
    # NOTE: that the same phi-function is used independent of age. This assumes a 
    # an afine/proportional damage evolution for different ages. 
    #
    age = Int(28, input = True)

    # composite E-modulus 
    #
    E_ = Float(28e5, input = True)

    # Poisson's ratio 
    #
    nu = Float(0.2, input = True)

    # @todo: for mats_eval the information of the unit cell should be used
    # in order to use the same number of microplanes and model version etc...
    #
    plate_mats = Property(Instance(MATS2DElastic),
                          depends_on = 'input_change')
    @cached_property
    def _get_plate_mats(self):

    # used for basic testing:
    # return MATS3DElastic( 
    # E = self.E_c,
    # nu = self.nu )
        return MATS2DElastic(
                                E = self.E_,
                                nu = self.nu,

                                )

    
    #-----------------
    # fets:
    #-----------------

    # use quadratic serendipity elements
    #
    plate_fets = Property(Instance(FETSEval),
                             depends_on = 'input_change')
    @cached_property
    def _get_plate_fets(self):
        return FETS2D58H20U(mats_eval = self.plate_mats)

    
    # used for basic testing:
    #        return FETS3D8H20U( mats_eval = self.mats_eval )
    #        return FETS3D8H( mats_eval = self.mats_eval )

    def peval(self):
        '''
        Evaluate the model and return the array of results specified
        in the method get_sim_outputs.
        '''
        U = self.tloop.eval()

        self.f_w_diagram_center.refresh()
        F_max = max(self.f_w_diagram_center.trace.ydata)

        u_center_top_z = U[ self.center_top_dofs ][0, 0, 2]
        return Float([ u_center_top_z, F_max ])

        fe_domain = Property(depends_on = '+ps_levels, +input')
    @cached_property
    
    def _get_fe_domain(self):
        return FEDomain()

    gt1_fe_level = Property(depends_on = '+ps_levels, +input')
    def _get_gt1_fe_level(self):
        return  FERefinementGrid(name = 'plate patch',
                                 fets_eval = self.plate_fets,
                                 domain = self.fe_domain)
    
    gt1_fe_grid = Property(Instance(FEGrid), depends_on = '+ps_levels, +input')
    @cached_property
    def _get_gt1_fe_grid(self):
        # only a quarter of the plate is simulated due to symmetry:
        fe_grid = FEGrid((coord_min = (-1., -1.),
                          coord_max = (1., 1.),
                          shape = (self.n_x, self.n_z),
                          level = self.plate_fe_level,
                          fets_eval = self.plate_fets,
                          geo_transform = self.gt1,
                          )
        return fe_grid
    
    
    #===========================================================================
    # Boundary conditions
    #===========================================================================
    bc_list = Property(depends_on = '+ps_levels, +input')
    @cached_property
    
    def _get_bc_list(self):
        gt1 = self.gt1_fe_grid
        
        #--------------------------------------------------------------
        # boundary conditions 
        #--------------------------------------------------------------
        # the x-axis correspons to the axis of symmetry along the longitudinal axis of the beam:
        bc_gt1_xy = BCSlice(var = 'u', value = 0., dims = [0,1],
                                 slice = gt1[0, :, 0, :])
        #--------------------------------------------------------------
        # loading
        #--------------------------------------------------------------
        
        # w_max = center displacement:
        w_max = -0.020 # [m]
        f_max = -0.010 / 0.10 # [MN/m]

        #--------------------------------------------------------------
        # nodal displacement applied at top at the center of the plate
        #--------------------------------------------------------------
        # NOTE: the entire symmetry axis (yz)-plane is moved downwards 
        # in order to avoid large indentations at the top nodes
        #
#        bc_center_w = BCSlice( var = 'u', value = w_max, dims = [2], slice = domain[-1, :, :, -1, :, :] )
#        bc_center_f = BCSlice( var = 'f', value = 1.0, dims = [2], slice = domain[-1, :, :, -1, :, :] )

        # apply displacement at all center node (line load)
        #
        bc_end = BCSlice(var = 'u', value = w_max, dims = [1],
                              slice = gt1[-1, -1, -1, -1])

        return [bc_gt1_xy,
                bc_end
#                              bc_center_f,
                ]


        
    tloop = Property(depends_on = 'input_change')
    @cached_property
    def _get_tloop(self):

        #--------------------------------------------------------------
        # ts 
        #--------------------------------------------------------------

        gt1 = self.gt1_fe_grid
        
        # center_top_line_dofs
        #
        ctl_dofs = gt1[-1,-1, -1,-1].dofs[:, :].flatten()

#        # right_bottom_line_dofs
#        #
#        ctl_dofs = domain[0, :, 0, 0, :, 0].dofs[:, :, 2].flatten()
#        print 'ctl_dofs', ctl_dofs

        self.ctl_dofs = np.unique(ctl_dofs)

        # center_dof (used for tracing of the displacement)
        #
        center_dof = self.ctl_dofs[0]

        # force-displacement-diagram 
        # 
        self.f_w_diagram_center = RTraceGraph(name = 'displacement (center) - reaction 2',
                                       var_x = 'U_k'  , idx_x = center_dof,

                                       # nodal displacement at center node
                                       #
#                                       var_y = 'F_int', idx_y = center_dof,

#                                       # line load
#                                       #
                                       var_y = 'F_int', idx_y_arr = self.ctl_dofs,

                                       record_on = 'update',
                                       transform_x = '-x * 1000', # %g * x' % ( fabs( w_max ),),
                                       # due to symmetry the total force sums up from four parts of the beam (2 symmetry axis):
                                       #
                                       transform_y = '-4000. * y')

        ts = TS(
                sdomain = self.fe_domain,
                bcond_list = self.bc_list,
                rtrace_list = [
                             self.f_w_diagram_center,
                             RTraceDomainListField(name = 'Displacement' ,
                                            var = 'u', idx = 0, warp = True),
                             RTraceDomainListField(name = 'Stress' ,
                                            var = 'sig_app', idx = 0, warp = True,
                                            record_on = 'update'),
                             RTraceDomainListField(name = 'Strain' ,
                                        var = 'eps_app', idx = 0, warp = True,
                                        record_on = 'update'),
                             RTraceDomainListField(name = 'Damage' ,
                                        var = 'omega_mtx', idx = 0, warp = True,
                                        record_on = 'update'),
                             RTraceDomainListField(name = 'IStress' ,
                                            position = 'int_pnts',
                                            var = 'sig_app', idx = 0,
                                            record_on = 'update'),
                             RTraceDomainListField(name = 'IStrain' ,
                                            position = 'int_pnts',
                                            var = 'eps_app', idx = 0,
                                            record_on = 'update'),
                              ]
                )

        # Add the time-loop control
        tloop = TLoop(tstepper = ts,

                       # allow only a low tolerance 
                       #
                       KMAX = 50,
                       tolerance = 5e-4,

#                       # allow a high tolerance 
#                       #
#                       KMAX = 50,
#                       tolerance = 0.001,

                       RESETMAX = 0,
                       debug = False,
#                       tline = TLine( min = 0.0, step = 0.05, max = 0.05 )
                       tline = TLine(min = 0.0, step = 0.05, max = 1.0)
                       )

        return tloop

    tloop.eval()

    # Put the whole thing into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp(ibv_resource = tloop)
    app.main()

if __name__ == '__main__':
    FoldedBondTest ()
    fbt = FoldedBondTest()
    fbt.configure_traits()
