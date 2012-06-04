
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

from e05_fets2D4q_B_mapping import FETS2D4Q
import numpy as np
import math

def run():
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
    L2 = -0.2
    alpha = math.pi / 2.0 / 3.0
    d = 0.01

    geo_r = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], dtype = 'f')
    
    def N_transform(r, X):
        cx = np.array(geo_r, dtype = 'float_')
        Nr = np.array([1 / 4. * (1 + r[:, 0] * cx[i, 0]) * (1 + r[:, 1] * cx[i, 1])
                      for i in range(0, 4) ])
        return np.dot(Nr.T, X)
  
    def gt1(points):
        
        X1 = np.array([[0, 0], [L2, 0], [L2, d], [0, d]], dtype = 'f')
        
        T = np.array([[ math.cos(alpha), math.sin(alpha)],
                      [ -math.sin(alpha), math.cos(alpha)]], dtype = 'f')
        
        return np.dot(N_transform(points, X1),T)
        
        
    def gt2(points):
        
        X2 = np.array([[0,0],[0,0],[0,d],[-d*math.sin(alpha),d*math.cos(alpha)]], dtype ='f')
        return N_transform(points, X2)
    
    def gt3(points):
        X3 = np.array([[0, 0], [L1, 0], [L1, d], [0, d]], dtype = 'f')
        return N_transform(points, X3)
        
    fe_domain = FEDomain()
    
    
    #fe_grid3.elem_X_map
    #return

    fe_rg1 = FERefinementGrid(name = 'rg1',
                              fets_eval = fets_eval,
                              domain = fe_domain)

    # Discretization
    fe_grid1 = FEGrid(coord_min =(-1.,-1.),
                      coord_max = (1., 1.),
                      shape = (1, 1),
                      fets_eval = fets_eval,
                      geo_transform = gt1,
                      level = fe_rg1)


    fe_rg2 = FERefinementGrid(name = 'rg2',
                              fets_eval = fets_eval,
                              domain = fe_domain)
    
    # Discretization
    fe_grid2 = FEGrid(coord_min =(-1.,-1.),
                      coord_max = (1., 1.),
                      shape = (1, 1),
                      fets_eval = fets_eval,
                      geo_transform = gt2,
                      level = fe_rg2)
    
    
    fe_rg3 = FERefinementGrid(name = 'rg3',
                              fets_eval = fets_eval,
                              domain = fe_domain)

    # Discretization
    fe_grid3 = FEGrid(coord_min = (-1.,-1.),
                      coord_max =  (1., 1.),
                      shape = (1, 1),
                      fets_eval = fets_eval,
                      geo_transform = gt3,
                      level = fe_rg3)


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
    
    '''
    bc_link23 = BCSlice(var = 'u',
                        value = 0.,
                        dims = [0, 1],
                        slice = fe_grid2[-1, :, -1, :],
                        link_coeffs = [1.0, 1.0],
                        link_dims = [0, 1],
                        link_slice = fe_grid3[0, :, 0, :])
    '''
    
    
    bc_load2 = BCSlice(var = 'u', value = -0.01, dims = [1],
                      slice = fe_grid2[-1, -1, -1, -1 ])

    mf = MFnLineArray(#xdata = arange(10),
                       ydata = np.array([0, 1, 2, 3]))

    right_dof = 2
    tstepper = TS(sdomain = fe_domain,
                   bcond_list = [bc_fixed,
                                 bc_link12,
                                 #bc_link23,
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
    run()
