
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
    run()
