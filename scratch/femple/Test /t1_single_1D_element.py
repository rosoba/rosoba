'''
Created on 10.10.2012

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

from ibvpy.api import \
    IBVModel, FEDomain, FERefinementGrid, FEGrid, BCSlice, BCDofGroup, \
    RTraceGraph, RTraceDomainListField, RTraceDomainListInteg, \
    TLoop, TStepper as TS, TLine

from mathkit.mfn import MFnLineArray

import numpy as np

from ibvpy.fets.fets1D import FETS1D2L

from ibvpy.mats.mats1D import \
    MATS1DDamage, MATS1DPlastic, MATS1DElastic

from ibvpy.fets.fets_eval import FETSEval
import math


class FoldedBondTest(IBVModel):
    '''Idealization of the test for the characterization of
    bond behavior within crease line of a folded plate.
    '''

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

    age = Int(28, input = True)

    #----------------------------------------------------------------------------------
    # material parameters
    #----------------------------------------------------------------------------------
    E = Float(50e3, input = True)
    Sigma_max_crack = Float(23, input = True)
    #-----------------
    # fets
    #-----------------

    fets = Property(Instance(FETSEval),
                          depends_on = 'E')
    @cached_property
    def _get_fets(self):
        return FETS1D2L(mats_eval = MATS1DPlastic(E = self.E,
                                                   sigma_y = self.Sigma_max_crack,
                                                   K_bar = 0.,
                                                   H_bar = 0.))

    geo_r = Array(value = [[-1, -1], [1, -1], [1, 1], [-1, 1]], dtype = 'f')

    def N_transform(self, r, X):
        Nr = np.array([[1. / 2 * r * (r - 1.), 1. - r * r, 1. / 2 * r * (r + 1.)]])
        return np.dot(Nr.T, X)

    def gt_tex(self, points):
        X1 = np.array([[-2, 3], [0, 3]], dtype = 'f')
        return self.N_transform(points, X1)

    fe_domain = Property
    @cached_property
    def _get_fe_domain(self):
        return FEDomain()

    fe_rg = Property
    @cached_property
    def _get_fe_rg(self):
        return FERefinementGrid(name = 'tex',
                                fets_eval = self.fets,
                                domain = self.fe_domain)

    fe_grid = Property
    @cached_property
    def _get_fe_grid(self):
        return FEGrid(coord_max = (1.,),
                          shape = (1,),
                          fets_eval = self.fets,
                          geo_transform = self.gt_tex,
                          level = self.fe_rg)

    tloop = Property()
    @cached_property
    def _get_tloop(self):
        self.fe_grid

        mf = MFnLineArray(xdata = np.array([0, 0.1, 0.6, 1], dtype = 'f'),
                          ydata = np.array([0, 0.4, -0.5, 1], dtype = 'f'))

        bc_fixed = BCSlice(var = 'u',
                           value = 0.,
                           dims = [0, 1],
                           slice = self.fe_grid[-1, -1])

        bc_load = BCSlice(var = 'u', value = 1, dims = [0], #time_function = mf.get_value,
                          slice = self.fe_grid[-1, -1, -1, -1 ])

        load_dofs = self.fe_grid[-1, -1 ].dofs[:, :, 1].flatten()

        tstepper = TS(sdomain = self.fe_domain,
                       bcond_list = [bc_fixed,
                                     bc_load,
                                      ],
             rtrace_list = [
                         RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                                   var_y = 'F_int', idx_y_arr = load_dofs,
                                   var_x = 'U_k', idx_x = load_dofs[-1],
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
                       tline = TLine(min = 0.0, step = 0.1, max = 1.0))

        return tloop

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


if __name__ == '__main__':

    fbt = FoldedBondTest()

    #fbt.tloop.eval()

    # Put the whole thing into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp(ibv_resource = fbt)
    app.main()

