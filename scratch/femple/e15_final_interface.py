'''
Created on 04.10.2012

@author: demian
'''

import sys
print sys.path

from etsproxy.traits.api import \
     Array, Bool, Callable, Enum, Float, HasTraits, Interface, implements, \
     Instance, Int, Trait, Str, Enum, Callable, List, TraitDict, Any, \
     on_trait_change, Tuple, WeakRef, Delegate, Property, cached_property, \
     Constant

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

from ibvpy.fets.fets2D import FETS2D4Q, FETS2D4Q8U


from e09_fets_crack import FETS1D5t2L4ULRH
from e09_fets_crack_without_T import FETS1D5t2L4ULRH_oT as FETS1D5t2L4ULRH_woT

from ibvpy.mats.mats1D import \
    MATS1DDamage, MATS1DPlastic, MATS1DElastic

from ibvpy.mats.mats1D5.mats1D5_bond import \
    MATS1D5Bond

from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
from ibvpy.mats.mats2D.mats2D_plastic.mats2D_plastic import MATS2DPlastic

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

    alpha = Float(math.pi / 2, input = True)

    #===========================================================================
    # Reinforcement properties
    #===========================================================================
    E_tex = Float(50e+3, input = True)
    A_tex = Float(0.89, input = True)
    n_tex = Int(9, input = True)

    K_tex = Property()
    def _get_K_tex(self):
        return self.E_tex * self.A_tex * self.n_tex

    #===========================================================================
    # Concrete properties
    #===========================================================================
    E_c = Float(33e+3, input = True)
    nu_c = Float(0.2, input = True)

    #===========================================================================
    # Interface properties
    #===========================================================================
    E_crack = Float(33e+3, input = True)
    nu_crack = Float(0.2, input = True)
    sig_max_crack = Float(6.0, input = True)
    eps_0 = Property
    def _get_eps_0(self):
        return self.sig_max_crack / self.E_crack
    eps_f = Float(.05, input = True)

    G_tex = Float(1000., input = True)
    u_max = Float(5, input = True)
    f_max = Float(10, input = True)

    #-----------------
    # fets
    #-----------------

    fets_tex = Property(Instance(FETSEval),
                          depends_on = 'E_tex')
    @cached_property
    def _get_fets_tex (self):
        return FETS1D5t2L4ULRH_woT(mats_eval = MATS1D5Bond(mats_phase1 = MATS1DElastic(E = self.E_tex),
                                                      mats_phase2 = MATS1DElastic(E = self.E_tex),
                                                      mats_ifopen = MATS1DElastic(E = 1e+12),
                                                      mats_ifslip = MATS1DPlastic(E = self.G_tex,
                                                                                  sigma_y = 0.1,
                                                                                  K_bar = self.K_tex,
                                                                                  H_bar = 0.)))
    fets_crack = Property(Instance(FETSEval),
                          depends_on = 'E_crack,nu_crack')
    @cached_property
    def _get_fets_crack (self):
        return FETS1D5t2L4ULRH(mats_eval = MATS1D5Bond(mats_phase1 = MATS1DElastic(E = 0.0),
                                                      mats_phase2 = MATS1DElastic(E = 0.0),
                                                      mats_ifopen = MATS1DDamage(E = self.E_crack,
                                                                                 epsilon_0 = self.eps_0,
                                                                                 epsilon_f = self.eps_f),
                                                      mats_ifslip = MATS1DElastic(E = 1e+5)))


    fets_plate = Property(Instance(FETSEval),
                                 depends_on = 'E_c,nu_c')
    @cached_property
    def _get_fets_plate(self):
        return FETS2D4Q8U(mats_eval = MATS2DElastic(
                                                  E = self.E_c,
                                                  nu = self.nu_c,
                                                  ))

    geo_r = Array(value = [[-1, -1], [1, -1], [1, 1], [-1, 1]], dtype = 'f')

    def N_transform(self, r, X):
        cx = np.array(self.geo_r, dtype = 'float_')
        Nr = np.array([1 / 4. * (1 + r[:, 0] * cx[i, 0]) * (1 + r[:, 1] * cx[i, 1])
                      for i in range(0, 4) ])
        return np.dot(Nr.T, X)

    def gt_crack(self, points):
        X1 = np.array([[0, 0], [6, 0], [6, 2], [0, 2]  ], dtype = 'f')
        T = np.array([[ math.cos(self.alpha), math.sin(self.alpha)],
                      [ -math.sin(self.alpha), math.cos(self.alpha)]], dtype = 'f')
        X1 = np.dot(X1, T)
        return self.N_transform(points, X1)

    def gt_tex_c(self, points):
        X1 = np.array([[-2, 2.5], [0, 2.5], [0, 3.5], [-2, 3.5]], dtype = 'f')
        return self.N_transform(points, X1)

    def gt_tex_p(self, points):
        X1 = np.array([[0, 2.5], [138, 2.5], [138, 3.5], [0, 3.5]], dtype = 'f')
        return self.N_transform(points, X1)

    def gt_plate(self, points):
        X = np.array([[0, 0], [138, 0], [138, 6], [0, 6]], dtype = 'f')
        return self.N_transform(points, X)


    fe_domain = Property
    @cached_property
    def _get_fe_domain(self):
        return FEDomain()

    fe_rg_crack = Property
    @cached_property
    def _get_fe_rg_crack(self):
        return FERefinementGrid(name = 'crack',
                                fets_eval = self.fets_crack,
                                domain = self.fe_domain)

    fe_grid_crack = Property
    @cached_property
    def _get_fe_grid_crack(self):
        return FEGrid(coord_min = (-1., -1.),
                          coord_max = (1., 1.),
                          shape = (8, 1),
                          fets_eval = self.fets_crack,
                          geo_transform = self.gt_crack,
                          level = self.fe_rg_crack)

    fe_rg_tex_c = Property
    @cached_property
    def _get_fe_rg_tex_c(self):
        return FERefinementGrid(name = 'tex_c',
                                fets_eval = self.fets_tex,
                                domain = self.fe_domain)

    fe_grid_tex_c = Property
    @cached_property
    def _get_fe_grid_tex_c(self):
        return FEGrid(coord_min = (-1., -1.),
                          coord_max = (1., 1.),
                          shape = (1, 1),
                          fets_eval = self.fets_tex,
                          geo_transform = self.gt_tex_c,
                          level = self.fe_rg_tex_c)

    fe_rg_tex_p = Property
    @cached_property
    def _get_fe_rg_tex_p(self):
        return FERefinementGrid(name = 'tex_p',
                                fets_eval = self.fets_tex,
                                domain = self.fe_domain)

    fe_grid_tex_p = Property
    @cached_property
    def _get_fe_grid_tex_p(self):
        return FEGrid(coord_min = (-1., -1.),
                          coord_max = (1., 1.),
                          shape = (5, 1),
                          fets_eval = self.fets_tex,
                          geo_transform = self.gt_tex_p,
                          level = self.fe_rg_tex_p)


    fe_rg_plate = Property
    @cached_property
    def _get_fe_rg_plate(self):
        return FERefinementGrid(name = 'plate',
                                fets_eval = self.fets_plate,
                                domain = self.fe_domain)

    fe_grid_plate = Property
    @cached_property
    def _get_fe_grid_plate(self):
        return FEGrid(coord_min = (-1., -1.),
                          coord_max = (1., 1.),
                          shape = (5, 8),
                          fets_eval = self.fets_plate,
                          geo_transform = self.gt_plate,
                          level = self.fe_rg_plate)

    tloop = Property()
    @cached_property
    def _get_tloop(self):
        #self.fe_grid_crack
        #self.fe_grid_tex_c
        self.fe_grid_tex_p
        self.fe_grid_plate
#        bc_fixed_crack = BCDofGroup(var = 'u', value = 0., dims = [0, 1],
#                                    get_dof_method = self.fe_grid_crack.get_top_dofs)
#                                    #slice = self.fe_grid_crack[ :, -1, :, -1])
#        bc_fixed_plate = BCSlice(var = 'u', value = 0., dims = [0, 1],
#                                    slice = self.fe_grid_plate[ 0, 0, 0, 0])
        bc_fixed_plate = BCSlice(var = 'u', value = 0., dims = [0, 1],
                                 slice = self.fe_grid_plate[ :, 0, :, 0])
#        bc_link_crack_plate = BCSlice(var = 'u',
#                                      value = 0.,
#                                      dims = [0, 1],
#                                      slice = self.fe_grid_crack[:, 0, :, 0],
#                                      link_coeffs = [1.0, 1.0],
#                                      link_dims = [0, 1],
#                                      link_slice = self.fe_grid_plate[0, : , 0, :],
#                                      )
#        bc_link_crack_tex_c = BCSlice(var = 'u',
#                                    value = 0.,
#                                    dims = [0, 1],
#                                    slice = self.fe_grid_crack[3, 0, :, -1],
#                                    link_coeffs = [1.0, 1.0],
#                                    link_slice = self.fe_grid_tex_c[0, 0, :, -1]
#                                    )

        bc_link_plate_tex_p = BCSlice(var = 'u',
                                    value = 0.,
                                    dims = [0, 1],
                                    slice = self.fe_grid_tex_p[:, -1, :, -1],
                                    link_coeffs = [1.0, 1.0],
                                    link_slice = self.fe_grid_plate[:, 3, :, -1]
                                    )

        bc_fix_tex_p = BCSlice(var = 'u',
                               value = 0.,
                               dims = [0, 1],
                               slice = self.fe_grid_tex_p[:, 0, :, 0],
                               )

        mf = MFnLineArray(xdata = np.array([0, 0.1, 0.6, 1], dtype = 'f'),
                          ydata = np.array([0, 0.4, -0.5, 1], dtype = 'f'))

        bc_load = BCSlice(var = 'u', value = -10.0, dims = [1], #time_function = mf.get_value,
                          slice = self.fe_grid_plate[-1, -1, -1, -1 ])
#        bc_load = BCSlice(var = 'u', value = 10.0, dims = [0], #time_function = mf.get_value,
#                          slice = self.fe_grid_tex_p[-1, 0, -1, 0 ])

        load_dofs = self.fe_grid_plate[-1, -1, -1, -1 ].dofs[:, :, 1].flatten()

        tstepper = TS(sdomain = self.fe_domain,
                      bcond_list = [#bc_fixed_crack,
                                    bc_fixed_plate,
                                    #bc_link_crack_plate,
                                    #bc_link_crack_tex_c,
                                    bc_link_plate_tex_p,
                                    bc_fix_tex_p,
                                    bc_load
                                   ],
             rtrace_list = [
                         RTraceGraph(name = 'Fi,right over u_right (iteration)' ,
                                   var_y = 'F_int', idx_y = load_dofs[0],
                                   var_x = 'U_k', idx_x = load_dofs[0],
                                   transform_x = '-x',
                                   transform_y = '-y',
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
                       tline = TLine(min = 0.0, step = 0.1, max = 1))

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

    fbt.tloop.eval()

    # Put the whole thing into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp(ibv_resource = fbt)
    app.main()

