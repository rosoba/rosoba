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

from ibvpy.api import \
    IBVModel, FEDomain, FERefinementGrid, FEGrid, BCSlice, BCDofGroup, \
    RTraceGraph, RTraceDomainListField, RTraceDomainListInteg, \
    TLoop, TStepper as TS, TLine

from mathkit.mfn import MFnLineArray

import numpy as np

from ibvpy.fets.fets2D import FETS2D4Q8U

#from ibvpy.fets.fets1D5 import FETS1D52L4ULRH

from e09_fets_crack import FETS1D5t2L4ULRH

from ibvpy.mats.mats1D import \
    MATS1DDamage, MATS1DPlastic, MATS1DElastic

from ibvpy.mats.mats1D5.mats1D5_bond import \
    MATS1D5Bond

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

    #===========================================================================
    # Discretization parameters
    #===========================================================================
    n_z = Int(1, desc = 'number of elements in the thickness direction')
    n_x = Int(1, desc = 'number of elements in the length direction of a plate')

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
    #-----------------------------------------------------------
    # material parameters:  E-modulus, Poisson's ratio
    #-----------------------------------------------------------
    E_concrete = Float(28e5, input = True)

    nu_concrete = Float(0.2, input = True)

    E_composite = Float(28e5, input = True)

    nu_composite = Float(0.2, input = True)

    E_steel = Float(28e5, input = True)

    nu_steel = Float(0.2, input = True)


    stiffness_crack = Float (34000 * 0.03 * 0.03, input = True)
    A_fiber = Float (1., input = True)
    E_fiber = Float (1., input = True)
    stiffness_fiber = 12


    tau_max_crack = Float (1.* math.sqrt(math.pi) * 2 * math.pi, input = True)
    G_crack = Float (1., input = True)
    u_max = Float (0.23, input = True)
    f_max = Float (0.2, input = True)





    #-----------------
    # fets:
    #-----------------






    # use quadratic serendipity elements
    #
    fets_crack = Property(Instance(FETSEval),
                          depends_on = 'E_crack,nu_crack')
    @cached_property
    def _get_fets_crack (self):
        return FETS1D5t2L4ULRH(mats_eval = MATS1D5Bond(mats_phase1 = MATS1DElastic(E = 0),
                                                      mats_phase2 = MATS1DElastic(E = 0),
                                                      mats_ifslip = MATS1DElastic(E = 1e+6),
                                                      mats_ifopen = MATS1DElastic(E = 1e+10)))

    fets_grouting = Property(Instance(FETSEval),
                                 depends_on = 'E_concrete,nu_concrete')
    @cached_property
    def _get_fets_grouting(self):
        return FETS2D4Q8U(mats_eval = MATS2DElastic(
                                                      E = self.E_concrete,
                                                      nu = self.nu_concrete,
                                                      ))


    fets_concrete = Property(Instance(FETSEval),
                                depends_on = 'E_,nu')
    @cached_property
    def _get_fets_concrete(self):
        return FETS2D4Q8U(mats_eval = MATS2DElastic(
                                                      E = self.E_concrete,
                                                      nu = self.nu_concrete,
                                                      ))


    geo_r = Array(value = [[-1, -1], [1, -1], [1, 1], [-1, 1]], dtype = 'f')

    def N_transform(self, r, X):
        cx = np.array(self.geo_r, dtype = 'float_')
        Nr = np.array([1 / 4. * (1 + r[:, 0] * cx[i, 0]) * (1 + r[:, 1] * cx[i, 1])
                      for i in range(0, 4) ])
        return np.dot(Nr.T, X)




    def gt_plate_left(self, points):
        X1 = np.array([[-self.L2 - 0.0001, 0], [0, 0], [0, self.d], [-self.L2 - 0.0001, self.d]], dtype = 'f')
        T = np.array([[ math.cos(self.alpha), math.sin(self.alpha)],
                      [ -math.sin(self.alpha), math.cos(self.alpha)]], dtype = 'f')
        X1 = np.dot(X1, T)
        return self.N_transform(points, X1)

    def gt_crack_left(self, points):
        X = np.array([[-0.0001, 0], [0, 0], [0, self.d], [-0.0001, self.d]], dtype = 'f')
        T = np.array([[ math.cos(self.alpha), math.sin(self.alpha)],
                      [ -math.sin(self.alpha), math.cos(self.alpha)]], dtype = 'f')
        X = np.dot(X, T)
        return self.N_transform(points, X)


    def gt_grouting(self, points):
        X2 = np.array([[0, 0], [self.h, 0], [self.h, self.d],
                       [-self.d * math.sin(self.alpha), self.d * math.cos(self.alpha)]], dtype = 'f')
        return self.N_transform(points, X2)

    def gt_crack_right(self, points):
        X = np.array([[self.h, 0], [self.h + 0.0001, 0], [self.h + 0.0001, self.d], [self.h, self.d]], dtype = 'f')
        return self.N_transform(points, X)


    def gt_plate_right(self, points):
        X3 = np.array([[self.h + 0.0001, 0], [self.L1 + self.h, 0],
                       [self.L1 + self.h, self.d], [self.h + 0.0001, self.d]], dtype = 'f')
        return self.N_transform(points, X3)



    fe_domain = Property
    @cached_property
    def _get_fe_domain(self):
        return FEDomain()

    fe_rg1 = Property
    @cached_property
    def _get_fe_rg1(self):
        return FERefinementGrid(name = 'rg1',
                                fets_eval = self.fets_concrete,
                                domain = self.fe_domain)

    fe_grid1 = Property
    @cached_property
    def _get_fe_grid1(self):
        return FEGrid(coord_min = (-1., -1.),
                          coord_max = (1., 1.),
                          shape = (self.n_x, self.n_z),
                          fets_eval = self.fets_concrete,
                          geo_transform = self.gt_plate_left,
                          level = self.fe_rg1)

    fe_crack_left = Property
    @cached_property
    def _get_fe_crack_left (self):
        return FERefinementGrid(name = 'rg4',
                                fets_eval = self.fets_crack,
                                domain = self.fe_domain)

    fe_grid_crack_left = Property
    @cached_property
    def _get_fe_grid_crack_left (self):
        return FEGrid(coord_min = (-1., -1.),
                          coord_max = (1., 1.),
                          shape = (self.n_x, self.n_z),
                          fets_eval = self.fets_crack,
                          geo_transform = self.gt_crack_left,
                          level = self.fe_crack_left)



    fe_rg2 = Property
    @cached_property
    def _get_fe_rg2(self):
        return FERefinementGrid(name = 'rg2',
                                fets_eval = self.fets_concrete,
                                domain = self.fe_domain)

    fe_grid2 = Property
    @cached_property
    def _get_fe_grid2(self):
        return FEGrid(coord_min = (-1., -1.),
                          coord_max = (1., 1.),
                          shape = (1, self.n_z),
                          fets_eval = self.fets_concrete,
                          geo_transform = self.gt_grouting,
                          level = self.fe_rg2)


    fe_crack_right = Property
    @cached_property
    def _get_fe_crack_right (self):
        return FERefinementGrid(name = 'rg5',
                                fets_eval = self.fets_crack,
                                domain = self.fe_domain)

    fe_grid_crack_right = Property
    @cached_property
    def _get_fe_grid_crack_right (self):
        return FEGrid(coord_min = (-1., -1.),
                          coord_max = (1., 1.),
                          shape = (self.n_x, self.n_z),
                          fets_eval = self.fets_crack,
                          geo_transform = self.gt_crack_right,
                          level = self.fe_crack_right)

    fe_rg3 = Property
    @cached_property
    def _get_fe_rg3(self):
        return FERefinementGrid(name = 'rg3',
                                fets_eval = self.fets_concrete,
                                domain = self.fe_domain)

    fe_grid3 = Property
    @cached_property
    def _get_fe_grid3(self):
        return FEGrid(coord_min = (-1., -1.),
                      coord_max = (1., 1.),
                      shape = (self.n_x, self.n_z),
                      fets_eval = self.fets_concrete,
                      geo_transform = self.gt_plate_right,
                      level = self.fe_rg3)




    tloop = Property()
    @cached_property
    def _get_tloop(self):

        self.fe_grid1
        self.fe_grid2
        self.fe_grid3
        self.fe_crack_left
        self.fe_crack_right
        bc_fixed_botton = BCSlice(var = 'u', value = 0., dims = [0, 1],
                           slice = self.fe_grid1[:, 0, :, 0])
        bc_fixed_top = BCSlice(var = 'u', value = 0., dims = [0, 1],
                               slice = self.fe_grid1[:, -1, :, -1])


        bc_link_1_crack_left = BCDofGroup(var = 'u',
                                      value = 0.,
                                      dims = [0, 1],
                                      get_dof_method = self.fe_grid1.get_right_dofs,
                                      #slice = self.fe_grid1[-1, :, -1, :],
                                      link_coeffs = [1.0],
                                      link_dims = [0, 1],
                                      get_link_dof_method = self.fe_grid_crack_left.get_bottom_dofs,
                                      #link_slice = self.fe_grid_crack_left[:, -1, :, -1]
                                      )

        bc_link_crack_left_2 = BCDofGroup(var = 'u',
                                      value = 0.,
                                      dims = [0, 1],
                                      get_dof_method = self.fe_grid_crack_left.get_top_dofs,
                                      #slice = self.fe_grid_crack_left[:, 0, :, 0],
                                      link_coeffs = [1.0],
                                      link_dims = [0, 1],
                                      get_link_dof_method = self.fe_grid2.get_left_dofs,
                                      #link_slice = self.fe_grid2[0, :, 0, :]
                                      )


        bc_link_2_crack_right = BCDofGroup(var = 'u',
                                       value = 0.,
                                       dims = [0, 1],
                                       get_dof_method = self.fe_grid2.get_right_dofs,
                                       #slice = self.fe_grid2[-1, :, -1, :],
                                       link_coeffs = [1.0],
                                       link_dims = [0, 1],
                                       get_link_dof_method = self.fe_grid_crack_right.get_bottom_dofs,
                                       #link_slice = self.fe_grid_crack_right[:, -1, :, -1]
                                       )

        bc_link_crack_right_3 = BCDofGroup(var = 'u',
                                        value = 0.,
                                        dims = [0, 1],
                                        get_dof_method = self.fe_grid_crack_right.get_top_dofs,
                                        #slice = self.fe_grid_crack_right[:, 0, :, 0],
                                        link_coeffs = [1.0],
                                        link_dims = [0, 1],
                                        get_link_dof_method = self.fe_grid3.get_left_dofs,
                                        #link_slice = self.fe_grid3[0, :, 0, :]
                                        )

        mf = MFnLineArray(xdata = np.array([0, 0.1, 0.6, 1], dtype = 'f'),
                          ydata = np.array([0, 0.4, -0.5, 1], dtype = 'f'))

        bc_load2 = BCSlice(var = 'u', value = -0.005, dims = [1], time_function = mf.get_value,
                          slice = self.fe_grid3[-1, -1, -1, -1 ])

        load_dofs = self.fe_grid3[-1, -1, :, -1 ].dofs[:, :, 1].flatten()

        tstepper = TS(sdomain = self.fe_domain,
                       bcond_list = [bc_fixed_botton,
                                     bc_fixed_top,
                                     bc_link_1_crack_left,
                                     bc_link_crack_left_2,
                                     bc_link_2_crack_right,
                                     bc_link_crack_right_3,
                                     bc_load2,
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
                       tline = TLine(min = 0.0, step = 1., max = 1.0))

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

    fbt = FoldedBondTest(n_x = 5, n_z = 3)

    fbt.tloop.eval()

    # Put the whole thing into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp(ibv_resource = fbt)
    app.main()

