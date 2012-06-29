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
    IBVModel, FEDomain, FERefinementGrid, FEGrid, BCSlice, \
    RTraceGraph, RTraceDomainListField, RTraceDomainListInteg, \
    TLoop, TStepper as TS, TLine

from mathkit.mfn import MFnLineArray

import numpy as np

from ibvpy.fets.fets2D import FETS2D4Q8U
from ibvpy.fets.fets1D5 import FETS1D52L4ULRH

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
    da = 2 * math.sqrt(math.pi)
    tau_max = 0.1 * da * math.pi
    G = 100
    u_max = 0.023
    f_max = 0.2
    
    
    E_concrete = Float(28e5, input = True)

    # Poisson's ratio 
    #
    nu_concrete = Float(0.2, input = True)

    # composite E-modulus 
    #
    E_steel = Float(28e5, input = True)

    # Poisson's ratio 
    #
    nu_steel = Float(0.2, input = True)

    #-----------------
    # fets:
    #-----------------

    # use quadratic serendipity elements
    #
    fets_grouting = Property(Instance(FETSEval),
                                 depends_on = 'E_,nu')
    @cached_property
    def _get_grouting(self):
        return FETS2D4Q8U(mats_eval = MATS2DElastic(
                                                      E = self.E_composite,
                                                      nu = self.nu_composite,
                                                      ))


    fets_concrete = Property(Instance(FETSEval),
                                depends_on = 'E_,nu')
    @cached_property
    def _get_fets_concrete(self):
        return FETS2D4Q8U(mats_eval = MATS2DElastic(
                                                      E = self.E_concrete,
                                                      nu = self.nu_concrete,
                                                      ))    

    fets_steel = Property(Instance(FETSEval),
                             depends_on = 'E_steel,nu_steel')
    @cached_property
    def _get_fets_steel(self):
        return FETS2D4Q8U(mats_eval = MATS2DElastic(
                                                      E = self.E_steel,
                                                      nu = self.nu_steel,
                                                      ))    


    geo_r = Array(value = [[-1, -1], [1, -1], [1, 1], [-1, 1]], dtype = 'f')
        
    def N_transform(self, r, X):
        cx = np.array(self.geo_r, dtype = 'float_')
        Nr = np.array([1 / 4. * (1 + r[:, 0] * cx[i, 0]) * (1 + r[:, 1] * cx[i, 1])
                      for i in range(0, 4) ])
        return np.dot(Nr.T, X)
      
    def gt_plate_left(self, points):
        X1 = np.array([[-self.L2, 0], [0, 0], [0, self.d], [-self.L2, self.d]], dtype = 'f')
        T = np.array([[ math.cos(self.alpha), math.sin(self.alpha)],
                      [ -math.sin(self.alpha), math.cos(self.alpha)]], dtype = 'f')
        X1 = np.dot(X1, T)
        return self.N_transform(points, X1)
        
        
    def gt_grouting(self, points):
        X2 = np.array([[0, 0], [self.h, 0], [self.h, self.d],
                       [-self.d * math.sin(self.alpha), self.d * math.cos(self.alpha)]], dtype = 'f')
        return self.N_transform(points, X2)
    
    def gt_plate_right(self, points):
        X3 = np.array([[self.h, 0], [self.L1 + self.h, 0],
                       [self.L1 + self.h, self.d], [self.h, self.d]], dtype = 'f')
        return self.N_transform(points, X3)
        
    def gt_buttstrap_bottom(self, points):
        X4 = np.array([[-self.L2, -self.d], [0, -self.d], [0, 0], [-self.L2, 0]], dtype = 'f')
        T = np.array([[ math.cos(self.alpha), math.sin(self.alpha)],
                      [ -math.sin(self.alpha), math.cos(self.alpha)]], dtype = 'f')
        X4 = np.dot(X4, T)
        print 'X4=',X4
        print 'X4erste Zeile',X4[0,0]
        print 'X4erster Wert', X4[0,0],[0,0]
        return self.N_transform(points, X4)
    
    def gt_buttstrap_top(self, points):
        X5 = np.array([[-self.L2, self.d], [0, self.d],
                       [0, 2 * self.d], [-self.L2, 2 * self.d]], dtype = 'f')
        T = np.array([[ math.cos(self.alpha), math.sin(self.alpha)],
                      [ -math.sin(self.alpha), math.cos(self.alpha)]], dtype = 'f')
        X5 = np.dot(X5, T)
        return self.N_transform(points, X5)
        
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
    
    
    fe_rg4 = Property
    @cached_property
    def _get_fe_rg4(self):
        return FERefinementGrid(name = 'rg4',
                                  fets_eval = self.fets_steel,
                                  domain = self.fe_domain)
    
    fe_grid4 = Property
    @cached_property
    def _get_fe_grid4(self):
        return FEGrid(coord_min = (-1., -1.),
                          coord_max = (1., 1.),
                          shape = (self.n_x, 1),
                          fets_eval = self.fets_steel,
                          geo_transform = self.gt_buttstrap_bottom,
                          level = self.fe_rg4)
        
    fe_rg5 = Property
    @cached_property
    def _get_fe_rg5(self):
        return FERefinementGrid(name = 'rg5',
                                  fets_eval = self.fets_steel,
                                  domain = self.fe_domain)
    
    fe_grid5 = Property
    @cached_property
    def _get_fe_grid5(self):
        return FEGrid(coord_min = (-1., -1.),
                          coord_max = (1., 1.),
                          shape = (self.n_x, 1),
                          fets_eval = self.fets_steel,
                          geo_transform = self.gt_buttstrap_top,
                          level = self.fe_rg5)
    
    tloop = Property()
    @cached_property
    def _get_tloop(self):
    
        self.fe_grid1
        self.fe_grid2
        self.fe_grid3
        self.fe_grid4
        self.fe_grid5
        print 'count dofs', self.fe_domain.n_dofs
       
        bc_fixed = BCSlice(var = 'u', value = 0., dims = [0, 1],
                           slice = self.fe_grid4[:, 0, :, 0])
        
        bc_link14 = BCSlice(var = 'u',
                            value = 0.,
                            dims = [0, 1],
                            slice = self.fe_grid4[:, -1, :-1, -1],
                            link_coeffs = [1.0, 1.0],
                            link_dims = [0, 1],
                            link_slice = self.fe_grid1[:, 0, :-1, 0])
      
        bc_link15 = BCSlice(var = 'u',
                            value = 0.,
                            dims = [0, 1],
                            slice = self.fe_grid5[:, 0, :-1, 0],
                            link_coeffs = [1.0, 1.0],
                            link_dims = [0, 1],
                            link_slice = self.fe_grid1[:, -1, :-1, -1])
         
    
        bc_link12 = BCSlice(var = 'u',
                            value = 0.,
                            dims = [0, 1],
                            slice = self.fe_grid1[-1, :, -1, :],
                            link_coeffs = [1.0, 1.0],
                            link_dims = [0, 1],
                            link_slice = self.fe_grid2[0, :, 0, :])
        
        bc_link23 = BCSlice(var = 'u',
                            value = 0.,
                            dims = [0, 1],
                            slice = self.fe_grid2[-1, :, -1, :],
                            link_coeffs = [1.0, 1.0],
                            link_dims = [0, 1],
                            link_slice = self.fe_grid3[0, :, 0, :])
       
        
        mf = MFnLineArray(xdata = np.array([0, 0.1, 0.6, 1], dtype = 'f'),
                          ydata = np.array([0, 0.4, -0.5, 1], dtype = 'f'))
    
        bc_load2 = BCSlice(var = 'u', value = -0.05, dims = [1], time_function = mf.get_value,
                          slice = self.fe_grid3[-1, -1, -1, -1 ])
    
        load_dofs = self.fe_grid3[-1, -1, :, -1 ].dofs[:, :, 1].flatten()
        
        tstepper = TS(sdomain = self.fe_domain,
                       bcond_list = [bc_fixed,
                                     bc_link14,
                                     bc_link15,
                                     bc_link12,
                                     bc_link23,
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

    fbt = FoldedBondTest()

    fbt.tloop.eval()

    # Put the whole thing into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp(ibv_resource = fbt)
    app.main()

