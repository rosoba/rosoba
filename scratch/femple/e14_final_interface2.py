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

from ibvpy.fets.fets1D5 import FETS1D52L4ULRH

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

    E_crack = Float(33e9, input = True)
    nu_crack = Float(0.2, input = True)
    Sigma_max_crack = Float(4e6, input = True)
    
    E_tex = Float(50e3, input = True)
    nu_tex = Float(0.2, input = True)
    
    E_platte = Float(33071625e3, input = True)
    nu_platte = Float(0.2, input = True)
            
    #-----------------
    # fets
    #-----------------
    
    fets_tex = Property(Instance(FETSEval),
                          depends_on = 'E_crack,nu_crack')
    @cached_property
    def _get_fets_tex(self):
        return FETS1D5t2L4ULRH(mats_eval = MATS1D5Bond(mats_phase1 = MATS1DElastic(E = 0.0),
                                                      mats_phase2 = MATS1DElastic(E = self.E_tex),
                                                      mats_ifopen = MATS1DElastic(E = 1e12),
                                                      mats_ifslip = MATS1DPlastic(E = 30000,
                                                                                  sigma_y = self.Sigma_max_crack,
                                                                                  K_bar = 0.,
                                                                                  H_bar = 0.)))
    fets_crack = Property(Instance(FETSEval),
                          depends_on = 'E_crack,nu_crack')
    @cached_property
    def _get_fets_crack(self):
        return FETS1D5t2L4ULRH(mats_eval = MATS1D5Bond(mats_phase1 = MATS1DElastic(E = self.E_crack),
                                                      mats_phase2 = MATS1DElastic(E = self.E_crack),
                                                      mats_ifslip = MATS1DElastic(E = self.E_crack),
                                                      mats_ifopen = MATS1DPlastic(E = 30000,
                                                          sigma_y = self.Sigma_max_crack,
                                                          K_bar = 0.,
                                                          H_bar = 0.)))
    
            
    fets_platte = Property(Instance(FETSEval),
                                 depends_on = 'E_concrete,nu_concrete')
    @cached_property
    def _get_fets_platte(self):
        return FETS2D4Q8U(mats_eval = MATS2DElastic(
                                                      E = self.E_platte,
                                                      nu = self.nu_platte,
                                                      ))    

    geo_r = Array(value = [[-1, -1], [1, -1], [1, 1], [-1, 1]], dtype = 'f')
        
    def N_transform(self, r, X):
        cx = np.array(self.geo_r, dtype = 'float_')
        Nr = np.array([1 / 4. * (1 + r[:, 0] * cx[i, 0]) * (1 + r[:, 1] * cx[i, 1])
                      for i in range(0, 4) ])
        return np.dot(Nr.T, X)
    ''' 
    def gt_crack_top(self, points):
        X1 = np.array([[0,18], [0,10],[2, 10], [2,18]], dtype = 'f')
        return self.N_transform(points, X1)
    '''
    def gt_crack(self, points):
        X1 = np.array([[0, 0], [0, 6], [-2, 6], [-2, 0], ], dtype = 'f')
        return self.N_transform(points, X1)

    def gt_tex(self, points):
        X1 = np.array([[-2, 3], [0, 3], [0, 4], [-2, 4]], dtype = 'f')
        return self.N_transform(points, X1)
    
    def gt_platte(self, points):
        X = np.array([[0, 0], [138, 0], [138, 6], [0, 6]], dtype = 'f')
        return self.N_transform(points, X)
            
    fe_domain = Property
    @cached_property
    def _get_fe_domain(self):
        return FEDomain()
        
    fe_rg_crack = Property
    @cached_property
    def _get_fe_rg_crack(self):
        return FERefinementGrid(name = 'crack top',
                                fets_eval = self.fets_crack,
                                domain = self.fe_domain)
    
    fe_grid_crack = Property
    @cached_property 
    def _get_fe_grid_crack(self): 
        return FEGrid(coord_min = (-1., -1.),
                          coord_max = (1., 1.),
                          shape = (2, 1),
                          fets_eval = self.fets_crack,
                          geo_transform = self.gt_crack,
                          level = self.fe_rg_crack)
    
    fe_rg_tex = Property
    @cached_property
    def _get_fe_rg_tex(self):
        return FERefinementGrid(name = 'tex',
                                fets_eval = self.fets_tex,
                                domain = self.fe_domain)
    
    fe_grid_tex = Property
    @cached_property 
    def _get_fe_grid_tex(self): 
        return FEGrid(coord_min = (-1., -1.),
                          coord_max = (1., 1.),
                          shape = (1, 1),
                          fets_eval = self.fets_tex,
                          geo_transform = self.gt_tex,
                          level = self.fe_rg_tex)
    
         
    fe_rg_platte = Property
    @cached_property
    def _get_fe_rg_platte (self):
        return FERefinementGrid(name = 'platte',
                                fets_eval = self.fets_platte,
                                domain = self.fe_domain)
    
    fe_grid_platte = Property
    @cached_property
    def _get_fe_grid_platte (self):
        return FEGrid(coord_min = (-1., -1.),
                          coord_max = (1., 1.),
                          shape = (1, 2),
                          fets_eval = self.fets_platte,
                          geo_transform = self.gt_platte,
                          level = self.fe_rg_platte)
    
    tloop = Property()
    @cached_property
    def _get_tloop(self):
        self.fe_grid_crack
        self.fe_grid_tex
        self.fe_grid_platte

        print self.fe_grid_crack[:, -1, :, -1].dofs
        
        bc_fixed_crack = BCSlice(var = 'u', value = 0., dims = [0, 1],
#                                    get_dof_method = self.fe_grid_crack.get_top_dofs)
                                     slice = self.fe_grid_crack[ :, -1, :, -1])
        
        bc_link_platte_crack = BCSlice(var = 'u',
                                      value = 0.,
                                      dims = [0, 1],
                                      #get_dof_method = self.fe_grid_crack.get_bottom_dofs,
                                      slice = self.fe_grid_crack[:, 0, :, 0],
                                      link_coeffs = [1.0, 1.0],
                                      link_dims = [0, 1],
                                      #get_link_dof_method = self.fe_grid_platte.get_left_dofs,
                                      link_slice = self.fe_grid_platte[0, :, 0, :]
                                      )
                   
        bc_link_tex_crack = BCSlice(var = 'u',
                                      value = 0.,
                                      dims = [0, 1],
                                      slice = self.fe_grid_tex[:, 0, :, 0],
                                      link_coeffs = [1.0, 1.0],
                                      link_dims = [0, 1],
                                      link_slice = self.fe_grid_crack[0, :, -1, :]
                                      )
                   
        mf = MFnLineArray(xdata = np.array([0, 0.1, 0.6, 1], dtype = 'f'),
                          ydata = np.array([0, 0.4, -0.5, 1], dtype = 'f'))
    
        bc_load_top = BCSlice(var = 'u', value = 1, dims = [0], time_function = mf.get_value,
                          slice = self.fe_grid_platte[-1, -1, -1, -1 ])
        bc_load_bottom = BCSlice(var = 'u', value = -1, dims = [0], time_function = mf.get_value,
                          slice = self.fe_grid_platte[-1, 0, -1, 0 ])
    
        load_dofs = self.fe_grid_platte[-1, -1, :, -1 ].dofs[:, :, 1].flatten()
        
        tstepper = TS(sdomain = self.fe_domain,
                       bcond_list = [bc_fixed_crack,
                                     bc_link_tex_crack,
                                     bc_link_platte_crack,
                                     bc_load_top,
                                     bc_load_bottom,
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

    fbt = FoldedBondTest(n_x = 5, n_z = 3)

    fbt.tloop.eval()

    # Put the whole thing into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp(ibv_resource = fbt)
    app.main()

