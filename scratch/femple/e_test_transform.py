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

from ibvpy.fets.fets2D import FETS2D4Q8U, FETS2D4Q

#from ibvpy.fets.fets1D5 import FETS1D52L4ULRH

from e_09_fets_crack import FETS1D5t2L4ULRH

from ibvpy.mats.mats1D import MATS1DElastic

from ibvpy.mats.mats1D5.mats1D5_bond import MATS1D5Bond
    
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
    L1 = Float(1.0, desc = 'Length of the left part', input = True)
    L2 = Float(1.0, desc = 'Length of the right part', input = True)
    alpha = Float(0.0, desc = 'Fold angle', input = True)
    d = Float(1.0, desc = 'thickness of the plate', input = True)
    h = Float(0.1, desc = 'width of the bond', input = True)

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

    #-----------------
    # fets:
    #-----------------
    
    fets_bond = Property(Instance(FETSEval),
                          depends_on = '+input')
    @cached_property
    def _get_fets_bond (self):
        return FETS1D5t2L4ULRH(mats_eval = MATS1D5Bond(mats_phase1 = MATS1DElastic(E = 0),
                                                      mats_phase2 = MATS1DElastic(E = 0),
                                                      mats_ifslip = MATS1DElastic(E = 10.0),
                                                      mats_ifopen = MATS1DElastic(E = 1.0)),
                              alpha = self.alpha)
       
    fets_concrete = Property(Instance(FETSEval),
                                depends_on = '+input')
    @cached_property
    def _get_fets_concrete(self):
        return FETS2D4Q(mats_eval = MATS2DElastic(
                                                  E = 100,
                                                  nu = 0.0,
                                                  ),
                           alpha = self.alpha)    

    geo_r = Array(value = [[-1, -1], [1, -1], [1, 1], [-1, 1]], dtype = 'f')

    T_mtx = Property(depends_on = 'alpha')
    @cached_property
    def _get_T_mtx(self):
        return np.array([[ math.cos(self.alpha), math.sin(self.alpha)],
                         [ -math.sin(self.alpha), math.cos(self.alpha)]], dtype = 'f')
        
    def N_transform(self, r, X):
        cx = np.array(self.geo_r, dtype = 'float_')
        Nr = np.array([1 / 4. * (1 + r[:, 0] * cx[i, 0]) * (1 + r[:, 1] * cx[i, 1])
                      for i in range(0, 4) ])
        return np.dot(Nr.T, X)
    
    def gt_plate_left(self, points):
        X1 = np.array([[0, 0], [self.d, 0], [self.d, self.L1], [0, self.L1]], dtype = 'f')
        X1 = np.dot(X1, self.T_mtx)
        return self.N_transform(points, X1)
    
    def gt_bond(self, points):
        X1 = np.array([[0, self.L1],
                       [self.d, self.L1],
                       [self.d, self.L1 + self.h],
                       [0, self.L1 + self.h]],
                      dtype = 'f')
        X1 = np.dot(X1, self.T_mtx)
        return self.N_transform(points, X1)    
    
    def gt_plate_right(self, points):
        X1 = np.array([[0, self.L1 + self.h],
                       [self.d, self.L1 + self.h],
                       [self.d, self.L1 + self.h + self.L2],
                       [0, self.L1 + self.h + self.L2]], dtype = 'f')
        X1 = np.dot(X1, self.T_mtx)
        return self.N_transform(points, X1)
    
    fe_domain = Property(depends_on = '+input')
    @cached_property
    def _get_fe_domain(self):
        return FEDomain()
        
    fe_rg1 = Property(depends_on = '+input')
    @cached_property
    def _get_fe_rg1(self):
        return FERefinementGrid(name = 'rg1',
                                fets_eval = self.fets_concrete,
                                domain = self.fe_domain)
    
    fe_grid1 = Property(depends_on = '+input')
    @cached_property 
    def _get_fe_grid1(self): 
        return FEGrid(coord_min = (-1., -1.),
                          coord_max = (1., 1.),
                          shape = (self.n_x, self.n_z),
                          fets_eval = self.fets_concrete,
                          geo_transform = self.gt_plate_left,
                          level = self.fe_rg1)
    
    fe_rg2 = Property(depends_on = '+input')
    @cached_property
    def _get_fe_rg2 (self):
        return FERefinementGrid(name = 'rg2',
                                fets_eval = self.fets_bond,
                                domain = self.fe_domain)
    
    fe_grid2 = Property(depends_on = '+input')
    @cached_property
    def _get_fe_grid2(self):
        return FEGrid(coord_min = (-1., -1.),
                          coord_max = (1., 1.),
                          shape = (self.n_x, 1),
                          fets_eval = self.fets_bond,
                          geo_transform = self.gt_bond,
                          level = self.fe_rg2)
    
    fe_rg3 = Property(depends_on = '+input')
    @cached_property
    def _get_fe_rg3(self):
        return FERefinementGrid(name = 'rg3',
                                fets_eval = self.fets_concrete,
                                domain = self.fe_domain)
    
    fe_grid3 = Property(depends_on = '+input')
    @cached_property
    def _get_fe_grid3(self):
        return FEGrid(coord_min = (-1., -1.),
                      coord_max = (1., 1.),
                      shape = (self.n_x, self.n_z),
                      fets_eval = self.fets_concrete,
                      geo_transform = self.gt_plate_right,
                      level = self.fe_rg3)
    
    tloop = Property(depends_on = '+input')
    @cached_property
    def _get_tloop(self):
        self.fe_grid1
        self.fe_grid2
        self.fe_grid3
        print 'count dofs', self.fe_domain.n_dofs
       
        bc_fixed = BCSlice(var = 'u', value = 0., dims = [0, 1],
                           slice = self.fe_grid3[:, -1, :, -1])

        bc_link_plate_left_crack = BCDofGroup(var = 'u',
                                              value = 0.,
                                              dims = [0, 1],
                                              get_dof_method = self.fe_grid1.get_top_dofs,
                                              link_coeffs = [1.0], # , 1.0],
                                              link_dims = [0, 1],
                                              get_link_dof_method = self.fe_grid2.get_bottom_dofs
                                              )
    
        bc_link_crack_plate_right = BCDofGroup(var = 'u',
                                              value = 0.,
                                              dims = [0, 1],
                                              get_dof_method = self.fe_grid2.get_top_dofs,
                                              link_coeffs = [1.0], # , 1.0],
                                              link_dims = [0, 1],
                                              get_link_dof_method = self.fe_grid3.get_bottom_dofs
                                              )
    
        bc_load = BCSlice(var = 'u', value = -1.0, dims = [1], # time_function = mf.get_value,
                          integ_domain = 'local',
                          slice = self.fe_grid1[:, 0, :, 0 ])
    
        tstepper = TS(sdomain = self.fe_domain,
                       bcond_list = [bc_fixed,
                                     bc_link_plate_left_crack,
                                     bc_link_crack_plate_right,
                                     bc_load,
                                    ],
        rtrace_list = [
                         RTraceDomainListField(name = 'Displacement' ,
                                        var = 'u', idx = 0,
                                        record_on = 'update',
                                        warp = True),
                       RTraceGraph(name = 'F(u)' ,
                               var_y = 'F_int', idx_y = 0,
                               var_x = 'U_k', idx_x = 0),
                       RTraceGraph(name = 'F(v)' ,
                               var_y = 'F_int', idx_y = 1,
                               var_x = 'U_k', idx_x = 1),
                    ]
                )
    
        # Add the time-loop control
        tloop = TLoop(tstepper = tstepper, KMAX = 10, tolerance = 1e-4,
                      debug = False,
                       tline = TLine(min = 0.0, step = 1., max = 1.0))

        return tloop

if __name__ == '__main__':

    fbt = FoldedBondTest(
                         alpha = -math.pi / 2.0
                         )

    fbt.tloop.eval()

    # Put the whole thing into the simulation-framework to map the
    # individual pieces of definition into the user interface.
    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    app = IBVPyApp(ibv_resource = fbt)
    app.main()
