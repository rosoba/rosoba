#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Jul 10, 2010 by: rch


from etsproxy.traits.api import \
    HasTraits, HasStrictTraits, Float, Property, cached_property, Instance, \
    Int, Button

from etsproxy.traits.ui.api import \
    View, Item, Group

import numpy as np
import os.path as path

from ibvpy.api import \
    IBVModel, TStepper, TLoop, BCSlice, \
    RTraceDomainListField

from ibvpy.mats.mats3D.mats3D_elastic.mats3D_elastic import \
    MATS3DElastic

from ibvpy.fets.fets_eval import \
    FETSEval
from ibvpy.fets.fets3D.fets3D8h import \
    FETS3D8H
from ibvpy.mesh.fe_grid import \
    FEGrid
from ibvpy.mesh.fe_refinement_grid import \
    FERefinementGrid
from ibvpy.mesh.fe_domain import \
    FEDomain

class Hypar(HasStrictTraits):

    lx = Float(10.0)
    ly = Float(2.0)
    lz = Float(0.1)

    z_base = Float(0.0)
    z_plus = Float(0.4)
    z_minus = Float(0.3)

    def __call__(self, X_arr):

        # derive the function

        xi, yi, zi = X_arr.T

        # number of local grid points for each coordinate direction
        # values must range between 0 and 1

        x = xi * self.lx
        y = yi * self.ly
        z = zi * self.lz
        a = self.lx ** 2 / self.z_minus
        b = self.ly ** 2 / self.z_plus
        z += y ** 2 / b - x ** 2 / a
        z += self.z_base

        return np.c_[ x, y, z ]

class HyparModel(IBVModel):
    '''SFB - Demontrator model specification.
    '''

    shape_x = Int(10, auto_set=False, enter_set=True, input=True)
    '''Number of elements in x direction
    '''

    shape_y = Int(4, auto_set=False, enter_set=True, input=True)
    '''Number of elements in x direction
    '''

    shape_z = Int(1, auto_set=False, enter_set=True, input=True)

    # dimensions of the shell structure
    length_x = Float(10.0, auto_set=False, enter_set=True, input=True)  # [m]
    length_y = Float(2.0, auto_set=False, enter_set=True, input=True)  # [m]
    low_thickness = Float(0.06, auto_set=False, enter_set=True, input=True)  # , ps_levels=(1.0, 2.0, 4)) # [m]
    mid_thickness = Float(0.38, auto_set=False, enter_set=True, input=True)  # , ps_levels=(1.0, 2.0, 4)) # [m]
    top_thickness = Float(0.06, auto_set=False, enter_set=True, input=True)  # , ps_levels=(1.0, 2.0, 4)) # [m]

    z_plus = Float(1.0, input=True)
    z_minus = Float(1.0, input=True)

    E_c = Float(28700, auto_set=False, enter_set=True,)  # [MPa]
    E_pur = Float(8700, auto_set=False, enter_set=True,)  # [MPa]
    nu = Float(0.2, auto_set=False, enter_set=True,)  # [-]

    mats_c = Instance(MATS3DElastic, input=True)
    def _mats_c_default(self):
        return MATS3DElastic(E=self.E_c, nu=self.nu)

    mats_pur = Instance(MATS3DElastic, input=True)
    def _mats_pur_default(self):
        return MATS3DElastic(E=self.E_pur, nu=self.nu)

    fets_c = Instance(FETSEval, transient=True, input=True)
    def _fets_c_default(self):
        return FETS3D8H(mats_eval=self.mats_c)

    fets_pur = Instance(FETSEval, transient=True, input=True)
    def _fets_pur_default(self):
        return FETS3D8H(mats_eval=self.mats_pur)

    fe_domain = Property(Instance(FEDomain), depends_on='+input')
    @cached_property
    def _get_fe_domain(self):

        low_base = 0.0
        mid_base = low_base + self.low_thickness
        top_base = mid_base + self.mid_thickness

        z_plus = self.z_plus
        z_minus = self.z_minus

        geo_hypar_low = Hypar(lx=self.length_x, ly=self.length_y,
                              lz=self.low_thickness, z_base=low_base,
                              z_minus=z_minus, z_plus=z_plus)
        geo_hypar_mid = Hypar(lx=self.length_x, ly=self.length_y,
                              lz=self.mid_thickness, z_base=mid_base,
                              z_minus=z_minus, z_plus=z_plus)
        geo_hypar_top = Hypar(lx=self.length_x, ly=self.length_y,
                              lz=self.top_thickness, z_base=top_base,
                              z_minus=z_minus, z_plus=z_plus)

        node_map = [[ -1., -1., -1.],
                   [  1., -1., -1.],
                   [  1., 1., -1.],
                   [ -1., 1., -1.],
                   [ -1., -1., 1.],
                   [  1., -1., 1.],
                   [  1., 1., 1.],
                   [ -1., 1., 1.],
                   ]
        vtk_cells = [[0, 1, 3, 2, 4, 5, 7, 6]]

        fets_c = FETS3D8H(mats_eval=MATS3DElastic(E=10),
                           geo_r=node_map,
                           vtk_r=node_map,
                           vtk_cells=vtk_cells,
                           )

        fets_c = self.fets_c
        fets_pur = self.fets_pur

        fe_domain = FEDomain()
        fe_rgrid_low = FERefinementGrid(domain=fe_domain, fets_eval=fets_c)
        fe_grid_low_shell = FEGrid(level=fe_rgrid_low,
                                     coord_min=(-1.0, -1.0, 0.0),
                                     coord_max=(1.0, 1.0, 1.0),
                                     geo_transform=geo_hypar_low,
                                     shape=(self.shape_x, self.shape_y, self.shape_z),
                                     fets_eval=fets_c)
        fe_rgrid_mid = FERefinementGrid(domain=fe_domain, fets_eval=fets_pur)
        fe_grid_mid_shell = FEGrid(level=fe_rgrid_mid,
                                     coord_min=(-1.0, -1.0, 0.0),
                                     coord_max=(1.0, 1.0, 1.0),
                                     geo_transform=geo_hypar_mid,
                                     shape=(self.shape_x, self.shape_y, self.shape_z),
                                     fets_eval=fets_pur)
        fe_rgrid_top = FERefinementGrid(domain=fe_domain, fets_eval=fets_c)
        fe_grid_top_shell = FEGrid(level=fe_rgrid_top,
                                     coord_min=(-1.0, -1.0, 0.0),
                                     coord_max=(1.0, 1.0, 1.0),
                                     geo_transform=geo_hypar_top,
                                     shape=(self.shape_x, self.shape_y, self.shape_z),
                                     fets_eval=fets_c)

        return fe_domain

    tloop = Property(depends_on='+input')
    @cached_property
    def _get_tloop(self):

        print self.fe_domain.n_dofs

        subdomains = self.fe_domain.subdomains

        fe_grid_low = subdomains[0].fe_subgrids[0]
        fe_grid_mid = subdomains[1].fe_subgrids[0]
        fe_grid_top = subdomains[2].fe_subgrids[0]

        low_grid_top_slice = fe_grid_low[:, :, -1, :, :, -1]
        mid_grid_low_slice = fe_grid_mid[:, :, 0, :, :, 0]
#
        mid_grid_top_slice = fe_grid_mid[:, :, -1, :, :, -1]
        top_grid_low_slice = fe_grid_top[:, :, 0, :, :, 0]

        left_support_slice = fe_grid_low[0, :, 0, 0, :, 0]
        right_support_slice = fe_grid_low[-1, :, 0, -1, :, 0]

        ts = TStepper(sdomain=self.fe_domain,
                      bcond_list=[
                                  BCSlice(var='u', slice=low_grid_top_slice,
                                       link_slice=mid_grid_low_slice, value=0.0,
                                       dims=[0, 1, 2], link_dims=[0, 1, 2],
                                       link_coeffs=[1., 1., 1.]),
                                  BCSlice(var='u', slice=mid_grid_top_slice,
                                          link_slice=top_grid_low_slice, value=0.0,
                                          dims=[0, 1, 2], link_dims=[0, 1, 2],
                                          link_coeffs=[1., 1., 1.]),
                                  BCSlice(var='u', slice=left_support_slice, value=0.0,
                                          dims=[1, 2]),
                                  BCSlice(var='u', slice=right_support_slice, value=0.0,
                                          dims=[1, 2]),
                                  BCSlice(var='u', slice=fe_grid_low[0, 0, 0, 0, 0, 0], value=0.0,
                                          dims=[0]),
                                  # LC1: dead load [MN/m^3]
                                  BCSlice(var='f', value= -0.0224 , dims=[2],
                                          integ_domain='global',
                                          slice=fe_grid_top[:, :, :, :, :, :]),
                                  BCSlice(var='f', value= -0.0224 , dims=[2],
                                          integ_domain='global',
                                          slice=fe_grid_low[:, :, :, :, :, :]),
                               ],
                      rtrace_list=[
                                   RTraceDomainListField(name='max principle stress' ,
                                                         idx=0,
                                                         var='max_principle_sig', warp=True,
                                                         record_on='update',),
                                   RTraceDomainListField(name='Epsilon 0' ,
                                                         var='eps0_app',
                                                         record_on='update'),
                                   RTraceDomainListField(name='Stress' ,
                                                         var='sig_app', idx=0, warp=False,
                                                         record_on='update',),
                                   ]
                      )
        tloop = TLoop(tstepper=ts, tolerance=1e-4,)
        return tloop

    generate_abaqus_input = Button()
    '''Generate input file for abaqus computation
    '''
    def _generate_abaqus_input_fired(self):
        self._generate_abaqus_input()

    def _generate_abaqus_input(self):
        subdomains = self.fe_domain.subdomains

        fe_grid_low = subdomains[0].fe_subgrids[0]
        fe_grid_mid = subdomains[1].fe_subgrids[0]
        fe_grid_top = subdomains[2].fe_subgrids[0]

        for idx, fe_grid in enumerate([fe_grid_low, fe_grid_mid, fe_grid_top]):
            geo_grid = fe_grid.geo_grid
            elem_node_map = geo_grid.cell_node_map + 1
            enum_arr = np.arange(len(elem_node_map)) + 1
            fn_nd = path.join(sim_model.dir, 'nodes_%d.dat' % idx)
            fn_el = path.join(sim_model.dir, 'elems_%d.dat' % idx)
            print 'saving data in %s' % sim_model.dir
            np.savetxt(fn_nd, geo_grid.cell_grid.point_X_arr, '%10.5f', delimiter=',')
            np.savetxt(fn_el, np.hstack([enum_arr[:, np.newaxis], elem_node_map]), '%d', delimiter=',')

    traits_view = View(Item('E_c'),
                       Item('E_pur'),
                       Item('length_x'),
                       Item('length_y'),
                       Item('low_thickness', label='Thickness of the lower shell'),
                       Item('mid_thickness', label='Thickness of the middle shell'),
                       Item('top_thickness', label='Thickness of the upper shell'),
                       Group(
                       Item('shape_x'),
                       Item('shape_y'),
                       Item('shape_z'),
                       label='number of elements'
                       ),
                       Item('z_minus'),
                       Item('z_plus'),
                       Item('generate_abaqus_input', show_label=False),
                       height=0.4,
                       width=0.4,
                       resizable=True)

if __name__ == '__main__':

    sim_model = HyparModel(shape_x=10,
                           shape_y=5,
                           shape_z=1,
                           z_plus=0.5,
                           z_minus=0.5,
                           length_x=8.0,
                           length_y=2.0,
                           low_thickness=0.06,
                           mid_thickness=0.20,
                           top_thickness=0.06,
                           )

    do = 'generate'

    if do == 'generate':
        sim_model._generate_abaqus_input()

    elif do == 'edit':
        sim_model.configure_traits()

    elif do == 'eval':

        sim_model.tloop.eval()
        from ibvpy.plugins.ibvpy_app import IBVPyApp
        app = IBVPyApp(ibv_resource=sim_model)
        app.main()
