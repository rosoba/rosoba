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
    HasTraits, Float, Property, cached_property, Instance, \
    Int

import numpy as np
import os.path as path

from ibvpy.api import \
    IBVModel

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

class Hypar(HasTraits):

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

    n_dofs_xy = Int(16, ps_levels=(10, 30, 3))
    '''Number of elements in all dims
    '''

    shape_xy = Property
    '''Number of
    '''
    def _get_shape_xy(self):
        if self.fets == self.fe_linear:
            return self.n_dofs_xy - 1
        elif self.fets == self.fe_quad_serendipity \
            or self.fets == self.fe2d5_quad_serendipity \
            or self.fets == self.fe_quad_lagrange:
            return int(self.n_dofs_xy / 2)
        else:
            raise ValueError

    shape_z = Int(1, ps_levels=(1, 3, 1))

    # dimensions of the shell structure
    length_x = Float(5.0, ps_levels=(4.0, 5.0, 3)) # [m]
    length_y = Float(1.0, ps_levels=(4.0, 5.0, 3)) # [m]
    low_thickness = Float(0.06, ps_levels=(1.0, 2.0, 4)) # [m]
    mid_thickness = Float(0.18, ps_levels=(1.0, 2.0, 4)) # [m]
    top_thickness = Float(0.06, ps_levels=(1.0, 2.0, 4)) # [m]

    E = Float(28700) # [MPa]
    nu = Float(0.2) # [-]

    fets = Instance(FETSEval)
    def _fets_default(self):
        return self.fe_linear

    mats = Instance(MATS3DElastic)
    def _mats_default(self):
        return MATS3DElastic(E=self.E, nu=self.nu)

    fe_linear = Instance(FETSEval, transient=True)
    def _fe_linear_default(self):
        return FETS3D8H(mats_eval=self.mats)

    fe_domain = Property(Instance(FEDomain))
    @cached_property
    def _get_fe_domain(self):

        low_base = 0.0
        mid_base = low_base + self.low_thickness
        top_base = mid_base + self.mid_thickness

        geo_hypar_low = Hypar(lx=self.length_x, ly=self.length_y, lz=self.low_thickness, z_base=low_base)
        geo_hypar_mid = Hypar(lx=self.length_x, ly=self.length_y, lz=self.mid_thickness, z_base=mid_base)
        geo_hypar_top = Hypar(lx=self.length_x, ly=self.length_y, lz=self.top_thickness, z_base=top_base)

        fets_4u = FETS3D8H(mats_eval=MATS3DElastic(E=10),
                           geo_r=np.array([[ -1., -1., -1.],
                                           [  1., -1., -1.],
                                           [  1., 1., -1.],
                                           [ -1., 1., -1.],
                                           [ -1., -1., 1.],
                                           [  1., -1., 1.],
                                           [  1., 1., 1.],
                                           [ -1., 1., 1.],
                                           ], dtype='f')
                           )

        fe_domain = FEDomain()
        fe_rgrid = FERefinementGrid(domain=fe_domain, fets_eval=fets_4u)
        fe_grid_low_shell = FEGrid(level=fe_rgrid,
                                     coord_min=(-1.0, -1.0, -1.0),
                                     coord_max=(1.0, 1.0, 1.0),
                                     geo_transform=geo_hypar_low,
                                     shape=(self.shape_xy, self.shape_xy, self.shape_z),
                                     fets_eval=fets_4u)
        fe_grid_mid_shell = FEGrid(level=fe_rgrid,
                                     coord_min=(-1.0, -1.0, -1.0),
                                     coord_max=(1.0, 1.0, 1.0),
                                     geo_transform=geo_hypar_mid,
                                     shape=(self.shape_xy, self.shape_xy, self.shape_z),
                                     fets_eval=fets_4u)
        fe_grid_top_shell = FEGrid(level=fe_rgrid,
                                     coord_min=(-1.0, -1.0, -1.0),
                                     coord_max=(1.0, 1.0, 1.0),
                                     geo_transform=geo_hypar_top,
                                     shape=(self.shape_xy, self.shape_xy, self.shape_z),
                                     fets_eval=fets_4u)

        return fe_domain

if __name__ == '__main__':

    sim_model = HyparModel(n_dofs_xy=10,
                           shape_z=1,
                           length_x=10.0,
                           length_y=3.0,
                           low_sh_thickness=0.06,
                           mid_sh_thickness=0.18,
                           top_sh_thickness=0.06,
                           )

    do = 'generate_mesh'

    if do == 'generate_mesh':
        subdomains = sim_model.fe_domain.subdomains
        fe_grid = subdomains[0].fe_subgrids[0]
        geo_grid = fe_grid.geo_grid

        fn_nd = path.join(sim_model.dir, 'nodes.dat')
        fn_el = path.join(sim_model.dir, 'elems.dat')

        print 'saving data in %s', sim_model.dir
        np.savetxt(fn_nd, geo_grid.cell_grid.point_X_arr, '%10.5f')
        np.savetxt(fn_el, geo_grid.cell_node_map + 1, '%d')
