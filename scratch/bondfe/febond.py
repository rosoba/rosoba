'''
Created on Sep 12, 2015

@author: rch
'''

from ibvpy.fets.fets2D import FETS2D4Q
from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
from ibvpy.mesh.fe_grid import FEGrid
from mathkit.matrix_la.sys_mtx_assembly import SysMtxAssembly
import numpy as np


if __name__ == '__main__':

    #=========================================================================
    # Material matrix
    #=========================================================================
    mats_eval = MATS2DElastic(E=28e3, nu=0.2)
    D_el = mats_eval.D_el
    n_s = D_el.shape[0]

    #=========================================================================
    # Element definition
    #=========================================================================
    fets_eval = FETS2D4Q(mats_eval=mats_eval)
    n_geo_r, n_dim_geo = fets_eval.geo_r.shape
    n_dof_r, n_dim_dof = fets_eval.dof_r.shape
    n_ip = fets_eval.n_gp
    n_el_dofs = n_dof_r * n_dim_dof

    r_ip = fets_eval.ip_coords[:, :-1].T
    # [ i ]
    w_ip = fets_eval.ip_weights
    # [ d, n ]
    geo_r = fets_eval.geo_r.T
    # [ d, n, i ]
    dNr_geo = geo_r[
        :, :, None] * (1 + np.flipud(r_ip)[:, None, :] * np.flipud(geo_r)[:, :, None]) / 4.0
    # [ i, n, d ]
    dNr_geo = np.einsum('dni->ind', dNr_geo)

    #=========================================================================
    # Discretization
    #=========================================================================

    # Number of elements
    n_e_x = 2
    n_e_y = 1
    # length
    L_x = 20.0
    L_y = 10.0
    # [ r, i ]
    domain = FEGrid(coord_max=(L_x, L_y),
                    shape=(n_e_x, n_e_y),
                    fets_eval=fets_eval)
    n_e = domain.n_active_elems
    n_dofs = domain.n_dofs
    # element array with nodal coordinates
    # [ n_e, n_geo_r, n_dim_geo ]
    elem_x_map = domain.elem_X_map
    # [ n_e, n_dof_r, n_dim_dof ]
    elem_dof_map = domain.elem_dof_map

    # [ n_e, n_ip, n_geo_r, n_dim_geo ]
    J_mtx = np.einsum('ind,enf->eidf', dNr_geo, elem_x_map)
    J_inv = np.linalg.inv(J_mtx)
    J_det = np.linalg.det(J_mtx)

    # shape function for the unknowns are identical with the geomeetric
    # approximation
    dNr = dNr_geo

    # [ n_e, n_ip, n_dof_r, n_dim_dof ]
    dNx = np.einsum('eidf,inf->eind', J_inv, dNr)

    B = np.zeros((n_e, n_ip, n_dof_r, n_s, n_dim_dof), dtype='f')
    B_n_rows, B_n_cols, dNx_d = [0, 1, 2, 2], [0, 1, 0, 1], [0, 1, 1, 0]
    B[:, :, :, B_n_rows, B_n_cols] = dNx[:, :, :, dNx_d]

    #=========================================================================
    # System matrix
    #=========================================================================
    K = np.einsum('i,einsd,st,eimtf,ei->endmf', w_ip, B, D_el, B, J_det)
    K_mtx = SysMtxAssembly()
    K_mtx.add_mtx_array(K.reshape(-1, n_el_dofs, n_el_dofs), elem_dof_map)

    #=========================================================================
    # Load vector
    #=========================================================================

    R = np.zeros((n_dofs,), dtype='float_')
    R[[8, 10]] = 1.0
    K_mtx.register_constraint(a=0)
    K_mtx.register_constraint(a=1)
    K_mtx.register_constraint(a=2)
    u = K_mtx.solve(R)
    print 'u', u
