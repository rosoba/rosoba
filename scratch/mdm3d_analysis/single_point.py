'''
Created on Oct 2, 2013

@author: rch
'''

from ibvpy.mats.mats2D.mats2D_explorer_bcond import BCDofProportional
from ibvpy.api import RTraceGraph, TLine, BCDof
from ibvpy.mats import \
    MATS3DMicroplaneDamage, PhiFnStrainSoftening
from ibvpy.mats.mats_explore import MATSExplore
from ibvpy.mats.mats3D.mats3D_explore import MATS3DExplore
import math

if __name__ == '__main__':

    ec = {
          # overload the default configuration
#          'bcond_list'  : [ BCDofProportional(max_strain=0.0002, alpha_rad=math.pi / 6.0) ],
          'bcond_list'  : [ BCDof(var='u', dof=5, value=0.0002) ],
          'rtrace_list' : [
               RTraceGraph(name='stress - strain',
                           var_x='eps_app', idx_x=0,
                           var_y='sig_app', idx_y=0,
                           update_on='update'),
                        ],
          }

    mats_eval = MATS3DMicroplaneDamage(
                                    regularization=False,
                                    # n_mp=30,
                                    nu=0.2,
                                    E=30e+3,
                                    elastic_debug=False,
                                    # stress_state='plane_stress',
                                    symmetrization='product-type',
                                    model_version='compliance',
                                    phi_fn=PhiFnStrainSoftening(G_f=0.001117,
                                                                f_t=2.8968),
                                    )

    print 'D2_e', mats_eval.D2_e

    mats_eval.regularization = False

    mats_explore2d = MATS3DExplore(
                                 mats_eval=mats_eval,
                                 explorer_config=ec,
                                 )

    print 'ssss', mats_explore2d.mats_eval.regularization

    me = MATSExplore(dim=mats_explore2d)
    me.tloop.tline = TLine(min=0.0, step=1.0 / 1.0, max=1.0)
    u = me.tloop.eval()
    me.tloop.rtrace_list[0].configure_traits()
