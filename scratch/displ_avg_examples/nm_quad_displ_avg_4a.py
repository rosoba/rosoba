'''
Created on Jul 21, 2010

@author: jakub
'''
from ibvpy.api import \
    TStepper as TS, RTraceGraph, TLoop, \
    TLine, BCDof, IBVPSolve as IS, DOTSEval, FEDomain, FERefinementGrid, \
    FEGrid, BCSlice, RTraceDomainListField

import numpy as np

# from apps.scratch.jakub.mlab.mlab_trace import RTraceDomainListField
from ibvpy.mats.mats2D.mats2D_sdamage.mats2D_sdamage import MATS2DScalarDamage
from ibvpy.mats.mats2D.mats2D_elastic.mats2D_elastic import MATS2DElastic
from ibvpy.mats.mats2D.mats2D_cmdm.mats2D_cmdm import MATS2DMicroplaneDamage, PhiFnStrainSoftening
from ibvpy.mats.mats2D.mats2D_sdamage.strain_norm2d import Euclidean, Mazars, Rankine
from ibvpy.mats.mats_proxy import MATSProxy
from ibvpy.fets.fets2D.fets2D4q import FETS2D4Q
from ibvpy.fets.fets2D.fets2D4q9u import  FETS2D4Q9U
from ibvpy.fets.fets2D.fets2D4q8u import FETS2D4Q8U
from numpy import array, cos, sin, pi, sqrt, deg2rad, arctan
from mathkit.mfn import MFnLineArray
import matplotlib.pyplot as plt

from ibvpy.cntl.displ_avg.rt_nonlocal_averaging import \
    RTNonlocalAvg, QuarticAF

def app():
    avg_radius = 7
    thickness = 50.0  # mm
    md = MATS2DScalarDamage(E=30.0e3,
                            nu=0.2,
                            epsilon_0=1.0e-4,
                            epsilon_f=8.0e-4,
                            # epsilon_f = 12.0e-4, #test doubling the e_f
                            stress_state="plane_stress",
                            stiffness="secant",
                            # stiffness  = "algorithmic",
                            strain_norm=Rankine())

    mdm = MATS2DMicroplaneDamage(E=30.0e3,
                                 nu=0.2,
                                model_version='compliance',
                                phi_fn=PhiFnStrainSoftening(
                                                              G_f=0.1,
                                                              f_t=3.0,
                                                              md=0.0,
                                                              h=2. * avg_radius))

    # mdm = MATS2DElastic()

    fets_eval = FETS2D4Q9U(mats_eval=mdm)

    fe_domain = FEDomain()

    fe_rgrid = FERefinementGrid(name='fe_grid1', fets_eval=fets_eval, domain=fe_domain)

    n_half = 15
    n_el = n_half * 2 + 1
    # Discretization
    fe_grid = FEGrid(coord_max=(200.0, 200.0, 0.),
                      shape=(n_el, n_el),
                      fets_eval=fets_eval,
                      level=fe_rgrid)

    # numer of elements to deactivate from each side
    n_dact = 3 * n_half / 10
    print n_dact



    for i in range(0, n_dact):
        fe_grid.deactivate((i, n_half))

    for j in range(n_el - n_dact, n_el):
        fe_grid.deactivate((j, n_half))

    Ps_fn = MFnLineArray(xdata=array([0, 0.5, 1 ]),
                       ydata=array([0, 1., 1 ]))

    P_fn = MFnLineArray(xdata=array([0, 0.5, 1 ]),
                       ydata=array([0, 0., 1 ]))

    # averaging function
    avg_processor = RTNonlocalAvg(avg_fn=QuarticAF(radius=avg_radius,
                                                        correction=True))

    bc_slice_up = fe_grid[:, -1, :, -1]
    print 'BC Up', bc_slice_up.elems.flatten()

    bc_slice_left = fe_grid[0, n_half + 1:, 0, :]
    print 'BC Left', bc_slice_left.elems.flatten()

    bc_slice_right = fe_grid[-1, 0:n_half, -1, :]
    print 'BC Right', bc_slice_right.elems.flatten()

    bc_slice_down = fe_grid[:, 0, (1, 2), 0]
    print 'BC Down', bc_slice_down.elems.flatten()

    bc_loading = fe_grid[0, 0, 0, 0]

    loading_dof_ps = bc_loading.dofs.flatten()[0]
    loading_dof_p = bc_loading.dofs.flatten()[1]
    print 'loading_dof_ps', loading_dof_ps
    print 'loading_dof_p', loading_dof_p
    print 'dofs right', bc_slice_right.dofs.flatten()
    print 'dofs down', bc_slice_down.dofs.flatten()
    redundant_dofs_left = fe_grid[0:n_dact, n_half, :-1, 1]
    print 'redundant dofs left' , redundant_dofs_left.dofs.flatten()
    redundant_dofs_right = fe_grid[n_el - n_dact:n_el, n_half, 1:, 1]
    print 'redundant dofs right' , redundant_dofs_right.dofs.flatten()


    aa_x = np.hstack((loading_dof_ps, np.unique(bc_slice_down.dofs[:, :, 0].flatten()), np.unique(bc_slice_right.dofs[:, :, 0].flatten())))
    aa_y = np.hstack((loading_dof_p, np.unique(bc_slice_down.dofs[:, :, 1].flatten()), np.unique(bc_slice_right.dofs[:, :, 1].flatten())))

    print 'aa_x', aa_x
    print 'aa_y', aa_y

    link_right = BCSlice(var='u', value=0., dims=[0, 1],
                                   slice=bc_slice_right,
                                   link_slice=bc_loading,
                                   link_dims=[0, 1],
                                   link_coeffs=[1.])

    link_down = BCSlice(var='u', value=0., dims=[0, 1],
                                   slice=bc_slice_down,
                                   link_slice=bc_loading,
                                   link_dims=[0, 1],
                                   link_coeffs=[1.])

    ts = TS(sdomain=fe_grid,
            u_processor=avg_processor,
             bcond_list=[
                        # constraint for all left dofs in y-direction:
                        link_right, link_down,
                        # constraint for all left dofs in y-direction:
                        BCSlice(var='u', slice=bc_slice_up, dims=[0, 1], value=0.),
                        BCSlice(var='u', slice=bc_slice_left, dims=[0, 1], value=0.),
                        BCSlice(var='u', slice=redundant_dofs_left, dims=[0, 1], value=0.),
                        BCSlice(var='u', slice=redundant_dofs_right, dims=[0, 1], value=0.),
                        BCSlice(var='f', slice=bc_loading, dims=[0],
                                time_function=Ps_fn.get_value,
                                value= -100),
                        BCSlice(var='u', slice=bc_loading, dims=[1],
                                time_function=P_fn.get_value,
                                value= -0.05),

                        ],


             rtrace_list=[
                            RTraceGraph(name='Ps' ,
                                      var_y='F_int', idx_y_arr=aa_x,
                                      var_x='U_k', idx_x=loading_dof_ps,
                                      transform_x='-x', transform_y='-y*%g' % thickness,
                                      record_on='update'),
                            RTraceGraph(name='P' ,
                                      var_y='F_int', idx_y_arr=aa_y,
                                      var_x='U_k', idx_x=loading_dof_p,
                                      transform_x='-x', transform_y='-y*%g' % thickness,
                                      record_on='update'),
                            RTraceDomainListField(name='Deformation' ,
                                           var='eps_app', idx=0,
                                           record_on='update'),
                            RTraceDomainListField(name='Deformation' ,
                                           var='sig_app', idx=0,
                                           record_on='update'),
                            RTraceDomainListField(name='Displacement' ,
                                           var='u', idx=1,
                                           record_on='update',
                                           warp=True),
                            RTraceDomainListField(name='fracture_energy' ,
                                           var='fracture_energy', idx=0,
                                           record_on='update',
                                           warp=True),
                            RTraceDomainListField(name='Damage' ,
                                        var='omega_mtx', idx=0, warp=True,
                                        record_on='update'),
                        ]
            )

    # Add the time-loop control
    #
    tl = TLoop(tstepper=ts,
                tolerance=5.0e-5,
                KMAX=200,
                tline=TLine(min=0.0, step=.05, max=1.0))

    tl.eval()

    with open("4a.txt") as f:
        data = f.read()

    data = data.split('\n')

    x = [row.split(' ')[0] for row in data]
    y = [row.split(' ')[1] for row in data]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    print ts.rtrace_list
    P = ts.rtrace_list[1].trace
    ax.set_title("Load-displacement diagramm")
    ax.set_xlabel(r'$\delta$ (mm)')
    ax.set_ylabel('P (N)')
    ax.plot(x, y, 'r--', label='Test Load Path 4b')
    ax.plot(P.xdata, P.ydata, c='b', label='Simulation (MDM)')
    leg = ax.legend()
    plt.show()
#    # Put the whole stuff into the simulation-framework to map the
#    # individual pieces of definition into the user interface.
#    #
    from ibvpy.plugins.ibvpy_app import IBVPyApp
    ibvpy_app = IBVPyApp(ibv_resource=ts)
    ibvpy_app.main()


if __name__ == '__main__':
    app()
